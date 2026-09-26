#!/usr/bin/env python3
"""Bounds-aware model topology and strict gene-grounded nutrient-route audit.

Never reads sample abundance, ANCOM results or case/control labels. Genes are
verified against the input FASTA after CarveMe's SBML identifier sanitisation.
Non-seed/export-route status is NOT proof of secretion or feasible synthesis.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, re, sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET
import networkx as nx

NS = {'s': 'http://www.sbml.org/sbml/level3/version1/core',
      'f': 'http://www.sbml.org/sbml/level3/version1/fbc/version2'}
F = '{' + NS['f'] + '}'
METHOD = 'bounds_aware_compartment_graph_source_SCC_max5_v1'
PHYLOMINT_COMMIT = 'a394c6dd8ac53a90994f91b6219a0dd1b27c44f5'
WHITELIST_CANDIDATES = {
 'organic_acids': ['ac','but','ppa','lac__D','lac__L','succ','pyr','for','fum','mal__L','cit','akg','bhb','h4but'],
 'amino_acids': ['ala__L','arg__L','asn__L','asp__L','cys__L','gln__L','glu__L','gly','his__L','ile__L','leu__L','lys__L','met__L','phe__L','pro__L','ser__L','thr__L','trp__L','tyr__L','val__L','orn','4abut'],
 'B_vitamins': ['thm','ribflv','nac','ncam','pnto__R','pydx','pydam','pydxn','btn','fol','cbl1','adocbl']}

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def save_json(path, obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def write_tsv(path, rows, fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(rows[0])
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'wt',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)

def read_model(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rb') as f: root=ET.parse(f).getroot()
    mets={}
    for m in root.findall('.//s:species',NS):
        ident=m.get('id');comp=m.get('compartment','').removeprefix('C_')
        bare=ident.removeprefix('M_')
        if not bare.endswith('_'+comp): raise ValueError(f'Unmapped metabolite compartment: {ident}/{comp}')
        mets[ident]={'id':ident,'base_id':bare[:-(len(comp)+1)],'name':m.get('name',''),
                    'compartment':comp,'formula':m.get(F+'chemicalFormula',''),
                    'annotations':sorted({v for e in m.iter() for k,v in e.attrib.items() if k.endswith('resource')})}
    params={e.get('id'):float(e.get('value')) for e in root.findall('.//s:parameter',NS)}
    objectives={e.get(F+'reaction') for e in root.findall('.//f:fluxObjective',NS)}
    reactions=[]
    for r in root.findall('.//s:reaction',NS):
        sub=[e.get('species') for e in r.findall('s:listOfReactants/s:speciesReference',NS)]
        prod=[e.get('species') for e in r.findall('s:listOfProducts/s:speciesReference',NS)]
        lower=params.get(r.get(F+'lowerFluxBound'));upper=params.get(r.get(F+'upperFluxBound'))
        if lower is None or upper is None: raise ValueError(f'Missing FBC2 bounds {r.get("id")}')
        genes=sorted({e.get(F+'geneProduct') for e in r.findall('.//f:geneProductRef',NS)})
        reactions.append({'id':r.get('id'),'name':r.get('name',''),'substrates':sub,'products':prod,
                          'lower':lower,'upper':upper,'genes':genes,'objective':r.get('id') in objectives})
    return mets,reactions

def source_seeds(graph,max_size=5):
    weights={};components=[]
    for component in nx.strongly_connected_components(graph):
        if len(component)>max_size:continue
        # Match PhyloMint's singleton rule exactly (self-loop singleton not seed).
        if len(component)==1:
            node=next(iter(component))
            if graph.in_degree(node)!=0:continue
        elif any(a not in component for node in component for a,_ in graph.in_edges(node)):continue
        components.append(sorted(component))
        for node in component:weights[node]=1/len(component)
    return dict(sorted(weights.items())),sorted(components)

def exclusion(r):
    rid=r['id'].removeprefix('R_')
    if r['objective'] or 'biomass' in (rid+' '+r['name']).lower():return 'biomass_or_objective'
    if rid.startswith(('DM_','SK_','sink_','demand_')):return 'demand_or_sink'
    if not r['substrates'] or not r['products']:return 'boundary_exchange_or_pseudoreaction'
    if r['lower']==0 and r['upper']==0:return 'blocked_by_bounds'
    return None

def build_graph(mets,reactions,stock=False):
    graph=nx.DiGraph();edges=[];excluded=[]
    for r in reactions:
        reason=None if stock else exclusion(r)
        if reason:excluded.append({'reaction_id':r['id'],'reason':reason});continue
        dirs=[(r['substrates'],r['products'],'forward')] if stock else []
        if not stock and r['upper']>0:dirs.append((r['substrates'],r['products'],'forward'))
        if not stock and r['lower']<0:dirs.append((r['products'],r['substrates'],'reverse'))
        for sub,prod,direction in dirs:
            for a in sub:
                for b in prod:
                    graph.add_edge(a,b)
                    edges.append({'source':a,'target':b,'reaction_id':r['id'],'direction':direction})
    return graph,edges,excluded

def input_gene_ids(proteins):
    ids=[line[1:].split()[0] for line in Path(proteins).read_text().splitlines() if line.startswith('>')]
    mapped=['G_'+re.sub('[^A-Za-z0-9_]','_',i) for i in ids]
    if len(set(mapped))!=len(mapped):raise ValueError('Protein identifiers collide after CarveMe sanitisation')
    return set(mapped),len(ids)

def route_evidence(mets,reactions,known_genes,whitelist):
    # Each edge conserves a validated compound ID across different compartments;
    # biochemical conversions and generic EX reactions never create transport.
    strict=defaultdict(nx.DiGraph);loose=defaultdict(nx.DiGraph)
    evidence=[];boundary=defaultdict(lambda:{'in':[],'out':[]})
    transport_ids=set();supported_ids=set();exchange_ids=set()
    for r in reactions:
        subs,prods=r['substrates'],r['products'];participants=set(subs+prods)
        comps={mets[x]['compartment'] for x in participants}
        if len(participants)==1 and (not subs or not prods):
            m=next(iter(participants))
            if mets[m]['compartment']=='e':
                exchange_ids.add(r['id']);base=mets[m]['base_id']
                if (bool(subs) and r['lower']<0) or (bool(prods) and r['upper']>0):boundary[base]['in'].append(r['id'])
                if (bool(subs) and r['upper']>0) or (bool(prods) and r['lower']<0):boundary[base]['out'].append(r['id'])
        if len(comps)<2 or exclusion(r):continue
        transport_ids.add(r['id'])
        grounded=bool(r['genes']) and all(g in known_genes for g in r['genes'])
        if grounded:supported_ids.add(r['id'])
        directions=[]
        if r['upper']>0:directions.append((subs,prods))
        if r['lower']<0:directions.append((prods,subs))
        for src,dst in directions:
            for a in src:
                for b in dst:
                    ma,mb=mets[a],mets[b];base=ma['base_id']
                    if base!=mb['base_id'] or ma['compartment']==mb['compartment'] or base not in whitelist:continue
                    loose[base].add_edge(a,b,reaction=r['id'])
                    if grounded:strict[base].add_edge(a,b,reaction=r['id'])
                    evidence.append({'metabolite_id':base,'source':a,'target':b,'reaction_id':r['id'],
                                     'gene_supported':grounded,'genes':';'.join(r['genes']),
                                     'missing_input_genes':';'.join(g for g in r['genes'] if g not in known_genes)})
    routes={}
    for base in sorted(whitelist):
        variants={m['compartment']:mid for mid,m in mets.items() if m['base_id']==base}
        record={'inward':False,'outward':False,'unsupported_inward':False,'unsupported_outward':False,
                'inward_route':[],'outward_route':[],'exchange_in':boundary[base]['in'],'exchange_out':boundary[base]['out']}
        for direction,start,end,key in [('inward','e','c','in'),('outward','c','e','out')]:
            if start not in variants or end not in variants or not boundary[base][key]:continue
            a,b=variants[start],variants[end]
            for name,graphs in [('strict',strict),('loose',loose)]:
                g=graphs[base]
                if a in g and b in g and nx.has_path(g,a,b):
                    if name=='strict':
                        record[direction]=True;path=nx.shortest_path(g,a,b)
                        record[direction+'_route']=[{'source':u,'target':v,'reaction_id':g[u][v]['reaction']} for u,v in zip(path,path[1:])]
                    else:record['unsupported_'+direction]=True
        routes[base]=record
    return routes,evidence,transport_ids,supported_ids,exchange_ids

def freeze_whitelist(reference,outdir):
    path=outdir/'metabolite_whitelist.tsv';lockpath=outdir/'extraction_lock.json'
    if path.exists():
        lock=json.loads(lockpath.read_text())
        if sha256(path)!=lock['whitelist_sha256'] or sha256(reference)!=lock['reference_sha256']:raise ValueError('Frozen whitelist/reference hash drift')
        return list(csv.DictReader(path.open(),delimiter='\t')),lock
    mets,_=read_model(reference);rows=[];excluded=[];order=0
    for category,ids in WHITELIST_CANDIDATES.items():
        for ident in ids:
            variants=[m for m in mets.values() if m['base_id']==ident]
            if not variants:
                excluded.append({'metabolite_id':ident,'reason':'absent_from_pinned_bacterial_universe_not_assumed_biological_absence'});continue
            names=sorted({m['name'] for m in variants});formulas=sorted({m['formula'] for m in variants})
            if len(names)!=1 or len(formulas)!=1:raise ValueError(f'Conflicting compartment identity {ident}: {names}/{formulas}')
            order+=1;rows.append({'metabolite_id':ident,'metabolite_name':names[0],'category':category,
                                 'validated':'TRUE','display_order':order,'formula':formulas[0],
                                 'compartments':';'.join(sorted({m['compartment'] for m in variants})),
                                 'annotation_uris':';'.join(sorted({a for m in variants for a in m['annotations']}))})
    write_tsv(path,rows);write_tsv(outdir/'excluded_whitelist_candidates.tsv',excluded,['metabolite_id','reason'])
    lock={'locked_utc':datetime.now(timezone.utc).isoformat(),'method':METHOD,'max_scc_size':5,
          'phyloMint_commit':PHYLOMINT_COMMIT,'whitelist_sha256':sha256(path),'reference_sha256':sha256(reference),
          'selection_used_case_control_or_ANCOM':False,'strict_transport':'same compound e->p->c or e->c path; every step nonempty GPR with all IDs in input proteins; compatible default exchange bound',
          'seeds':'source SCC <=5, confidence 1/size; singleton must have indegree0; infer before whitelist; retain compartment identity',
          'primary_graph':'all substrate-product edges in allowed bound directions; excludes biomass/objective, demand/sink/boundary, bound-blocked; no currency-metabolite pruning',
          'base_join':'validated literal BiGG ID; keep stereoisomers; strict provider requires cytosolic graph node and no seed variant of same base',
          'roles':'all compartment-node seeds, not restricted whitelist; matrix row focal recipient/initiator A; column partner B',
          'claim':'predicted graph requirement/nonseed potential with evidenced route, never measured or FBA-verified secretion'}
    save_json(lockpath,lock)
    return rows,lock

def extract_one(row,whitelist,outdir,lock):
    cid=row['catalog_id'];path=Path(row['model_local_path']);proteins=Path(row['genes_local_path'])
    complete=json.loads(Path(row['complete_local_path']).read_text());qc=json.loads(Path(row['qc_local_path']).read_text())
    for key in ('model.xml','genes.faa','qc.json'):
        expected=complete['output_sha256'].get(key)
        if expected is None or sha256(path.parent/key)!=expected:raise ValueError(f'Completion hash mismatch {cid}/{key}')
    if not qc.get('technical_qc_pass'):raise ValueError('Technical reconstruction QC failed')
    dest=outdir/'models'/cid
    cache_identity={'model_sha256':sha256(path),'proteins_sha256':sha256(proteins),
                    'qc_sha256':sha256(Path(row['qc_local_path'])),'whitelist_sha256':lock['whitelist_sha256'],
                    'extractor_sha256':sha256(__file__),'reference_sha256':lock['reference_sha256']}
    cache=dest/'cache_complete.json'
    if cache.is_file():
        previous=json.loads(cache.read_text())
        if previous.get('input_identity')==cache_identity and all((dest/name).is_file() and sha256(dest/name)==value for name,value in previous.get('output_sha256',{}).items()) and previous.get('output_sha256'):
            payload=json.loads((dest/'graph.json').read_text())
            return payload,previous['qc_metrics']
    mets,rxns=read_model(path);known,ngenes=input_gene_ids(proteins)
    graph,edges,excluded=build_graph(mets,rxns);weights,components=source_seeds(graph)
    stock,_,_=build_graph(mets,rxns,stock=True);stock_weights,_=source_seeds(stock)
    routes,evidence,trans,supported,exchanges=route_evidence(mets,rxns,known,whitelist)
    base_seeds={mets[x]['base_id'] for x in weights};base_nodes={mets[x]['base_id'] for x in graph}
    cytosol={m['base_id'] for i,m in mets.items() if i in graph and m['compartment']=='c'}
    topo_req=sorted(base_seeds & whitelist);topo_prov=sorted((base_nodes-base_seeds)&whitelist)
    req=[x for x in topo_req if routes[x]['inward']]
    prov=[x for x in topo_prov if x in cytosol and routes[x]['outward']]
    stock_bases={mets[x]['base_id'] for x in stock}
    stock_seed_bases={mets[x]['base_id'] for x in stock_weights}
    stock_cytosol={mets[x]['base_id'] for x in stock if mets[x]['compartment']=='c'}
    stock_req=sorted(x for x in stock_seed_bases&whitelist if routes[x]['inward'])
    stock_prov=sorted(x for x in (stock_bases-stock_seed_bases)&whitelist&stock_cytosol if routes[x]['outward'])
    payload={'catalog_id':cid,'schema_version':1,'graph_method':METHOD,'requirements':req,'providers':prov,
             'topology_only_requirements':topo_req,'topology_only_providers':topo_prov,
             'stock_requirements':stock_req,'stock_providers':stock_prov,
             'seed_weights':weights,'seeds':sorted(weights),'nonseeds':sorted(set(graph)-set(weights)),
             'source_components':components,'metabolites':mets,'graph_nodes':sorted(graph),
             'strict_route_evidence':routes,'stock_phyloMint_seed_weights':stock_weights,
             'stock_phyloMint_graph_nodes':sorted(stock),'whitelist_sha256':lock['whitelist_sha256'],
             'model_sha256':sha256(path),'proteins_sha256':sha256(proteins),'input_protein_count':ngenes}
    dest.mkdir(parents=True,exist_ok=True)
    save_json(dest/'graph.json',payload)
    write_tsv(dest/'edges.tsv.gz',edges,['source','target','reaction_id','direction'])
    write_tsv(dest/'excluded_reactions.tsv',excluded,['reaction_id','reason'])
    write_tsv(dest/'transport_evidence.tsv',evidence,['metabolite_id','source','target','reaction_id','gene_supported','genes','missing_input_genes'])
    metrics={'catalog_id':cid,'genes':ngenes,'model_genes':qc['genes'],'reactions':len(rxns),'metabolites':len(mets),
             'exchange_reactions':len(exchanges),'transport_reactions':len(trans),'gene_supported_transport_reactions':len(supported),
             'exchange_fraction':len(exchanges)/len(rxns),'transport_fraction':len(trans)/len(rxns),
             'gene_supported_transport_fraction':len(supported)/len(rxns),
             'n_requirements':len(req),'n_providers':len(prov),'n_topology_only_requirements':len(topo_req),
             'n_stock_requirements':len(stock_req),'n_stock_providers':len(stock_prov),
             'n_seeds':len(weights),'n_stock_seeds':len(stock_weights),'n_graph_nodes':len(graph),
             'technical_qc_pass':True,'graph_status':'complete','graph_method':METHOD}
    save_json(cache,{'input_identity':cache_identity,'qc_metrics':metrics,
                    'output_sha256':{name:sha256(dest/name) for name in ['graph.json','edges.tsv.gz','excluded_reactions.tsv','transport_evidence.tsv']}})
    return payload,metrics

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--inventory',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True);ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--stage',choices=['pilot','partial','full'],default='pilot');ap.add_argument('--freeze-only',action='store_true')
    args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    white,lock=freeze_whitelist(args.reference,args.outdir)
    if args.freeze_only:print(json.dumps(lock,indent=2));return
    rows=list(csv.DictReader(args.inventory.open(),delimiter='\t'))
    if args.stage=='pilot':rows=[r for r in rows if r['model_source']=='pilot_v1']
    planned=len(rows);pending=[]
    if args.stage=='partial':
        available=[]
        for row in rows:
            if all(Path(row[key]).is_file() for key in ['model_local_path','genes_local_path','qc_local_path','complete_local_path']):available.append(row)
            else:pending.append(row['catalog_id'])
        rows=available
    result={};metrics=[];failed=[]
    for row in rows:
        cid=row['catalog_id']
        try:
            data,qc=extract_one(row,{w['metabolite_id'] for w in white},args.outdir,lock)
            result[cid]={k:data[k] for k in ['catalog_id','requirements','providers','topology_only_requirements','topology_only_providers','stock_requirements','stock_providers','model_sha256','whitelist_sha256','graph_method']}
            metrics.append(qc);print(f'{cid}\tU={len(data["requirements"])}\tP={len(data["providers"])}',flush=True)
        except Exception as error:
            failed.append({'catalog_id':cid,'error':repr(error)});print(f'FAILED {cid}: {error}',file=sys.stderr,flush=True)
    suffix='.partial' if args.stage=='partial' else ''
    if metrics:write_tsv(args.outdir/('model_graph_qc'+suffix+'.tsv'),metrics)
    save_json(args.outdir/('model_support_sets'+suffix+'.json'),result)
    save_json(args.outdir/('extraction_audit'+suffix+'.json'),{'stage':args.stage,'planned_models':planned,'completed_models':len(result),
              'failed_models':failed,'pending_local_models':pending,'complete':len(result)==planned,'nonempty_strict_requirements':sum(bool(v['requirements']) for v in result.values()),
              'nonempty_strict_providers':sum(bool(v['providers']) for v in result.values()),
              'whitelist_sha256':lock['whitelist_sha256'],'extractor_sha256':sha256(__file__),
              'networkx_version':nx.__version__,'graph_method':METHOD,'no_donor_outcomes_read':True})
    if failed:sys.exit(2)

if __name__=='__main__':main()
