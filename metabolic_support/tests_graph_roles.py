import unittest
import numpy as np
import networkx as nx
from extract_model_graphs import source_seeds,build_graph,route_evidence
from compute_ecological_roles import interaction_matrices,bh,role_codes,pcoa_euclidean,paired_summary,exchange_profiles

class GraphTests(unittest.TestCase):
    def test_source_scc_weights(self):
        g=nx.DiGraph([('a','b'),('b','a'),('b','c')]);w,_=source_seeds(g)
        self.assertEqual(w,{'a':.5,'b':.5})
        g.add_edge('d','a');self.assertEqual(source_seeds(g)[0],{'d':1.})
    def test_singleton_self_loop_matches_stock(self):
        g=nx.DiGraph([('x','x')]);self.assertEqual(source_seeds(g)[0],{})
    def test_bound_direction_and_boundary_exclusion(self):
        r={'id':'R_x','name':'x','substrates':['a'],'products':['b'],'lower':-10,'upper':0,'objective':False,'genes':[]}
        self.assertEqual(set(build_graph({},[r])[0].edges),{('b','a')})
        self.assertEqual(set(build_graph({},[r],stock=True)[0].edges),{('a','b')})
        r['products']=[];self.assertFalse(build_graph({},[r])[0].edges)
    def test_generic_boundary_not_transport(self):
        mets={f'M_a_{c}':{'base_id':'a','compartment':c} for c in ['e','p','c']}
        def rxn(name,sub,prod,genes):return {'id':name,'name':name,'substrates':sub,'products':prod,'lower':0,'upper':1000,'genes':genes,'objective':False}
        ex=rxn('R_EX_a_e',['M_a_e'],[],[]);ex['lower']=-1000
        ep=rxn('R_ep',['M_a_e'],['M_a_p'],[]);pc=rxn('R_pc',['M_a_p'],['M_a_c'],['G_g'])
        routes,*_=route_evidence(mets,[ex,ep,pc],{'G_g'},{'a'})
        self.assertFalse(routes['a']['inward']);self.assertTrue(routes['a']['unsupported_inward'])
        ep['genes']=['G_g'];routes,*_=route_evidence(mets,[ex,ep,pc],{'G_g'},{'a'})
        self.assertTrue(routes['a']['inward']);self.assertFalse(routes['a']['outward'])
        ep['genes']=['G_unknown'];routes,*_=route_evidence(mets,[ex,ep,pc],{'G_g'},{'a'})
        self.assertFalse(routes['a']['inward'])
    def test_matrix_direction(self):
        gs=[{'seed_weights':{'a':1,'b':1},'graph_nodes':['a','b','c']},{'seed_weights':{'a':1},'graph_nodes':['a','b']}]
        c,k=interaction_matrices(gs);self.assertEqual(c[0,1],.5);self.assertEqual(c[1,0],1)
        self.assertEqual(k[0,1],.5);self.assertEqual(k[1,0],0)
    def test_zero_denominators_are_undefined(self):
        gs=[{'seed_weights':{},'graph_nodes':['a']},{'seed_weights':{'b':1},'graph_nodes':['b']}]
        c,k=interaction_matrices(gs)
        self.assertTrue(np.isnan(c[0,1]));self.assertTrue(np.isnan(k[0,1]));self.assertTrue(np.isnan(k[1,0]))
    def test_all_zero_paired_differences(self):
        a=np.ones((4,4));np.fill_diagonal(a,np.nan);mean,p,n,_=paired_summary(a)
        np.testing.assert_array_equal(mean,np.zeros(4));np.testing.assert_array_equal(p,np.ones(4))
        np.testing.assert_array_equal(n,np.full(4,3))
    def test_profiles_preserve_all_zero_whitelist_columns(self):
        gs=[{'seed_weights':{'M_a_e':1},'graph_nodes':['M_a_e'],'metabolites':{'M_a_e':{'base_id':'a'}}},
            {'seed_weights':{},'graph_nodes':['M_a_c'],'metabolites':{'M_a_c':{'base_id':'a'}}}]
        x,columns=exchange_profiles(gs,[{'metabolite_id':'a'},{'metabolite_id':'absent'}])
        self.assertEqual(columns,['uptake__a','uptake__absent','provision__a','provision__absent'])
        np.testing.assert_array_equal(x,[[1,0,0,0],[0,0,1,0]])
    def test_quadrants_and_fdr(self):
        np.testing.assert_allclose(bh([.01,.04,.03]),[.03,.04,.04])
        np.testing.assert_array_equal(role_codes(np.array([1,-1,1,-1,1]),np.array([1,1,-1,-1,1]),np.zeros(5),np.array([0,0,0,0,.2])),[0,1,2,3,4])
    def test_pcoa_distances(self):
        x=np.array([[0,0],[1,0],[1,1]],float);coords,_,_=pcoa_euclidean(x)
        np.testing.assert_allclose(np.linalg.norm(coords[0]-coords[1]),1)

if __name__=='__main__':unittest.main()
