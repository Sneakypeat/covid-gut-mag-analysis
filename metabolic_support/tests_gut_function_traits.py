import unittest
from build_gut_function_traits import trait_values, numeric_count, state, direction, all_evidence, any_evidence


class TraitTests(unittest.TestCase):
    def base(self):
        return {g:0 for g in ["GH13","GH33","GH29","GH95","GH20","GH35","K00634","K00929","K00625","K00925"]}

    def test_pairs_require_both(self):
        d=self.base();d["K00634"]=1
        self.assertFalse(trait_values(d,{"GH13"},[])["butyrate_terminal_pair"])
        d["K00929"]=2
        self.assertTrue(trait_values(d,{"GH13"},[])["butyrate_terminal_pair"])

    def test_acetate_pair(self):
        d=self.base();d["K00625"]=1;d["K00925"]=1
        self.assertTrue(trait_values(d,{"GH13"},[])["acetate_terminal_pair"])

    def test_mucin_requires_three_functions_not_one(self):
        d=self.base();d["GH33"]=1
        self.assertFalse(trait_values(d,{"GH13"},[])["mucin_gh_repertoire"])
        d["GH95"]=1;d["GH35"]=1
        self.assertTrue(trait_values(d,{"GH13"},[])["mucin_gh_repertoire"])

    def test_transport_missing_is_unknown(self):
        self.assertIsNone(trait_values(self.base(),{"GH13"},None)["strict_acid_uptake"])
        self.assertFalse(trait_values(self.base(),{"GH13"},[])["strict_acid_uptake"])
        self.assertTrue(trait_values(self.base(),{"GH13"},["succ"])["strict_acid_uptake"])

    def test_overlap_is_allowed(self):
        d={k:1 for k in self.base()}
        self.assertTrue(all(trait_values(d,{"GH13"},["ac"]).values()))

    def test_missing_markers_not_zero(self):
        self.assertIsNone(trait_values({}, {"GH13"}, [])["plant_backbone_cazyme"])
        self.assertEqual(state(None),"unknown")
        self.assertEqual(state(False),"evidence_not_detected")
        self.assertTrue(any_evidence([True,None]))
        self.assertFalse(all_evidence([False,None]))

    def test_counts_strict(self):
        self.assertIsNone(numeric_count("")); self.assertEqual(numeric_count("2"),2)
        for bad in ["-1","0.1","nan"]:
            with self.assertRaises(ValueError):numeric_count(bad)

    def test_ancom_overlay_only(self):
        self.assertEqual(direction(None),"Not tested")
        self.assertEqual(direction({"diff_GroupCase":"FALSE","lfc_GroupCase":"100"}),"Not significant")
        self.assertEqual(direction({"diff_GroupCase":"TRUE","lfc_GroupCase":"-2"}),"Depleted")


if __name__ == "__main__":
    unittest.main()
