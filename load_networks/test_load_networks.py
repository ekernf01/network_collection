PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import os
import unittest
import pandas as pd
import scanpy as sc 
os.chdir(PROJECT_PATH + "network_collection")

import sys
import importlib
sys.path.append("load_networks")
import load_networks
importlib.reload(load_networks)
os.environ["GRN_PATH"] = "networks"

class TestSimpleFunctions(unittest.TestCase):
    def test_simple_loaders(self):
      self.assertIsInstance(load_networks.load_grn_metadata(), pd.DataFrame)
      self.assertIsInstance(load_networks.list_subnetworks("gtex_rna"), list)
      self.assertIsInstance(load_networks.load_grn_by_subnetwork("gtex_rna", "Adipose_Subcutaneous.parquet"), pd.DataFrame)

class TestLightNetwork(unittest.TestCase):
    def test_LightNetwork(self):
      self.assertIsInstance(
        load_networks.LightNetwork("gtex_rna"), 
        load_networks.LightNetwork
        )
      self.assertIsInstance(
        load_networks.LightNetwork("gtex_rna", ["Adipose_Subcutaneous.parquet"]), 
        load_networks.LightNetwork
        )
      demo_df = pd.DataFrame({'regulator':'not_a_regulator', 'target':'not_a_target', 'weight':-1, "cell_type":"psc"}, index = [0])
      self.assertIsInstance(
        load_networks.LightNetwork(df = demo_df), 
        load_networks.LightNetwork
        )
      self.assertIsInstance(
        load_networks.LightNetwork("gtex_rna", df = demo_df), 
        load_networks.LightNetwork
        )
      self.assertIsInstance(
        load_networks.LightNetwork("gtex_rna"), 
        load_networks.LightNetwork
        )
      self.assertRaises(ValueError, load_networks.LightNetwork("gtex_rna").save, "temp_file")
      self.assertTrue(load_networks.LightNetwork("gtex_rna").save("temp_file178023571364.parquet") is None)
      os.unlink("temp_file178023571364.parquet")
      self.assertIsInstance(load_networks.LightNetwork("gtex_rna").get_all(), pd.DataFrame)
      self.assertIsInstance(load_networks.LightNetwork("gtex_rna").get_regulators("GAPDH"), pd.DataFrame)
      self.assertIsInstance(load_networks.LightNetwork("gtex_rna").get_all_regulators(), set)
      self.assertIsInstance(load_networks.LightNetwork("gtex_rna").get_all_one_field("cell_type"), set)
      self.assertTrue("psc" in load_networks.LightNetwork(df=demo_df).get_all_one_field("cell_type"))
      self.assertIsInstance(load_networks.LightNetwork("gtex_rna").get_num_edges(), int)
      self.assertTrue(all("GAPDH"        == load_networks.LightNetwork("gtex_rna").get_regulators("GAPDH")["target"]))
      self.assertTrue(all("not_a_target" == load_networks.LightNetwork(df=demo_df).get_regulators("not_a_target")["target"]))

if __name__=="__main__":
  unittest.main()