PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import os
import unittest
import pandas as pd
import scanpy as sc 
import numpy as np
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


test_network_co_format = pd.DataFrame({"peak": np.nan, "gene": ["AATF", "ALX3", "MYOD1"], "AATF": [1.0,0,1], "ALX3": [1.0,1, 0], "MYOD1": [0.0,0,1]})
class TestNetworkHandlingCO(unittest.TestCase):

  def test_makeNetworkDense(self):
      pd.testing.assert_frame_equal(
          test_network_co_format.copy(), 
          load_networks.makeNetworkDense(
              load_networks.makeNetworkSparse(test_network_co_format, 0.0)
          )
      )
      X = load_networks.makeNetworkDense(
          load_networks.makeRandomNetwork(density=1,  TFs=["AATF", "ALX3", "MYOD1"], targetGenes = ["AATF", "ALX3", "MYOD1"])
      )
  
  def test_pivotNetworkWideToLong(self):
      tallified = load_networks.pivotNetworkWideToLong(
              pd.DataFrame({
                  "gene_short_name": ["a", "b", "c"],
                  "peak": "blerghhh",
                  "A": [1, 0, 0],
                  "B": [1, 1, 0],
                  "C": [1, 0, 1],
              })
          )
      tallified.index = range(5)
      pd.testing.assert_frame_equal(
          tallified,
          pd.DataFrame({
              "regulator": ['A', 'B', 'B', 'C', 'C'],
              "target":    ['a', 'a', 'b', 'a', 'c'],
              "weight": 1,
          })
      )

if __name__=="__main__":
  unittest.main()
