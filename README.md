# config_resistance
Will be public repository that computes effective resistance from the weighted graph Laplacian

This code is written to run on the data files generated from the code in the config_generate repository

# Lloyd_and_resistance_parfor_v2.m

Main script to use the generated point clouds and compute the effective resistance (and voltage at each node, current at each edge).

  **find_corners_adjacency_A** - finds the NE node and makes this node 1, the SW node and makes this node N, then connects nodes based on the Delaunay triangulation.  

  **find_corners_adjacency_B** - finds the NW node and makes this node 1, the SE node and makes this node N, then connects nodes based on the Delaunay triangulation.  

  **compute_voltage_Adj** - defines resistance as the Euclidean distance between connected nodes and applies a current across the network from node 1 to node N

Within, can define parameters for comparison to experimental results:

  xarea: the cross sectional area of one beam of the network in cm^2

  tmp: the resistance of each beam of the network in mOhm

  Returns the effective resistance across the entire network, as well as the voltage at each node and the current across each edge
