# config_resistance
Given a set of (x,y) coordinates of a point cloud in 2D, this code creates a network metamaterial by connection points via the Delaunay triangulation in a bounding box.  It then computes the effective resistance across either diagonal of the material, as well as the option to return the voltage at each node and the current flow along each edge.

This code is written to run on the data files generated from the code in the config_generate repository.

# Lloyd_and_resistance_parfor_v2.m

Main script to use the generated point clouds and compute the effective resistance (and voltage at each node, current at each edge).

  **find_corners_adjacency_A** - finds the NE node and makes this node 1, the SW node and makes this node N, then connects nodes based on the Delaunay triangulation.  

  **find_corners_adjacency_B** - finds the NW node and makes this node 1, the SE node and makes this node N, then connects nodes based on the Delaunay triangulation.  

Inspiration for the algorithm to create the network Adjacency matrix was taken from 
from https://people.sc.fsu.edu/~jburkardt/presentations/voronoi_neighbors.pdf

  **compute_voltage_Adj** - computes the weighted adjacency matrix and computes effective resistance, voltages, and currents.  Defines resistance as the Euclidean distance between connected nodes and applies a current across the network from node 1 to node N.  Within, can define parameters for comparison to experimental results:

  xarea: the cross-sectional area of one beam of the network in cm^2

  tmp: the resistance of each beam of the network in mOhm

# License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
