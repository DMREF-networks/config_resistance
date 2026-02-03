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

  **find_corners_adjacency_voronoi** - can take the place of `find_corners_adjacency_A` to construct a Voronoi tessellation of the point cloud instead of a Delaunay triangulation.

# Edge flip test

For a given network, identify two neighboring triangles that form a convex quadrilateral structure with their common edge being one of the two quadrilateral diagonals.  Switch the diagonal for the other diagonal and record the resulting change in effective resistance measured from the NE to SW corners and the change in total effective resistance.  

  **edge_flip_test.m** - takes in the adjacency matrix and performs the edge flip test.

  **isConvex_from_coords.m** - uses the (x,y) coordinates of the nodes in a quadrilateral and determines if it is convex.  

  **compute_Rs.m** - an updated version of `compute_voltage_Adj` that returns the effective resistance between any two nodes, as well as the total effective resistance.

  **draw_flip_heatmap_edges_only.m** - plotting routine to plot the result of the edge_flip_test on the network.  Calls the colormap defined in `purple_cmap.m`.

# License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

If you use this code, please also cite Obrero et al., (2025) Electrical transport in tunably disordered metamaterials, *Phys. Rev. E* 112, 035505
https://doi.org/10.1103/6bph-n6zj

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
