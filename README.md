# CS4710

## CS for Bioinformatics
Various programming projects useful for computational biology applications.


### Project 1: String Manipulations
Concatonating DNA sequences, finding percentages of of hydrophobic amino acids in a protein sequence, reverse complement of user-input DNA sequence.

### Project 2: PDB parser – Centroid and distances
Using 3D protein structures to calculate centroid, Euclidean distances, and plotting a historgram of the results.

### Project 3: Protein-Protein Interfaces
Computes protein-protein interfaces and analyzes them using an input of a PDB file of a protein, the identifiers of two chains, and a threshold value. Protein-protein interfaces are determined by computing the Euclidean distance and comparing it to the threshold value.

### Project 4: Computing Properties of Protein-Protein Interaction (PPI) Networks 
A PPI network is a graph where nodes represent proteins and an edge between two nodes represents interacting proteins (either physically or functionally). This program computes the degree distribution and clustering coefficient of a PPI network and compares them to random graphs.

### Project 5: Node Centrality
Computes the closeness centrality property of a Protein-Protein Interaction (PPI) network. The closeness-centrality of node x is defined as the reciprocal of this sum:
![alt text](https://wikimedia.org/api/rest_v1/media/math/render/svg/9acd2f94743dc11086eddc94880b9a6e07774548)
where d(y, x) is the length of the shortest path between x and y.
