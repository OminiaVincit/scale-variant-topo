# scale-variant-topo
Analysis code for "Scale-variant Topological Information for Characterizing Complex Networks" manuscript.

Preprint at https://arxiv.org/abs/1811.03573

Authors:

* Quoc Hoan Tran: zoro@biom.t.u-tokyo.ac.jp
* Van Tuan Vo: tuan@biom.t.u-tokyo.ac.jp
* Yoshihiko Hasegawa: hasegawa@biom.t.u-tokyo.ac.jp

Hasegawa Lab, Department of Information and Communication Engineering,
Graduate School of Information Science and Technology,
The University of Tokyo

## Generate synthetic network data

We describe source code used in our experiments to generate network data from network models.

### gen-networks/gen_girvan_newman.m (Matlab)
Generate networks from Girvan-Newman network model.
We use the network with 128 nodes partitioned into four communities of size 32 and the average degree of 16.
In this model, we focus on the ratio r between the probability of 
inter-community links (p_out) and intra-community links (p_in).
We vary r=p_out/p_in=0.01, 0.02,...,1.0 and generate 10 random realizations of the network for each value of r.

### gen-networks/gen_WS.m (Matlab)
Generate networks from Watts-Strogatz model.
We use the network with 128 nodes and rewiring probability beta=0.00,0.01,...,1.0, 
with 10 random realizations for each value of beta.

### gen-networks/gen_erdos_renyi.m (Matlab)
Generate networks from Erdos-Renyi model.
We set the number of nodes as 128, vary the pair-link probability p_link=0.020,0.021,...,0.1, 
and generate 10 random realizations of the networks for each value of p_link.

### gen-networks/gen_SP.m (Matlab)
Generate networks from Sales-Pardo model.
We generates hierarchical networks with 640 nodes and three community structures nested in one another.
There are 64 communities of 10 nodes at the small-scale level, 
16 communities of 40 nodes at the medium-scale level 
and 4 communities of 160 nodes at the large-scale level.
There are two parameters; the average degree k_bar, 
which represents the density of the network, 
and rho, which quantifies the separations between the three scales.
We keep k_bar=16, vary the parameter rho=0.05,0.10,...,2.0 
to generate 10 random realizations of the networks at each value of rho.

## Real-world network data

For the classification task, we test on ten real-world network datasets of 
bioinformatics and social networks. 
We save these datasets into two format-types *_adj.mat (adjacency matrix) and *.graphml file
(see folder real-networks/{dataname}).

The first dataset is MUTAG; 
it consists of 188 networks of chemical compounds divided into two classes according to 
their mutagenic effect on a bacterium, 
i.e. mutagenic aromatic and heteroaromatic nitro compounds.

BZR, COX2, DHFR datasets consist of 405 ligands for the benzodiazepine receptor (BZR), 
467 cyclooxygenase-2 inhibitors (COX2) and 756 inhibitors 
of the dihydrofolate reductase (DHFR) come with 3D coordinates.
In each of these datasets, the chemical compounds are divided into 
active and inactive compounds that we need to classify.

FRANKENSTEIN is a modified version of BURSI
dataset made of 4337 molecules with 2401 mutagens and 1936 nonmutagens to classify.
Each molecule is represented as a small network whose vertices are 
labelled by the chemical atom symbol and edges by the bond type.

PROTEINS is a dataset of 1113 proteins represented as networks, 
where the nodes are secondary structure elements and the edges represent 
the neighbourhood within the 3-D structure or along the amino acid chain. 
We classify the proteins into enzymes and non-enzymes.

NCI1 and NCI109 are datasets of chemical compounds with 4110 and 4127 compounds, respectively. 
We classify the compounds with the ability to suppress or inhibit 
the growth of a panel of human tumour cell lines.

IMDB--BINARY is a dataset of 1000 movie collaboration ego-networks for 
each actor (actress) in Action and Romance genres from IMDB.com.
For each network, two nodes representing actors or actresses 
are connected if they appear in the same movie.
The task is to identify whether a given ego-network of 
an actor (actress) belongs to the Action or Comedy genre.

IMDB--MULTI is the multi-class version of IMDB--BINARY 
and contains 1500 ego-networks belonging to Comedy, Romance and Sci-Fi genres.

## Compute diffusion distance matrix from the network

###  generate_diffusion_mat.py (Python)

Generate diffusion distance matrix from network adjency matrix.

    Options:
            --matpath : input directory consisting of the collection of adjency matrices 
                (nodes_*_index_{}_adj.mat --> for ER, GN, SP, WS, Drosophila-Melanogaster networks 
                nodes_*_index_{}_adj.nse  --> for LFR, LFR-H networks )
            --outpath : output directory to save diffusion matrices
                (*.npy file with matrix size of K x N x N, K: number of tau, N: number of nodes)
            --normflag : normalize (=1) the diffusion matrix by its maximum element or not
            --avgflag  : average (=1) the diffusion matrix for K values of tau or not
            --idxbg    : start index of adjency matrices file to divide in multi processing
            --idxed    : end index of adjency matrices file to divide in multi processing
            --taubg    : start value of tau
            --taued    : end value of tau
            --interval : sampling interval for tau from taubg to taued
            --nproc    : number of processes using in this program 

    Example: python generate_diffusion_mat.py 
                        --matpath get-networks\mat\GN-net --outpath distmat\GN-net 
                        --normflag 0  --avgflag 0  --idxbg 1      --idxed 1
                        --taubg 1.0   --taued 10.0 --interval 1.0 --nproc 8

## Compute three-dimensional persistence diagrams and persistence kernels for networks

### Compile by Visual Studio 2015

### Release/PersistenceRunner_D64.exe: calculate persistence diagrams from point clouds

### Release/DiffusionRunner_D64.exe: calculate persistence diagrams from diffusion distance matrix

### Release/DiagramDistanceRunner_D64.exe: calculate kernel for persistence diagrams



## Compute graph kernel for networks
### graph-kernels/gen_igraph.R  (R)
Generate igraph database files from .graphml files.
Example: 

    genigraph("MUTAG")


to convert all *.graphml files in data/MUTAG/graphml to list of graphs: jointddf_MUTAG.RData.

### graph-kernels/calculate_kernel_igraph.R
Calculate graph-kernels from igraph database files.
Example: 

    datname <- "MUTAG"
    calckernel()

## Compute common features for networks
### make_common_features.py (Python)

Calculate common features from network adjency matrix


    16 common features:
        avg_deg_centrality, avg_clustering_coeff, avg_shortest_paths, 
        avg_eccentricity, avg_node_closeness_centrality, 
        avg_node_betweeness_centrality, avg_edge_betweeness_centrality, 
        avg_local_efficiency, avg_triangles,
        density, transitivity, diameter, radius, 
        assortiative_coef, global_efficiency, nums_connected_parts
    
    Usage:   python make_common_features.py -d dataname -i infolder -o outfolder
    Example: python make_common_features.py -d MUTAG -i real-networks -o common-features

