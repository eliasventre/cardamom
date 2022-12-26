# CARDAMOM

CARDAMOM is an executable Gene Regulatory Network (GRN) inference method, adapted for timestamped scRNA-seq dataset. The algorithm consists in fitting the parameters of a mechanistic model of gene expression: the simulation of the model, once calibrated, allows to reproduce the dataset used for the inference. The method has been introduced in [1]. It has been benchmarked among other GRN inference tools and applied on a real dataset in [2]. 

# Dependencies

CARDAMOM depends on the standard scientific libraries _numpy_ and _scipy_, as well as _matplotlib_. The library _numba_ accelerates the inference method as well as the simulations. The package _harissa_ is used for the function "simulate_data". The package _umap-learn_ is used for the function "visualize_data". Finally the package alive_progress is used to show the progression of the function "simulate_data".

They can be installed using pip:

#### pip install harissa
#### pip install umap-learn
#### pip install alive_progress

# Tutorial

## Structure of the directories

The user must create a separate directory named "myproject". The two examples given here correspond to the dataset of the directory "Network4" generated with one of the networks used for the benchmark of [2], and the dataset of the directory "Semrau" corresponding to the network inferred from an experimental dataset collected in vitro on single mouse embryonic stem cells induced to differentiate by all-trans retinoic acid addition in [3].

The directory "myproject" must contain 3 directories:

- #### cardamom:
It can be empty before the inference, and will contain the results of the inference, after running the script "infer_network".

- #### Data: 
It must contain 2 files, one named "panel.txt" which contains the dataset for the inference, and one named "panel_genes.txt" which contains the names of the genes.

- #### Rates: 
It must contain a file "degradation_rates.txt" which will be used for simulating the model, when running the script  "simulate_data".

- #### Results: 
It can be empty before the inference, and will contain the results of the simulations, after running the script "visualize_data".

## Structure of the data

The dataset must be in the .txt format. The first line must corresponds to the timepoints at which the cells are sampled, and the first column to the numero of each gene. Then, each line represents the mRNAs counts associated to a gene for each cell at each timepoint. Note that the second line corresponds to the Stimulus, which is set to 0 at t=0h and to 1 at t > 0h (see [1] Section 5.1 for the details).

## 1- Calibrating the model from a reference dataset

Run the following script for calibrating the model from the file "myproject/Data/panel.text":

#### python infer_network.py -i [myproject]

The output are the following files:

 "myproject/cardamom/basal_t.npy": matrix of size (nt X ng X 1); containing the basal parameters of the GRN inferred by CARDAMOM at erach timepoint.

 "myproject/cardamom/inter_t.npy": matrix of size (nt X ng X ng); containing the GRN inferred by CARDAMOM at each timepoint.

 "myproject/cardamom/basal.npy": matrix of size (ng X 1); containing the basal parameters of the GRN inferred by CARDAMOM. Note that basal = basal_t[-1].

 "myproject/cardamom/inter.npy": matrix of size (ng X ng); containing the GRN inferred by CARDAMOM. Note that inter = inter_t[-1].

 "myproject/cardamom/kmin.npy": vector of size (ng X 1); containing the minimal burst rates frequency for each gene.

 "myproject/cardamom/kmax.npy": vector of size (ng X 1); containing the maximal burst rates frequency for each gene.

 "myproject/cardamom/bet.npy": vector of size (ng X 1); containing the scaling of the burst sizes for each gene.

"myproject/cardamom/data_bool.npy": matrix of the same size as the data used for the inference, that is used for initializing the simulations.
 
Here, ng denotes the number of genes (including the stimulus), and nt the number of timepoints in the data.
 

## 2- Simulate a dataset from an inferred network

Run the following script for simulating the model from the parameters stored in the directory myproject/cardamom:

#### python simulate_data.py -i [myproject]

The output is the file "myproject/Data/panel_simulated.text".

## 3- Compare the simulations to the reference dataset

Run the following script for comparing the UMAP representations between the dataset "myproject/Data/panel.text" and "myproject/Data/panel_simulated.text":

#### python visualize_data.py -i [myproject]

The outputs are the file "Marginals.pdf" and "UMAP.pdf" that can be found in ./[myproject]/Results.

## References

[1] E. Ventre. “Reverse engineering of a mechanistic model of gene expression using metastability and temporal dynamics”. In: In Silico Biology 14 (2021), pp. 89–113.

[2] E. Ventre, U. Herbach et al. "One model fits all: combining inference and simulation of gene regulatory networks". In: BioRxiv (2022).

[3] S. Semrau et al. “Dynamics of lineage commitment revealed by single-cell transcriptomics of differentiating embryonic stem cells”. In: Nat Commun 8 (2017), pp. 1–16.
