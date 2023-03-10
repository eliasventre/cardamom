# CARDAMOM

CARDAMOM is an executable Gene Regulatory Network (GRN) inference method, adapted for timestamped scRNA-seq dataset. The algorithm consists in fitting the parameters of a mechanistic model of gene expression: the simulation of the model, once calibrated, allows to reproduce the dataset used for the inference. The method has been introduced in [1]. It has been benchmarked among other GRN inference tools and applied on a real dataset in [2]. 

# Dependencies

The package CARDAMOM depends on the standard scientific libraries _numpy_ and _scipy_. It also depends on the library _numba_ which accelerates the inference method as well as the simulations. The package _harissa_ is used for the function "simulate_data". The package _alive_progress_ is used to show the progression of the function "simulate_data". Finally, the package _umap-learn_ and _matplotlib_ are used for the function "visualize_data".

They can be installed using pip:

#### pip install harissa
#### pip install umap-learn
#### pip install alive_progress

# Structure of the directories

The user must create a separate directory named "myproject". The example given here corresponds to the dataset of the directory "test/Network4" generated with one of the networks used for the benchmark of [2].

The directory "myproject" must contain 3 directories:

- #### cardamom:
It can be empty before the inference, and will contain the results of the inference, after running the script "infer_network".

- #### Data: 
It must contain 2 files, one named "panel.txt" which contains the dataset for the inference, and one named "panel_genes.txt" which contains the names of the genes.

- #### Rates: 
It must contain a file "degradation_rates.txt" which will be used for simulating the model, when running the script  "simulate_data".

- #### Results: 
It can be empty before the inference, and will contain the results of the simulations, after running the script "visualize_data".

# Before using CARDAMOM

Prior to analyzing an experimental scRNA-seq data with CARDAMOM, We recommend that you perform the following steps before launching CARDAMOM on the count table:

	1.	Select a list of relevant genes,
	2.	Define degradation rates for the corresponding mRNAs and proteins,
	3.	Generate the count table.
 

## 1- Select a list of relevant genes
This is a critical and not so easy task. 
You should aim to obtain a list in the hundred range. CARDAMOM can perform with more but the resulting output GRN will be less and less interpretable as the number of genes grows. For the time being, we recommend starting with say 50 genes, and add more if needed (i.e. if the overall data generation is improved).

The selection step can be performed using any combination of the following:
1.1 Select genes known in the litterature to be important for the process being studied
1.2 Select HVG using custom algorithms like scran or Seurat
1.3  Select ???dynamaically important??? genes using
1.4 Select differentially expressed genes using a basic t-test between two time points.
1.5 Select most delta-entropic genes between two time points.

This selected list of genes will have to be incremented by adding on the first line a ???gene??? called ???Stimulus???. It should take a zero value for the 0h time point and 1 for all the others (see below). 

This gene list will then have to be saved in the myproject/Data folder un the name : panel_genes.txt


## 2- Attribute degradation rates for the corresponding mRNAs and proteins
You can find half-lives for mRNAs and proteins in human in the two papers:
	???	Blumberg et al. 2021): https://doi.org/10.1186/s12915-021-00949-x 
	???	Li et al. 2021).??: https://doi.org/10.1016/j.molcel.2021.09.015 
Only ???short-lived??? proteins half-lives are estimated in the Li paper. In any case, the ???maximum??? half life of proteins will be set by the cell cycle duration. So make sure that no protein has a larger half-life that the cell cycle duration.

A compiled table can be found 
at : https://osf.io/4hqt9/?view_only=23288f5b09274a858cc32009c5a0fe78

You will then have to generate the degradation rates from the half-lives using the following formula: degradation_rate=log(2)/half-life

This degradation rate list will then have to be saved in the myproject/Rates folder un the name : degradation_rates.txt.

## 3- Generate the count table

CARDAMOM will need a count table with cells as columns and genes as rows, which must be in the .txt format. The first row must corresponds to the timepoints at which the cells are sampled, and the first column to the numero of each gene. Then, each line represents the mRNAs counts associated to a gene for each cell at each timepoint. Note that the second line corresponds to the Stimulus, which is set to 0 at t=0h and to 1 at t > 0h (see [1] Section 5.1 for the details). CARDAMOM will be expecting integer values (number of mRNAs molecules).

# Tutorial

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


# References

[1] E. Ventre. ???Reverse engineering of a mechanistic model of gene expression using metastability and temporal dynamics???. In: In Silico Biology 14 (2021), pp. 89???113.

[2] E. Ventre, U. Herbach et al. "One model fits all: combining inference and simulation of gene regulatory networks". In: BioRxiv (2022).

[3] S. Semrau et al. ???Dynamics of lineage commitment revealed by single-cell transcriptomics of differentiating embryonic stem cells???. In: Nat Commun 8 (2017), pp. 1???16.
