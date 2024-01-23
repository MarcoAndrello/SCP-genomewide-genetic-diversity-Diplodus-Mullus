# Genetic-SCP-Diplodus-Mullus
Evaluating different approaches to integrate genome-wide genetic diversity in spatial conservation prioritization

This repository contains code to perform all the analyses presented in the manuscript. The analysis is subdivided in the following steps, each corresponding to a code file and detailed below, to be executed in this order:

01 - Create_pus.R

This code creates the planning units (PUs): 5203 square PUs (10 km side) spanning the Mediterranean Sea over coastal areas (<200 m depth).
It extracts the presence/absence of the two species from the FishMed database.
It identifies the PUs already protected according to the marine reserves identified by Abecasis et al (2023) and calculates the surface area protected (Table S1).
It imports the conservation costs of Mazor et al (2014).


02 - Run_PCA.R

This script performs the principal component analysis (PCA) on the genomic datasets of Boulanger et al (2022): 8068 SNPs for D. sargus and 2753 SNPs for M. surmuletus; and plots the first two axes of the PCA (Figure S5).
For each axis of the PCA, it interpolates the scores of the sampled sites on a raster spanning the study region using inverse distance weighting interpolation.


03 - Compare_spatial_interpolation_PCA.R

This script compares two spatial interpolation models for PCA scores: (i) nearest-neighbor interpolation and (ii) inverse distance weighting interpolation using k-fold cross validation (Figure S9).


04 - Find_number_clusters.R

This scripts finds the optimal number of clusters for the multidimensional discrete genetic cluster approach, by applying several clustering indices to the data (Table S2).


05 - Evaluation_extension_current_Mr.R

This script defines a raptr problem with one multidimensional genetic space for each of the two species using one demand point for each planning unit where the species is present.
It evaluates the amount and space held by the current set of marine reserves.
It formulates and solveS (using prioritizr) a prioritization problem to extend the current set of marine reserves to cover 15% of the species' ranges.


MI FERMO QUI:
MODIFICARE RUN SCENARIOS RAPTR PER FAR CORRISPONDERE I NOMI DEI FILES DEI RISULTATI



06 - Run_scenarios_raptr.R

This script defines the planning problems for raptr and solves them.
It defines the 13 problems (10 single-PCA and 3 multi-PCA problems) by defining the attribute spaces, placing the demand points and setting the targets, then solves the problems with raptr using Gurobi.
The 10 single-PCA problems consider 4 methods to place demand points (quantile interval, equal interval, natural breaks and standard deviation) and 3 numbers of demand points (only 1 for the standard deviation).
The 3 multi-PCA problems consider 3 different numbers of demand points, placed using a hypervolume approach.


03a - Scenarios_prioritizr.R

This script defines the planning problems for prioritizr and solves them. It first runs the k-means analysis to find the optimal number of clusters for the multi-PCA approaches (Figure S1).
It defines the 12 problems (10 single-PCA and 2 multi-PCA problems) by defining the conservation features and setting the targets, then solves the problems with prioritizr using Gurobi.
The 10 single-PCA problems consider 4 classification methods (quantile interval, equal interval, natural breaks and standard deviation) and 3 numbers of classes (only 1 for the standard deviation).
The 2 multi-PCA problems consider 2 different numbers of clusters.





04a - Analyze_scenarios_prioritizr.R

This script analyses the results of the planning problems solved with prioritizr. It calculates:
 - Jaccard distances between solutions and their statistical significance through a permutation approach
 - "amount held" for each conservation feature across all problems and solutions
 - total conservation cost (number of selected PUs).
It plots Figure 1a and 1b, and Figure 2a and 2b


04b - Analyze_scenarios_raptr.R

This script analyses the results of the planning problems solved with raptr. It calculates
 - Jaccard distances between solutions and their statistical significance through a permutation approach
 - "space held" for each genetic space across all problems and solutions
 - total conservation cost (number of selected PUs).
 - maximum target for each genetic space
It plots Figure 1c and 1d, and Figure 2c and 2d


05 - Evaluate_current_MR_network.R

This script calculates
 - the "amount held" for the conservation features of each planning problem of prioritizr
 - the "space held" for the genetic spaces of each planning problem of raptr
by the current system of marine reserves of the Mediterranean Sea.
It plots Figure 3.


06 - Extend_current_MR_network.R

This script defines a planning problem with prioritizr to extend the current system of marine reserves to reach the 15% targets of representation for genome-wide genetic diversity, using the 12 planning problems of prioritizr.
It calculates total conservation costs of the 12 solutions and plots the consensus map of priority sites (Figure 4)
It compares the solutions with that of a species-only SCP.



In the /data folder, the following files are needed to prepare the data:


01 - Build_sampling_dataframe.R

This script produces the “cell_sampling.RData” file containing sampling coordinates (and other data) for each sampling cell


02 - Merge_neutral_adaptive.R

This script merges the original datasets of neutral and adaptive loci (from Boulanger et al 2022) into a single dataset, for each of the two species


In the /function folder, the split.taxon.R file contains two functions

split.taxon, used to split a layer into distinct conservation features using single PCA axes (in prioritizr) or to find the coordinates of demand points on single PCA axes (In raptr) using different methods and number of features/demand points.

Split.taxon.multi to split a layer into distinct conservation features using multiple PCA axes (in prioritizr) using a k-means algorithm

