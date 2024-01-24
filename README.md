Genetic-SCP-Diplodus-Mullus

# Evaluating different approaches to integrate genome-wide genetic diversity in spatial conservation prioritization

This repository contains code to perform all the analyses presented in the manuscript. The analysis is subdivided in the following steps, each corresponding to a code file and detailed below, to be executed in this order:

### 01 - Create_pus.R

This code creates the planning units (PUs): 5203 square PUs (10 km side) spanning the Mediterranean Sea over coastal areas (<200 m depth).
It extracts the presence/absence of the two species from the FishMed database.
It identifies the PUs already protected according to the marine reserves identified by Abecasis et al (2023) and calculates the surface area protected (Table S1).
It imports the conservation costs of Mazor et al (2014).

### 02 - Run_PCA.R

This script performs the principal component analysis (PCA) on the genomic datasets of Boulanger et al (2022): 8068 SNPs for D. sargus and 2753 SNPs for M. surmuletus; and plots the first two axes of the PCA (Figure S5).
For each axis of the PCA, it interpolates the scores of the sampled sites on a raster spanning the study region using inverse distance weighting interpolation.

### 03 - Compare_spatial_interpolation_PCA.R

This script compares two spatial interpolation models for PCA scores: (i) nearest-neighbor interpolation and (ii) inverse distance weighting interpolation using k-fold cross validation (Figure S9).

### 04 - Find_number_clusters.R

This scripts finds the optimal number of clusters for the multidimensional discrete genetic cluster approach, by applying several clustering indices to the data (Table S2).

### 05 - Evaluation_extension_current_Mr.R

This script defines a raptr problem with one multidimensional genetic space for each of the two species using one demand point for each planning unit where the species is present.
It evaluates the amount and space held by the current set of marine reserves.
It formulates and solves (using prioritizr) a prioritization problem to extend the current set of marine reserves to cover 15% of the species' ranges.

### 06 - Run_scenarios_raptr.R

This script defines the planning problems for raptr and solves them.

A first problem using 100% demand points, which is defined, but not solved, due to memory errors. This problem is the "gold standard" formulation and its genetic spaces and demand points are subsequently used to evaluate the solutions found in the other prioritization problems (see 08 - Analyze_scenarios.R")

Two approximate problems using 50% and 20% of the original sets of demand points, replicated 20 times.

### 07 - Run_scenarios_prioritizr.R

This script defines the planning problems for prioritizr and solves them. It consists of three parts:

1) Definition and solution of problems using unidimensional discrete conservation features: these problems consider 4 classification methods (quantile interval, equal interval, natural breaks and standard deviation) and 3 numbers of classes (only 1 for the standard deviation).

2) Partitioning around medoids: applies this multivariate clustering method to find clusters in the Multidimensional discrete genetic cluster approach 

3) Definition and solution of problems using multidimensional discrete genetic clusters: these problems consider 3 numbers of clusters (low, medium and high). 

### 08 - Analyze_scenarios.R

This script analyses the results of the prioritization problems. It calculates:

 - Jaccard distances between solutions and their statistical significance through a permutation approach;
 
 - genetic space held (relative to the genetic spaces and demand points of the "gold standard" approach) for across all problems and solutions;
 
 - total conservation cost.
 
 - Maps of selection frequency for each problem (Figure S10) and for the approximate problem with the best performances, i.e. raptr_50gs (Figure 5)


### 09 - Mapping

This script contains assorted codes to plot some of the maps shown in the paper, i.e.:

 - Figure S1: map of conservation cost
 
 - Figure S2: map of sampling points
 
 - Figure S6 and S7: spatially interpolated PCA scores for D. sargus and M. surmuletus
 
 - Figure S8: map of priority sites for protecting 15% of species ranges without explicit genetic objectives. 

<br/><br/>

## /data
In the /data folder, the following files are needed to prepare the data:

### 01 - Build_sampling_dataframe.R

This script produces the “cell_sampling.RData” file containing sampling coordinates (and other data) for each sampling cell

### 02 - Merge_neutral_adaptive.R

This script merges the original datasets of neutral and adaptive loci (from Boulanger et al 2022) into a single dataset, for each of the two species

<br/><br/>

## /function
In the /function folder, there is a single file containing one helper function:

### split.taxon.R
It contains the function split.taxon, used to split a layer into distinct conservation features using single PCA axes using different methods and number of features.
