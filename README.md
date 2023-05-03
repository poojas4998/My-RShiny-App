# BF591_Final_Project

#### Data Visualization (GSE64810)

### Data

The data used to build this R Shiny app was extracted from Gene Expression Omnibus (GEO). You can access the dataset here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810

This website utilizes the Huntington's Disease dataset from Labadorf et al., 2016, titled "mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals."

### Tab Information

1. The first tab provides information about the sample metadata. You can find the data in the data/sample_metadata.csv file in this repository.
2. The second tab uses the counts matrix to analyze the count structure and assist in gene filtering for visualizing diagnostic plots, heatmaps, and PCA. The data can be accessed from the data/norm_counts.csv file in this repository.
3. The third tab presents the results from the Differential Expression Analysis, including a summary table and volcano plots. The data can be found in the data/deseq_diff_exp_res.csv file in this repository.
4. The fourth tab enables individual gene expression visualization using the normalized counts and sample metadata. You can access the required data from the data/sample_metadata.csv and data/norm_counts.csv files in this repository.

#### Reference: 

Labadorf, A., Hoss, A., Lagomarsino, V., Latourelle, J., Hadzi, T., Bregu, J., MacDonald, M., Gusella, J., Chen, J., Akbarian, S., Weng, Z., and Myers, R. (2016). Correction: RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. PLOS ONE, 11(7), p.e0160295.
