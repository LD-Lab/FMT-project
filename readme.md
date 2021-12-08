# FMT_function20210415.R
- Run the script to import the function before all the analysis
## Including 
- data preprocessing function
- statistics function
- PCoA analysis
- alpha diversity analysis
- relative abundance analysis
- RNA analysis
- ploting function


# FMT20210415.R 
- import and preprocess refined data after Qiime2 analysis(metadata_all.csv, taxonomy_all.csv, tree_all.nw)
, LEFSE analysis(lefse_data_ctrl_O0_Y0-R_result.csv) or HUMANn2 analysis(meta_genome.csv),and other data
(FMT_DSS_DAI.csv, FMT_DSS_weight.csv, qPCR_data.csv)
- labeled each group with special color, including young and old group, FMT and non FMT group,FMT-DSS and FMT group
analysis including 
- QPCR line plot with significant test,
alpha diversity with Observed Species and Shannon
metagenome pathway analysis with heatmap 
- LEFSE results denoted the significant taxonomy with bar plot.
- weight change with line plot
- DAI change with line plot 


# FMT_relative20210415.R
- select the top richness of taxonomy in Genus based on all the data,sum up the lower level taxonomy to the upper level for performance,
absolute abundance plot based on QPCR.


# FMT_PCoA_20210415.R
- principle coordinate analysis between samples based on bray-curtis or wunifrac distance metric
- plot the sample points and lines based on treatment days  
- make distance analysis based on control day 0 to the treatment days.

# FMT_rna20210415.R
- input refined data from STAR and featureCounts.
(rna-all_counts.csv,rna-metadata.csv,feature_table_all.csv,rna_filter.csv)
obtain DEG_data_all.RData with DEG analysis with all refined data or selected data (DEG_data.RData)
- venn plot based on the DEGs of each group.
- volcano plot based on DEGs of each group.
- heatmap plot based on DEGs of each group.



