# Tumor detection by methylation and hemi-methylation of cfDNA



## Requirements
R version 4.3.0 (2023-04-21) and packages:
* keras_2.11.1
* GenomicRanges_1.52.0
* rtracklayer_1.60.0
* data.table_1.14.8
* openxlsx_4.2.5.2
* pheatmap_1.0.12
* RColorBrewer_1.1-3
* ggplotify_0.1.0
* cowplot_1.1.1
* pROC_1.18.2
* ggplot2_3.4.2
* ggrepel_0.9.3
* plyr_1.8.8
* dplyr_1.1.2
* BSgenome.Hsapiens.UCSC.hg19_1.4.3
* reshape2_1.4.4

## Installation
Copy this folder to your computer. It takes ~5 min.


## Patient and Sample information
* "56_validation_set.xlsx": the validation sample list.
* "215_training_set.xlsx": the training sample list.


## Models
* "QSEA_diff/215_training_set_other_nature_block_QSEA": The folder of identified DMR from training set.
* "DMBR_diff/215_training_set_other_nature_block_RPM1_each0.3_miss1": The folder of identified DHMR from training set.
* "DMR_model": Include cancer specific DMR models from 10 subsets of training samples.
* "DHMR_model": Include cancer specific DHMR models from 10 subsets of training samples.
* "Calibration_model": Include calibration model combined by DMR and DHMR models.

## Output
* "Calibration_prediction": The folder include prediction results of DMR, DHMR and calibration model.
* "./DMR_model_pred56_raw/GLMNET_model_test": Include prediction results summary by DMR model using GLMNET model.
* "./DHMR_model_pred56_raw/GLMNET_model_test": Include prediction results summary by DHMR model using GLMNET model.
* "./Calibrate_model_vali56_GLMNET_raw/calibration_test": Include prediction results summary by DMR+DHMR.


## Script
* "plasma_MCD.R": Run "Rscript plasma_MCD.R" can directly run the program.
The running may take 30 minutes depending on the computer configuration. The expceted result is the predicted probabilty for each sample and ROC curves shown as in the manuscript.


