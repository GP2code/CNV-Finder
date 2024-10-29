
# CNV-Finder Available Parameters

This guide walks you through the parameters used in the 3 major proceses carried out by CNV-Finder: ML data preprocessing, application of pre-trained models/training new models on the prepared data, and the creation of app-ready files for visualizations in the CNV-Finder app. Each parameter is configurable and includes default values.

## Table of Contents
- [ML Data Prep Parameters](#ml-data-prep-parameters)
  - [Interval and Region of Interest](#interval-and-region-of-interest)
  - [Analysis Configuration](#analysis-configuration)
  - [Input Files](#input-files)
  - [General Configurations](#general-configurations)
  - [Testing Configuration](#testing-configuration)
- [LSTM Model Application Parameters](#lstm-model-application-parameters)
  - [Data Files](#data-files)
  - [Model Training and Prediction](#model-training-and-prediction)
  - [Output](#output)
- [App Prep Parameters](#app-prep-parameters)
  - [Interval and Region of Interest](#interval-and-region-of-interest-1)
  - [Additional Testing and Metrics Files](#additional-testing-and-metrics-files)
  - [Model Probability Threshold](#model-probability-threshold)
  - [Output and CPUs](#output-and-cpus)
  - [App Preparation](#app-preparation)
- [Example Usage](#example-usage)

## ML Data Prep Parameters

### Interval and Region of Interest
- `--interval_name`: *(str)* Name for NDD-related gene region (default: `None`). If not found in reference file, it will be added to `CNV_modules/ref_files/custom_intervals.csv` for future use (e.g., 22q_small custom interval).
- `--interval_file`: *(str)* Path to gene or feature intervals file (default: `CNV_modules/ref_files/glist_hg38_intervals.csv`). File format: `[NAME, CHR, START, STOP]`, one line per interval.
- `--chrom`: *(str)* Chromosome for the region of interest (default: `None`).
- `--start`: *(int)* Starting position in base pairs (hg38) (default: `None`).
- `--stop`: *(int)* Stopping position in base pairs (hg38) (default: `None`).
- `--buffer`: *(int)* Size of kilobase window around each interval, in base pairs (default: `250000`).

### Analysis Configuration
- `--min_gentrain`: *(float)* Minimum GenTrain Score threshold (default: `0.2`).
- `--split_interval`: *(int)* Number of intervals to split each region (before overlap) (default: `5`).
- `--total_windows`: *(int)* Count of windows in each region per sample (with overlap) (default: `50`).

### Input Files
- `--bim_file`: *(str)* Path to PLINK `.bim` file after sample QC (default: `None`).
- `--pvar_file`: *(str)* Path to PLINK2 `.pvar` file after sample QC (default: `None`).
- `--master_file`: *(str)* Master key file for available samples in data release (default: `ref_files/master_key.txt`).
- `--training_ids`: *(str)* List of training set IDs with `IID` headers (default: `ref_files/training_set_IDs.csv`).
- `--testing_ids`: *(str)* List of testing set IDs with `IID` headers (default: `None`).

### General Configurations
- `--study_name`: *(str)* Cohort or subset from larger data release (default: `all`).
- `--metrics_path`: *(str)* Path to SNP metrics files in parquet format (default: `ref_files/snp_metrics`).
- `--out_path`: *(str)* Path for output reports (default: `None`). Suggested format: `{Cohort}_{Gene or Interval Name}`.
- `--cpus`: *(int)* Number of CPUs available for job (default: `8`).

### Testing Configuration
- `--create_testing`: *(boolean)* Create a new testing set with no overlaps in training set IDs (default: `True`).
- `--test_size`: *(int)* Number of samples in testing set (default: `500`).

## LSTM Model Application Parameters

### Data Files
- `--train_file`: *(str)* Path to training file with windowed samples (default: `None`).
- `--test_file`: *(str)* Path to testing file with windowed samples (default: `None`).

### Model Training and Prediction
- `--train`: *(boolean)* Train a new model when set to `True` (default: `True`).
- `--feature_names`: *(str)* Names of features for model training, separated by spaces (e.g., `avg_baf avg_mid_baf avg_lrr`).
- `--model_file`: *(str)* Path to saved model file (.keras) (default: `None`). Only .keras format is currently supported due to HPC limitations.
- `--predict`: *(boolean)* Generate predictions with the model, either newly trained or pre-saved (default: `True`).
- `--print_summary`: *(boolean)* Display the architecture of the ML model (default: `True`).

### Output
- `--out_path`: *(str)* Path to output reports (default: `None`). Suggested format: `Cohort_[Gene/Interval Name]`.

## App Prep Parameters

### Interval and Region of Interest
Parameters are the same as those under the **ML Data Prep** section.

### Additional Testing and Metrics Files
- `--test_set_ids`: *(str)* Path to file with testing set IDs and SNP metrics path generated in ML Data Prep section (default: `None`).
- `--test_set_windows`: *(str)* Path to file with windowed test set generated in ML Data Prep section (default: `None`).
- `--test_set_results`: *(str)* Path to file with model results on test set generated in LSTM Model Application section (default: `None`).

### Model Probability Threshold
- `--probability`: *(float)* Probability threshold for model predictions to include in app (default: `0.8`).

### Output and CPUs
- `--out_path`: *(str)* Path to app-ready report output (default: `None`). Suggested format: `Cohort_[Gene/Interval Name]`.
- `--cpus`: *(int)* Number of CPUs available for job (default: `8`).

### App Preparation
- `--make_app_ready`: *(boolean)* Create an app-ready file with all testing samples, necessary for app creation (default: `True`). 

--- 

For detailed usage of above parameters, consult the `run_pipeline.ipynb` notebook at https://github.com/nvk23/CNV-Finder.
