# CNV-Finder
 
## Overview
CNV-Finder is a novel pipeline integrating a Long Short-Term Memory network don SNP array data to expedite large-scale identification of CNVs within predefined genomic regions. Manuscript coming soon!

## Installation
### Clone the repository:

````
git clone https://github.com/nvk23/CNV_finder.git

cd CNV_finder
````

### Install the required packages:

````
pip install -r requirements.txt
````

## Running the Pipeline
Open run_pipeline.ipynb and run through each cell for the 3 major processes in the pipeline: data preparation for machine learning, application of pre-trained models onto ML-ready data, and creation of app-ready files for Streamlit visualizations. If HPC is unavailable, run the commands defined as "cmd" in terminal. 

### *Available Parameters*
#### **ML Data Prep**
--interval_name, type=str, default=None, Name for NDD-related Gene region. If not found in reference file, will be added to ref_files/custom_intervals.csv for future reference/file names (e.g. 22q_small custom interval).

--interval_file, type=str, default='ref_files/glist_hg38_intervals.csv', Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per interval. Autosomes only.

--chrom, type=str, default=None, Chromsome for region of interest.

--start, type=int, default=None, Starting position in base pairs (hg38).

--stop, type=int, default=None, Stopping position in base pairs (hg38).

--buffer, type=int, default=250000, Kilobase window around each interval, in bases.

--min_gentrain, type=float, default=0.2, Minimum GenTrain Score threshold.

--split_interval, type=int, default=5, Number of intervals to split region into (before overlap).

--total_windows, type=int, default=50, Count of windows in region per sample (with overlap).

--bim_file, type=str, default=None, PLINK .bim file following sample QC.

--pvar_file, type=str, default=None, PLINK2 .pvar file following sample QC.

--master_file, type=str, default='ref_files/master_key.txt', help='Master key for all available samples in data release (need IID and binary GDPR (0-no, 1-yes) column).

--study_name, type=str, default='all', Subset or cohort of interest from larger data release.

--metrics_path, type=str, default='ref_files/snp_metrics', help='Path to SNP metrics files/parquets with the format {metrics_path}/{sentrixbarcode}/snp_metrics_{sentrixbarcode}/Sample_ID={sample}.

--out_path, type=str, default=None, Path to output reports with suggested format {Cohort}_{Gene or Interval Name}.

--cpus, type=int, default=8, Number of CPUs available for the job.

--training_ids, type=str, default='ref_files/training_set_IDs.csv', List of IDs used in training set with headers IID. May also include the interval where CNV was found in each sample. Otherwise interval/chromosome with position range must be included in proper labels.

--testing_ids, type=str, default=None, List of IDs used in testing set with header IID.

--create_testing, type=boolean, default=True, Create new testing set with no overlaps in training set IIDs.

--test_size, type=int, default=500, Number of samples to include in testing set.

#### **LSTM Model Application**
--train_file, type=str, default=None, Path to training file with windowed samples.

--test_file, type=str, default=None, Path to testing file with windowed samples.

--train, type=boolean, default=True, Train new model when True.

--feature_names, type=str, Add feature names for model training separated by spaces (eg.avg_baf avg_mid_baf avg_lrr).

--model_file, type=str, default=None, Path to saved model file (.keras). HPC currently has issues handling other model save formats.

--predict, type=boolean, default=True, Generate prediction results with a newly trained or pre-saved model.

--print_summary, type=boolean, default=True, Display architecture of ML model.

--out_path, type=str, default=None, Path to output reports with suggested format Cohort_Gene or Interval Name.


#### **App Prep**
--interval_name, type=str, default=None, Name for NDD-related Gene region.

--interval_file, type=str, default='ref_files/glist_hg38_intervals.csv', Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per interval. Autosomes only.

--chrom, type=str, default=None, Chromsome for region of interest.

--start, type=int, default=None, Starting position in base pairs (hg38).

--stop, type=int, default=None, Stopping position in base pairs (hg38).

--buffer, type=int, default=250000, Kilobase window around each interval, in bases.

--min_gentrain, type=float, default=0.2, Minimum GenTrain Score threshold.

--bim_file, type=str, default=None, PLINK .bim file following sample QC.

--pvar_file, type=str, default=None, PLINK2 .pvar file following sample QC.

--test_set_ids, type=str, default=None, Path to file with testing set IIDs and SNP metrics path (headers: IID, snp_metrics_path).

--test_set_windows, type=str, default=None, Path to file with windowed test set used in training.

--test_set_results, type=str, default=None, Path to file with model results on test set.

--probability, type=float, default=0.8, Probability threshold of model predictions that must be met to be included in app.

--out_path, type=str, default=None, Path to output app-ready report with suggested format Cohort_Gene or Interval Name.

--cpus, type=int, default=8, Number of CPUs available for the job.

--make_app_ready, type=boolean, default=True, Create 1 app ready file including all testing samples with necessary info for app creation.
