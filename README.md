# CNV-Finder
 
## Overview
CNV-Finder is a novel pipeline integrating a Long Short-Term Memory network on SNP array data to expedite large-scale identification of CNVs within predefined genomic regions. Manuscript coming soon!

## Project Structure
```
CNV_finder/
│
├── app/
│   ├── data/
│   │   └── gp2_cohort
│   │   │   └── app
│   ├── Home.py
│   ├── variant_plots.py
├── modules/
│   ├── ref_files/
│   │   └── models/
│   │   │   └── final_del_5_50_combo4_lstm.keras
│   │   │   └── final_dup_10_70_combo6_lstm.keras
│   │   │   └── updated_del_5_50_combo4_lstm.keras
│   │   │   └── updated_dup_10_70_combo6_lstm.keras
│   │   │   └── prelim_del_5_50_combo4_lstm.keras
│   │   │   └── prelim_dup_10_70_combo6_lstm.keras
│   │   └── custom_intervals.csv
│   │   └── glist_hg38_interval.csv
│   │   └── training_set_IDs.csv
│   ├── cnv_finder/
│   │   └── data_methods.py
│   │   └── model_methods.py
│   ├── run_data_prep.py
│   ├── run_lstm_model.py
│   ├── run_app_prep.py
├── docs/
│   ├── parameter_guide.md
├── run_pipeline.ipynb
├── requirements.txt
├── README.md
```

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
Open the `run_pipeline.ipynb` notebook and sequentially run through each cell to perform the 3 major processes in the pipeline: ML data preprocessing, application of pre-trained models/training new models on the prepared data, and the creation of app-ready files for visualizations in the CNV-Finder app. If HPC is unavailable, run the commands defined as "cmd" in the terminal. 

### Available Parameters
For a more in-depth guide to the parameters available for each process, please read through the following documentation: `docs/parameter_guide.md`. 