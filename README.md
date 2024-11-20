# CNV-Finder

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/890546338.svg)](https://doi.org/10.5281/zenodo.14182563)

## Overview
CNV-Finder is a novel pipeline integrating a Long Short-Term Memory network on SNP array data to expedite large-scale identification of CNVs within predefined genomic regions. Manuscript coming soon!

## Installation
### Clone the repository:

````
git clone https://github.com/nvk23/CNV-Finder.git

cd CNV-Finder
````

### [Optional] Create a Conda environment:

````
conda create -n "cnv_finder" python=3.9.16 ipython

conda activate cnv_finder
````

### Install the required packages:

````
pip install -r requirements.txt
````

## Running the Pipeline
Open the `run_pipeline.ipynb` notebook and sequentially run through each cell to perform the 3 major processes in the pipeline: ML data preprocessing, application of pre-trained models/training new models on the prepared data, and the creation of app-ready files for visualizations in the CNV-Finder app.

### Available Parameters
For a more in-depth guide to the parameters available for each process, please read through the following documentation: `docs/parameter_guide.md`. 

## Project Structure
```
CNV_finder/
├── app/
│   ├── selections/
│   ├── Home.py
│   └── variant_plots.py
│
├── modules/
│   ├── cnv_finder/
│   │   ├── data_methods.py
│   │   └── model_methods.py
│   ├── run_data_prep.py
│   ├── run_lstm_model.py
│   └── run_app_prep.py
│
├── ref_files/
│   ├── models/
│   │   ├── final_del_5_50_combo4_lstm.keras
│   │   ├── final_dup_10_70_combo6_lstm.keras
│   │   ├── updated_del_5_50_combo4_lstm.keras
│   │   ├── updated_dup_10_70_combo6_lstm.keras
│   │   ├── prelim_del_5_50_combo4_lstm.keras
│   │   └── prelim_dup_10_70_combo6_lstm.keras
│   ├── custom_intervals.csv
│   ├── glist_hg38_interval.csv
│   └── training_set_IDs.csv
│
├── example_data/
│   ├── snp_metrics/
│   └── test_master_key.csv
│
├── testing/
│   ├── app_ready/
│   │   ├── cohort/
│   │   │   ├── final_del_model/
│   │   │   └── final_dup_model/
│   ├── del/
│   │   ├── cohort/
│   │   │   └── gene/
│   └── dup/
│       └── cohort/
│           └── gene/
│
├── docs/
│   └── parameter_guide.md
│
├── run_pipeline.ipynb
├── requirements.txt
└── README.md
```

## Software
|               Software              |      Version(s)     |                       Resource URL                       |       RRID      |                                               Notes                                               |
|:-----------------------------------:|:-------------------:|:--------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------:|
|               Python Programming Language              |      3.9.16     |        http://www.python.org/        | RRID:SCR_008394 |                Refer to requirements.txt for necessary packages                |
