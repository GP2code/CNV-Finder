# CNV-Finder: Keras Model Version

## Differences from the 'main' Branch Version
1. _Keras Model Support Only_

This branch includes the code used for preprint/manuscript preparation with the original CNV-Finder modules and Keras models developed in Python 3.9.16.

2. _Simplified Repo_

All explanations walking through CNV-Finder, its app-related features, and using its respective example data and NBA metadata hosted on Hugging Face can be found on the `main branch`!
 
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
conda create -n "cnv_finder" python=3.9.16

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
│   ├── NBA_metadata/
│   ├── custom_intervals.csv
│   ├── glist_hg38_interval.csv
│   └── training_set_IDs.csv
│
├── docs/
│   └── parameter_guide.md
│
├── requirements.txt
└── README.md
```

## Software
|               Software              |      Version(s)     |                       Resource URL                       |       RRID      |                                               Notes                                               |
|:-----------------------------------:|:-------------------:|:--------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------:|
|               Python Programming Language              |      3.9.16     |        http://www.python.org/        | RRID:SCR_008394 |                Refer to requirements.txt for necessary packages                |