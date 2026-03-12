# CNV-Finder: Alternate Version
 
## Overview
CNV-Finder is a novel pipeline integrating a Long Short-Term Memory network on SNP array data to expedite large-scale identification of CNVs within predefined genomic regions. Manuscript coming soon!

## Input Requirements

This pipeline requires **signal intensity files** containing two key features from genotyping arrays:

- **Log R Ratio (LRR)**
- **B Allele Frequency (BAF)**

We refer to these as **SNP (Single Nucleotide Polymorphism) metrics files**.

### Preparing Your Data

SNP metrics can be extracted from microarray metadata files such as Illumina's IDAT format. If you're working with IDAT files, see the [SNP Metrics repository](https://github.com/nvk23/SNP_metrics/tree/main) for instructions on generating these files using Illumina's IAAP-CLI tool.

### Optimizing Performance

For faster processing, we recommend using the Parquet format. See `example_data/snp_metrics` for reference files, and check out `run_pipeline.ipynb` for examples of how to load and inspect your data.

## Now let's begin CNV hunting!

### Clone the repository:

````
git clone https://github.com/nvk23/CNV-Finder.git

cd CNV-Finder
````

### [Optional] Create a Conda environment for Python 3.8+:

````
conda create -n "cnv_finder" python=3.11

conda activate cnv_finder
````

### Install the required packages:

````
pip install -r requirements.txt
````

## Running the Pipeline
Open the `run_pipeline.ipynb` notebook and sequentially run through each cell to perform the 3 major processes in the pipeline: ML data preprocessing, application of pre-trained models/training new models on the prepared data, and the creation of app-ready files for visualizations in the CNV-Finder app.

### ONNX Model Support

CNV-Finder uses [ONNX (Open Neural Network Exchange)](https://onnx.ai/) format for models, providing cross-platform compatibility and removing specific Python dependencies. Although original models were trained on Python 3.9.16, all pre-trained models are now compatible with Python 3.8+ through [ONNX Runtime](https://onnxruntime.ai/) and can run on various platforms.

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
│   │   ├── keras/  
│   │   │   ├── final_del_5_50_combo4_lstm.keras
│   │   │   ├── final_dup_10_70_combo6_lstm.keras
│   │   │   ├── updated_del_5_50_combo4_lstm.keras
│   │   │   ├── updated_dup_10_70_combo6_lstm.keras
│   │   │   ├── prelim_del_5_50_combo4_lstm.keras
│   │   │   ├── prelim_dup_10_70_combo6_lstm.keras
│   │   ├── ONNX/  
│   │   │   ├── final_del_5_50_combo4_lstm.onnx
│   │   │   ├── final_dup_10_70_combo6_lstm.onnx
│   │   │   ├── updated_del_5_50_combo4_lstm.onnx
│   │   │   ├── updated_dup_10_70_combo6_lstm.onnx
│   │   │   ├── prelim_del_5_50_combo4_lstm.onnx
│   │   │   ├── prelim_dup_10_70_combo6_lstm.onnx
│   ├── exons/
│   │   ├── PARK2_exons.csv
│   ├── NBA_metadata/
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
|               Python Programming Language              |      3.8+     |        http://www.python.org/        | RRID:SCR_008394 |                Refer to requirements.txt for necessary packages                |