# CNV-Finder

[![DOI](https://zenodo.org/badge/758269891.svg)](https://doi.org/10.5281/zenodo.19009011)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
 
## Overview
CNV-Finder is a novel pipeline integrating a Long Short-Term Memory network on SNP array data to expedite large-scale identification of CNVs within predefined genomic regions. Check out our [preprint](https://www.biorxiv.org/content/10.1101/2024.11.22.624040v2) for more project information!

## Input Requirements

This pipeline requires **signal intensity files** containing two key features from genotyping arrays:

- **Log R Ratio (LRR)**
- **B Allele Frequency (BAF)**

We refer to these as **SNP (Single Nucleotide Polymorphism) metrics files**.

### Preparing Your Data

SNP metrics can be extracted from microarray metadata files such as Illumina's IDAT format. If you're working with IDAT files, see the [SNP Metrics repository](https://github.com/nvk23/SNP_metrics/tree/main) for instructions on generating these files using Illumina's IAAP-CLI tool.

### Example Data

Example SNP metrics files are hosted on Hugging Face. You can browse and download them here:

**[CNV-Finder Dataset on Hugging Face](https://huggingface.co/datasets/nkuznet/CNV-Finder/tree/main)**

These files can be used to test the pipeline end-to-end or as a reference for formatting your own input data. For faster processing, we recommend using the Parquet format as shown by these example files. Refer to `run_pipeline.ipynb` for importing, loading, and inspecting this data.

Downloaded data should follow this Hive-partitioned directory structure:

```
example_data/
└── snp_metrics/
    ├── {barcode}/
    │   ├── {barcode}_{sample}/
    │   │   ├── chromosome=1/
    │   │   │   └── *.parquet
    │   │   ├── chromosome=2/
    │   │   ├── ...
    │   │   ├── chromosome=22/
    │   │   ├── chromosome=X/
    │   │   ├── chromosome=Y/
    │   │   └── chromosome=M/
```

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
│   ├── NBA_metadata/
│   │   ├── CHROM=1/
│   │   │   └── *.parquet
│   │   ├── CHROM=2/
│   │   ├── ...
│   │   ├── CHROM=22/
│   │   ├── CHROM=X/
│   │   ├── CHROM=Y/
│   │   └── CHROM=M/
│   ├── exons/
│   │   ├── PARK2_exons.csv
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
