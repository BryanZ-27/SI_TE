# SI_TE

## Repository Overview
This repository, `SI_TE`, contains the code and data used to support the paper "Translational Fidelity and Longevity are genetically linked". These resources aim to facilitate other researchers in reproducing the analyses and results presented in the paper, thereby promoting research progress in the relevant field.

## Repository Contents
### Code Files
- **R Scripts**: Multiple files with the `.R` extension, each responsible for different analysis tasks.
    - `1.model_simulation.R`: Used to simulate the models involved in the paper. By simulating different parameters and conditions, it validates relevant theories and hypotheses.
    - `2.Fig1_FigS1.R`, `4.Fig2_FigS2_FigS3.R`, etc.: These scripts are used to generate the figures in the paper.
    - `3.SI_QC.R`, `5.TE_QC.R`: Used for quality control of raw data.
    - `3_1.SRver5-1_function.R`, `3_2.SI_biorep_clean_function.R`, etc.: These scripts containing functions used in quality control of raw data.

### Data Files
- **All data is stored in the `data.zip` file.**
- **`.Rdata` Files**: Such as `2.cor_filter_P_data_200.Rdata`, `3.raw_survival_data.Rdata`, `8.chr11_2_genes.Rdata`, etc. These files store raw data, intermediate results or specific datasets from the analysis.
- **`.txt` Files**:
    - `8.lifespan_data.txt`: Contains the lifespan data of multiple individuals. Each line records the ID of an individual and the corresponding lifespan value, providing basic data for GWAS.
    - `8.sex_data.txt`: Stores the gender information of individuals. Each line contains the ID of an individual and the corresponding gender identifier, which is used for genotype-conditional association test.

## Code Usage Instructions
### Environment Requirements
Ensure that you have an R language environment installed on your system and have installed the relevant R packages required by the scripts. You can use the following command to install common R packages:
```R
install.packages(c("package1", "package2")) # Replace with the actual package names
```

### Running Steps
1. **Download the Repository**: Clone this repository to your local machine:
```bash
git clone https://github.com/BryanZ-27/SI_TE.git
```
2. **Run the Scripts**: Open the R environment and set the working directory to the directory where the repository is located:
```R
setwd("path/to/SI_TE")
```
Then run the corresponding R scripts as needed. For example, to run the model simulation script:
```R
source("1.model_simulation.R")
```
To run the figure generation script:
```R
source("2.Fig1_FigS1.R")
```

## Data Description
### `.Rdata` Files
These files can be loaded in the R environment using the `load()` function. For example:
```R
load("2.cor_filter_P_data_200.Rdata")
```

### `.txt` Files
You can use the `read.table()` or `read.csv()` functions to read them. For example, to read the lifespan data:
```R
lifespan_data <- read.table("8.lifespan_data.txt", header = TRUE)
```

## Contribution and Feedback
If you encounter any issues or have suggestions for improvement during use, please feel free to submit an issue or a pull request.

## Citation Information
If you use the code and data from this repository, please cite the paper "Translational Fidelity and Longevity are genetically linked".