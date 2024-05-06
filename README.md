# Welcome to SCIG: A tool for Specifying the Cell Identity Genes (CIGs) and their master transcription factors at single-cell levels, utilizing genetic sequence and RNA expression profiles
In every cell type, a unique set of cell identity genes (CIGs) plays a pivotal role in defining its specific characteristics. Alongside other gene categories like housekeeping genes and heat shock genes, cells express their CIGs, crucial for guiding cellular differentiation and the formation of tissues and organs in organisms. Thus, we have developed two logistic regression-based machine-learning methodologies for the identification of cell identity genes (CIGs) and their master transcription factors:

1. SCIG: SCIG: Specifying  the cell identity genes in a given cell (s) using RNA expression and genetic sequence profiles.
2. SCIGNet: Identifying the master transcription factors of cell identity genes network  (SCIG prediction result required).
   
# What constitutes the input data for SCIG and CIG-reg-Pred algorithms?

The input for SCIG consists of either raw read counts from bulk-seq or gene per cell unique molecular identifier (UMI) count matrix from single-cell sequencing data. In addition, SCIG automatically utilizes precomputed genetic sequence features during its prediction process.

For CIG-reg-Pred, the required inputs include the predicted cell identity scores obtained from the SCIG algorithm, along with gene regulatory network (GRN) information. It's worth noting that the GRN information is already integrated into the model.
# Installation
- Download the SCIG package from GitHub:
  ```sh
  git clone https://github.com/kulansam/SCIG.git
  ```
  ```sh
  cd SCIG
  ```
- Extract the required files in the data directory:
  ```sh
  cd data
  tar -xf human_mouse_GRN.zip
  cd ..
  ```
- Create a conda environment with Python:
  ```sh
  conda create -n scig (Python version >=3.9 preferable)
  ```
    ```sh
  conda activate scig
  ```
- Install required all Python packages
  ```sh
  pip install -r requirements.txt
  ```
-  install SCIG
  ```sh
  pip install .
  ```
# Tutorial 
- How to use SCIG and SCIGNet to uncover the cell identity genes (CIGs) and their master transcription factors using either bulk-RNA  or single-cell RNA sequencing profiles?

  Use the following command:
  
  ```sh
    python SCIG.py -organism <hs | mm> -assaytype <bulk | single> -inputtype <rawcount | tpm | umicount> -file <tab separated expression data file | cellranger output folder>
  ```
  ```sh
    For more information:
  
    -organism: Name of the organism.

                1. For humans, the value should be "hs".
                2. For mice, the value should be "mm".

    -assaytype: The input RNA expression data was quantified differently based on the level of analysis.

                1. For bulk or pseudo bulk data, the value is "bulk".
                2. For Single-cell level data, the value is "single".

    -inputtype: The format of the input data differs based on the analysis level.
                1. For bulk or pseudobulk data, the value is either "rawcount" or "tpm".
                2. For single-cell data, the value should be "umicount".
  
    -file: name of the input file name with the tab-separated file format.
                1. For bulk or pseudobulk data, the tab-separated file should contain the 'Genename' as the first column name and followed by expression values of cell type (s). Example: Genename<tab>celltype1<tab>..celltypen)
                2. For single-cell data, the directory path for the Cell Ranger output folder. It should contain the following files: barcodes.tsv, features.tsv, and matrix.mtx.

# Cite us
  ```sh
    Please cite our paper 'Kulandaisamy Arulsamy, Bo Xia,  Lili Zhang, Hong Chen & Kaifu Chen (2024). Machine Learning Uncovers Cell Identity Genes in Single Cells by Genetic Sequence Features
  ```
   ```sh
# Term of Usage
By accessing SCIG data, you agree to the following terms:

1. You agree not to share the SCIG data, whether modified or unmodified, with individuals outside your research group. This includes preventing access by unauthorized individuals and refraining from directly providing the data to others.

2. You agree not to develop another website or methods using the SCIG data without prior permission. Contact us for any such intentions.

3. You agree to appropriately cite the SCIG paper and its original contributions if utilized in your work.

4. You certify that you are authorized to accept these terms on behalf of your institution.

# Contact
For any queries, contact  kulandai28@gmail.com, or Kaifu.Chen@childrens.harvard.edu. Copy right @ Dr.Kaifu Chen lab, @ Boston Children's Hospital, Harvard Medical School
