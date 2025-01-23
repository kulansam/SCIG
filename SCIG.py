#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import os
import numpy as np
import anndata
import scanpy as sc
from scipy import io
import warnings
warnings.filterwarnings("ignore")
from scig_human_functions import *
from scig_mouse_functions import *


path = os.getcwd()
#print (path)
os.chdir(path)
def printHelp():
    
    print ("\nError found in your command. Kindly please go through the SCIG tutorial here:\n")
    print ("********************************************")
    print ("******************Welcome!******************")
    print ("********************************************\n")

    print ("SCIG has following two functional modules:\n")
    print ("i. SCIG: Specifying  the cell identity genes in a given cell (s) using RNA expression and genetic sequence profiles\n")
    print ("ii. SCIGNet: Identifying the master transcription factors of cell identity genes network  (SCIG prediction result required)\n")
    
    print ("******************Installation:******************\n")
    print ("To install SCIG and SCIGNet, run the following commands in the terminal/command prompt:\n")
    print ("1. Downloading the SCIG/SCIGNet package: enter 'git clone https://github.com/kulansam/CIGpred.git' \n")
    print ("2. Navigating into the SCIG folder: enter 'cd SCIG' \n")
    print ("3. Creating conda environment for SCIG: enter 'conda create -n scig' (Python>=3.9 preferable)\n")
    print ("4. Activating conda environment: enter 'conda activate scig' \n")
    print ("5. Install required packages: enter 'pip install -r requirements.txt' \n")
    print ("6. Compling the SCIG: enter 'pip install .' \n")

    print ("******************Tutorial:******************\n")

    print ("To run SCIG and SCIGNet, enter the following command under the 'src' directory of SCIG folder (Available to download the Github page)\n")
    print("  python SCIG.py -organism <hs> -assaytype <bulk> -inputtype <rawcount> -file <expression_data.txt>\n")
    print(" Here, \n -organism : Name of the organism.\n 1. For humans, the value should be either 'hs'.\n 2. For mice, the value should be 'mm'. \n")
    print("-assaytype: The input RNA expression data is would quantified either bulk or single-cell level:\n 1.  For bulk or pseudo bulk RNA data, the value is 'bulk'.\n 2. For Single-cell level data, the value is 'single'.\n\n")
    print("-inputtype: The format of the input data differs based on the analysis level.\n\n 1. For bulk or pseudobulk data, the value is either 'rawcount' or 'tpm'.\n 2. For single-cell data, the value should be umicount.\n\n")
    print(" if assaytype == bulk, the  -file: name of the input file name with the tab-separated file format \n(Example:/Users/Kulan/SCIG_input/input_rna_seq.txt, if file path has 'space' please use the backslash).\n")
    print ("The  'input_rna_seq.txt' file should contains the 'Genename' as first column name and followed by expression values of cell type (s). Example: Genename<tab>celltype1<tab>..celltypen \n")
    print(" if assaytype == single, the  -file: The directory path for the Cell Ranger output folder. It should contain the following files: barcodes.tsv, features.tsv, and matrix.mtx.\n\n")
    print ("******************Output:******************\n")

    print ("The output files will be write into user input file directory.\n")
    print("1. SCIG outputs the cell identity gene information in the file name that has combination of the user input file name with 'cig_pred_result.out' extension.\n In case of single-cell RNA seq data, the output file has the extension of '_cig_matrix_out.h5ad' \n")
    print("2. SCIGNet outputs of master transcription factors of cell identity genes in the file name that has combination of user input file name with 'REG_pred_result.out' extension.  \n")



    print ("***************************************************************************")
    print ("******************Thanks for using our SCIG for your research work*********")
    print ("***************************************************************************\n")
    print ("Please cite our paper 'Kulandaisamy Arulsamy, Bo Xia,  Lili Zhang, Hong Chen & Kaifu Chen (2024). Machine Learning Uncovers Cell Identity Genes in Single Cells by Genetic Sequence Features'\n")
    
    print ("For any queries, contact  kulandai28@gmail.com, or Kaifu.Chen@childrens.harvard.edu. Copy right @ Dr.Kaifu Chen lab, @ Boston Children's Hospital, Harvard Medical School\n")
#print(sys.argv)
#print(len(sys.argv))

if (len(sys.argv) ==9):
    parser = argparse.ArgumentParser(usage="...",description='', epilog="Chen lab, Boston Children’s Hospital/Harvard Medical School")
    parser.add_argument('-organism', dest='organism', type=str, help='Mention the organism name of a sample/data')
    parser.add_argument('-assaytype', dest='assaytype', type=str, help='Type of the data either bulk/pseudobulk or single-cell RNA expression data')
    parser.add_argument('-inputtype', dest='inputtype', type=str, help='Mention whether the expression values are in TPM or in raw count')
    parser.add_argument('-file', dest='file', type=str, help='Mention a filename that contains the expression values of genes')


    args = parser.parse_args()
    if args is not None:
        organism_name_pred = args.organism
        assaytype_pred = args.assaytype
        inputtype_pred = args.inputtype
        exp_filename_pred = args.file
        
    ##add sequence features and model training table
    if organism_name_pred =='hs':
        all_seq_features_human = pd.read_csv("../data/Requried_seq_features_cigpred_human.txt",sep="\t")
        gene_name_length_human=pd.read_csv("../data/gene_name_map_with_length_human.txt",sep="\t")
        training_table_human=pd.read_csv("../data/training_table_human.txt",sep="\t")   
        training_table_cigreg_human=pd.read_csv("../data/training_table_master_tf_cigs_human.txt",sep="\t")
        all_db_grn_human=pd.read_csv("../data/GRN_network_compiled_human.txt",sep="\t")
        tf_human=pd.read_csv("../data/tf_human_list.txt",sep="\t")
    elif organism_name_pred =='mm':
        all_seq_features_mouse = pd.read_csv("../data/Requried_seq_features_cigpred_mouse.txt",sep="\t")
        gene_name_length_mouse=pd.read_csv("../data/gene_name_map_with_length_mouse.txt",sep="\t")
        #gene_name_length_mouse['Genename'] = gene_name_length_mouse['Genename'].str.upper()
        training_table_mouse=pd.read_csv("../data/training_table_mouse.txt",sep="\t")
        training_table_cigreg_mouse=pd.read_csv("../data/training_table_master_tf_cigs_mouse.txt",sep="\t")
        all_db_grn_mouse=pd.read_csv("../data/all_db_GRN_combined_mouse.txt",sep="\t")
        tf_mouse=pd.read_csv("../data/tf_mouse_list_geneid.txt",sep="\t")

    
    if organism_name_pred == 'hs' and assaytype_pred == 'bulk' and exp_filename_pred is not None:
        print ("HUMAN",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        exp_matrix=pd.read_csv(exp_filename_pred,sep="\t")
        celltype_names=exp_matrix.columns.to_list()
        #exp_matrix.rename(columns = {str(celltype_names[0]): 'Genename'}, inplace = True)
        #print(exp_matrix)
        cig_pred_output_table=pesudobulk_cigpred_human(inputtype_pred,exp_matrix,all_seq_features_human,gene_name_length_human,training_table_human)
        first_column_symbol = cig_pred_output_table.pop('symbol')
        cig_pred_output_table.insert(0, 'symbol', first_column_symbol)
        cig_pred_output_table.to_csv(str(exp_filename_pred)+'-hs-SCIG_CIGs_result.out', index = False,sep="\t")
        
        cig_reg_pred_output_table=pesudobulk_cig_reg_pred_human(cig_pred_output_table,all_db_grn_human,tf_human,training_table_cigreg_human,celltype_names)
        print(cig_reg_pred_output_table.shape)
        cig_reg_pred_output_table.to_csv(str(exp_filename_pred)+'-hs-SCIG_MasterTFs_result.out', index = False,sep="\t")

    elif organism_name_pred == 'hs' and assaytype_pred == 'single' and inputtype_pred == 'umicount' and exp_filename_pred is not None:
        print ("HUMAN",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        barcodes = pd.read_csv(str(exp_filename_pred)+"barcodes.tsv",sep="\t",header=None)
        barcodes.columns=['barcode']
        features = pd.read_csv(str(exp_filename_pred)+"features.tsv",sep="\t",header=None)
        features.columns=['Geneid','gene_symbols','note']
        counts_mat = io.mmread(str(exp_filename_pred)+"matrix.mtx")
        counts_mat= counts_mat.toarray()
        counts_mat= np.matrix(counts_mat.transpose())
        # create anndata object
        adata = anndata.AnnData(counts_mat,obs=barcodes['barcode'].tolist(),var=features)
        adata.obs.index = barcodes['barcode'].tolist()
        adata.var.index = features['Geneid'].tolist()
        adata.var_names_make_unique()##REMOVE DUPLICATES
        adata_single_cell_cig=cig_pred_singlecell_human(adata,features,all_seq_features_human,training_table_human,exp_filename_pred)
        adata_single_cell_cig.T.write_h5ad(str(exp_filename_pred)+"_hs-cig_matrix_out.h5ad")

    elif organism_name_pred == 'mm' and assaytype_pred == 'bulk' and exp_filename_pred is not None:
        print ("MOUSE",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        exp_matrix=pd.read_csv(exp_filename_pred,sep="\t")
        celltype_names=exp_matrix.columns.to_list()
        print(exp_matrix.shape)
        cig_pred_output_table=pesudobulk_cigpred_mouse(inputtype_pred,exp_matrix,all_seq_features_mouse,gene_name_length_mouse,training_table_mouse)
        first_column_symbol = cig_pred_output_table.pop('symbol')
        cig_pred_output_table.insert(0, 'symbol', first_column_symbol)
        cig_pred_output_table.to_csv(str(exp_filename_pred)+'_mm-SCIG_CIGs_result.out', index = False,sep="\t")
        
        cig_reg_pred_output_table=pesudobulk_cig_reg_pred_mouse(cig_pred_output_table,all_db_grn_mouse,tf_mouse,training_table_cigreg_mouse,celltype_names)
        cig_reg_pred_output_table.to_csv(str(exp_filename_pred)+'-mm-SCIG_MasterTFs_result.out', index = False,sep="\t")

    elif organism_name_pred == 'mm' and assaytype_pred == 'single' and inputtype_pred == 'umicount' and exp_filename_pred is not None:
        print ("MOUSE",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        barcodes = pd.read_csv(str(exp_filename_pred)+"barcodes.tsv",sep="\t",header=None)
        barcodes.columns=['barcode']
        features = pd.read_csv(str(exp_filename_pred)+"features.tsv",sep="\t",header=None)
        features.columns=['Geneid','gene_symbols','note']
        counts_mat = io.mmread(str(exp_filename_pred)+"matrix.mtx")
        counts_mat= counts_mat.toarray()
        counts_mat= np.matrix(counts_mat.transpose())
        # create anndata object
        adata = anndata.AnnData(counts_mat,obs=barcodes['barcode'].tolist(),var=features)
        adata.obs.index = barcodes['barcode'].tolist()
        adata.var.index = features['Geneid'].tolist()
        adata.var_names_make_unique()##REMOVE DUPLICATES
        adata_single_cell_cig=cig_pred_singlecell_mouse(adata,features,all_seq_features_mouse,training_table_mouse,exp_filename_pred)
        adata_single_cell_cig.T.write_h5ad(str(exp_filename_pred)+"-mm_CIG_matrix_out.h5ad")

else:  
    printHelp()
    