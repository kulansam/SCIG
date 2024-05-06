#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import qnorm
import regex as re
import os
from sklearn.linear_model import LogisticRegression
from rpy2.robjects.packages import importr
# import R's "base" package
base = importr('base')
import rpy2
# import R's "utils" package
utils = importr('utils')
from rpy2.robjects.vectors import FloatVector
from sklearn.preprocessing import StandardScaler
from bioinfokit.analys import norm, get_data
import anndata
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")



tpm = norm()
def pesudobulk_cigpred_human (inputtype_pred,exp_filename_pred,alll_seq_features1,gene_name_length1,training_table_human):

    gene_exp0=exp_filename_pred
    gene_name_map=gene_name_length1
    alll_seq_features=alll_seq_features1  
    inputtype_pred=inputtype_pred
    training_table_human=training_table_human
    gene_exp0.index=gene_exp0['Genename']
    del gene_exp0['Genename']
    gene_exp0=gene_exp0.loc[~(gene_exp0==0).all(axis=1)]## REMOVE GENES that zero count across all the cell types/tissues
    
    gene_name_map=gene_name_map.drop_duplicates(subset=('Geneid'),keep='first')
    gene_exp_lengths = pd.merge(gene_name_map,gene_exp0, on=["Genename"], how="inner")
    gene_exp_lengths.index=gene_exp_lengths['Geneid']
    gene_exp_lengths=gene_exp_lengths.drop(['Geneid','Genename'], axis=1)
    gene_exp_lengths=gene_exp_lengths.fillna(0)
    
    if inputtype_pred =='rawcount':
        print("Un-normalised expression rawcounts are processing")
        #TPM (Transcripts per million) Calculation    
        tpm.tpm(df=gene_exp_lengths, gl='length')
        tpm_df_gene = tpm.tpm_norm
        tpm_df_gene['Geneid'] = tpm_df_gene.index
        tpm_df_gene1 = tpm_df_gene.reset_index(drop=True)
        tpm_df_gene1.index=tpm_df_gene1['Geneid']
    
        gene_exp=tpm_df_gene1.drop_duplicates()
        del gene_exp['Geneid']
        gene_exp=qnorm.quantile_normalize(gene_exp,axis=1)
        gene_exp = np.log2(gene_exp)
        gene_exp = gene_exp.replace({-np.inf: gene_exp[np.isfinite(gene_exp)].min().min()})
    
    else:
        print("SKIPPING the expression normalization step, Since the input expression values are already in TPM")
        gene_exp=gene_exp_lengths.drop_duplicates()
        del gene_exp['length']
        gene_exp=qnorm.quantile_normalize(gene_exp,axis=1)
        gene_exp = np.log2(gene_exp)
        gene_exp = gene_exp.replace({-np.inf: gene_exp[np.isfinite(gene_exp)].min().min()})

            
    cig_score_matrix = pd.DataFrame()
    ##ITERATE the loop for predicting CIGs in each cell type
    ten_celltype_train=gene_exp.columns.to_list()
    for i in ten_celltype_train:
        print("SCIG: Identifying the CIGs in ",i)
        ##GETTING GENE EXP DATA
        i1=i
        req_gene_exp=gene_exp.loc[:, gene_exp.columns.isin([i1])]
        req_gene_exp['Geneid']=req_gene_exp.index
        req_gene_exp=req_gene_exp.drop_duplicates()
        req_gene_exp.columns=['TPM_exon_exp_mean','Geneid']
        req_gene_exp.reset_index(drop=True, inplace=True)   
        #combine all seq features with gene expresion features
        combine_feat=pd.merge(req_gene_exp,alll_seq_features,on=['Geneid'],how='inner')
        combine_feat = combine_feat[['Geneid','TPM_exon_exp_mean','hpa_256_tiss_ntpmtau', 'gtex_db_tissuegini',
           'Total_utr_3_len', 'tss_dis', 'TF_motif_total_exons',
           'rna_motif_total_gene_norm_len_cv', 'utr3_mirna_count',
           'exon_cs_mean', 'tss_500bp_cs_median', 'V_codon_bias_sum',
           'P_codon_bias_max', 'Y_codon_bias_max', 'A_codon_bias_median',
           'cds_len', 'AT_cds_pct', 'AT_intron_pct', 'AT_tss_3kb_PTG_pct',
           'GC_tss_4kb_PTG_pct']]
        combine_feat=combine_feat.drop_duplicates()       
        combine_feat=combine_feat.set_index('Geneid')
        combine_feat=combine_feat.fillna(0)
        combine_feat.replace([np.inf, -np.inf], 0, inplace=True)
        ##cig_pred trained model  human
        X_train = training_table_human
        X_train = X_train.drop(columns=['label', 'Gene_name'])
        Y_train = training_table_human['label']
        cig_scaler = StandardScaler()
        X_train1 = cig_scaler.fit_transform(X_train)
        X_train = pd.DataFrame(X_train1, index=X_train.index, columns=X_train.columns)
        cig_loaded_model = LogisticRegression(penalty='l1',solver='liblinear')
        cig_loaded_model.fit(X_train,Y_train)
    
        pr = cig_scaler.fit_transform(combine_feat)
        pr = pd.DataFrame(pr, index=combine_feat.index, columns=combine_feat.columns)
        pred_cols = list(pr.columns.values)
        pred = pd.Series(cig_loaded_model.predict_proba(pr[pred_cols])[:,1])
        desc_fun_pred = pd.Series(cig_loaded_model.decision_function(pr[pred_cols]))
        desc_fun_pred = pd.DataFrame(desc_fun_pred)    
        desc_fun_pred.columns=[i1+'_distance']
        pred = pd.DataFrame(pred)
        pred.columns=[i1]
        gene_ids = combine_feat.index
        gene_ids = pd.DataFrame(gene_ids)
        cig_score=pd.concat([gene_ids, pred], axis=1)        
        cig_score = pd.concat([cig_score, desc_fun_pred], axis=1)
        stats = importr('stats')
        from scipy.stats import norm
        for l_mean in np.arange(-20, 1.15, 0.25):
            cig_score[i1+'_p_value'] = 1 - norm.cdf(cig_score[i1+'_distance'], loc=l_mean)
            cig_score[i1+'_FDR'] = stats.p_adjust(FloatVector(cig_score[i1+'_p_value'].tolist()),method='BH')
            if cig_score[cig_score[i1+'_FDR'] < 0.05].shape[0] < 600:
                break
        cig_score=cig_score.sort_values(by=[i1+'_distance'], ascending=[False])
        cig_score[i1+'_rank']=range(1, cig_score.shape[0] + 1)
        cig_score_matrix=pd.concat([cig_score_matrix, cig_score], axis=1)
    cig_score_matrix = cig_score_matrix.loc[:,~cig_score_matrix.T.duplicated(keep='first')]
    cig_score_matrix=pd.merge(cig_score_matrix,gene_name_map,on=['Geneid'],how='inner')
    cig_score_matrix = cig_score_matrix.rename(columns={'Genename': 'symbol'}) 
    del cig_score_matrix['length']
    gene_exp = gene_exp.add_suffix('_exp')
    cig_score_gene_exp=pd.merge(cig_score_matrix,gene_exp,on=['Geneid'],how='inner')
    cig_score_gene_exp=cig_score_gene_exp.round(3)
    return cig_score_gene_exp

def pesudobulk_cig_reg_pred_human (cig_pred_result,all_db_grn_human,tf_human,training_table_cigreg_human,celltype_names):
    all_db_grn=all_db_grn_human
    cig_score_gene_exp=cig_pred_result
    tf_human.columns=['Geneid','TF']
    celltype=celltype_names
    celltype.remove('Genename')
    n_flag=1
    mas_tf_pred_final=pd.DataFrame()
    for c_type in celltype:
        print("SCIGNet: Identifying the  Master TF of CIGs in ",c_type)
        c_type1=c_type
        cig_dist=c_type+'_distance'
        c_type=c_type+'_exp'
        
        #load master TF data
        mas_tf_cig_exp=pd.DataFrame()
        mas_tf_cig_exp_list=[]
        master_tf_cigscore=pd.merge(cig_score_gene_exp,tf_human,left_on=['symbol'],right_on=['TF'],how='inner')
        master_tf_cigscore['celltype']=c_type
        master_tf_cigscore=master_tf_cigscore[['TF','celltype',c_type1,c_type,cig_dist]]
        master_tf_cigscore.columns=['TF','celltype','Mas_CIG', 'Mas_EXP','Mas_Dist']

        #master_tf_child_edge=pd.merge(master_tf_cigscore,all_db_grn,left_on=['TF'],right_on=['source_genename'],how='inner').groupby(['TF']).size()
        master_tf_child_edge=pd.merge(master_tf_cigscore,all_db_grn,left_on=['TF'],right_on=['source_genename'],how='inner').groupby(['TF']).size()
        master_tf_parent_edge=pd.merge(master_tf_cigscore,all_db_grn,left_on=['TF'],right_on=['target_genename'],how='inner').groupby(['TF']).size()
        #get number of edges from GRN
        master_tf_child_parent_edge=pd.concat([master_tf_parent_edge,master_tf_child_edge], axis=1)
        master_tf_child_parent_edge=master_tf_child_parent_edge.fillna(0)
        master_tf_child_parent_edge.columns=['no_of_parent_edges','no_of_children_edges']
        
        master_tf_parent_edge=pd.merge(tf_human,all_db_grn,left_on=['TF'],right_on=['target_genename'],how='inner')
        master_tf_parent_edge=master_tf_parent_edge[['TF','source_genename']]
        master_tf_parent_edge=master_tf_parent_edge.drop_duplicates()
        master_tf_cigscore_parent=pd.merge(cig_score_gene_exp,master_tf_parent_edge,left_on=['symbol'],right_on=['source_genename'],how='inner')
        master_tf_cigscore_parent=master_tf_cigscore_parent[['TF',c_type1,c_type,cig_dist]]
        
        master_tf_parent_CIG_median=master_tf_cigscore_parent[['TF',c_type1,c_type,cig_dist]].groupby(['TF']).median()
        c_type11=str(c_type1)+'-CIG_median'    
        c_type112=str(c_type)+'-EXP_median'
        c_type113=str(cig_dist)+'-Dist_median'    
        master_tf_parent_CIG_median.columns=[c_type11,c_type112,c_type113]
        master_tf_parent_CIG_mean=master_tf_cigscore_parent[['TF',c_type1,c_type,cig_dist]].groupby(['TF']).mean()
        c_type11=str(c_type1)+'-CIG_mean'    
        c_type112=str(c_type)+'-EXP_mean'
        c_type113=str(cig_dist)+'-Dist_mean'    

        master_tf_parent_CIG_mean.columns=[c_type11,c_type112,c_type113]
        master_tf_parent_CIG_cv=(master_tf_cigscore_parent[['TF',c_type1,c_type,cig_dist]].groupby(['TF']).std())/(master_tf_cigscore_parent[['TF',c_type,c_type1]].groupby(['TF']).mean())   
        c_type11=str(c_type1)+'-CIG_CV'    
        c_type112=str(c_type)+'-EXP_CV'
        c_type113=str(cig_dist)+'-Dist_CV' 
        master_tf_parent_CIG_cv.columns=[c_type11,c_type112,c_type113]
        
        master_tf_child_edge=pd.merge(tf_human,all_db_grn,left_on=['TF'],right_on=['source_genename'],how='inner')
        master_tf_child_edge=master_tf_child_edge[['TF','target_genename']]
        master_tf_child_edge=master_tf_child_edge.drop_duplicates()
        master_tf_cigscore_child=pd.merge(cig_score_gene_exp,master_tf_child_edge,left_on=['symbol'],right_on=['target_genename'],how='inner')
        
        master_tf_child_CIG_mean=master_tf_cigscore_child[['TF',c_type1,c_type,cig_dist]].groupby(['TF']).mean()
        c_type11=str(c_type1)+'-CIG_mean'    
        c_type112=str(c_type)+'-EXP_mean'  
        c_type113=str(cig_dist)+'-Dist_mean'    

        master_tf_child_CIG_mean.columns=[c_type11,c_type112,c_type113]
        
        master_tf_child_CIG_median=master_tf_cigscore_child[['TF',c_type1,c_type,cig_dist]].groupby(['TF']).median()
        c_type11=str(c_type1)+'-CIG_median'    
        c_type112=str(c_type)+'-EXPmedian'    
        c_type113=str(cig_dist)+'-Distmedian' 
        master_tf_child_CIG_median.columns=[c_type11,c_type112,c_type113]
        
        master_tf_child_CIG_cv=(master_tf_cigscore_child[['TF',c_type1,c_type,cig_dist]].groupby(['TF']).std())/(master_tf_cigscore_child[['TF',c_type,c_type1]].groupby(['TF']).mean())
        c_type11=str(c_type1)+'-CIG_CV'    
        c_type112=str(c_type)+'-EXP_CV'   
        c_type113=str(cig_dist)+'-Dist_CV' 
        master_tf_child_CIG_cv.columns=[c_type11,c_type112,c_type113]

        master_tf_cig_final=pd.merge(tf_human,master_tf_child_parent_edge,on=['TF'],how='left')
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_cigscore,on=['TF'],how='left')    
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_parent_CIG_mean,on=['TF'],how='left')
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_parent_CIG_median,on=['TF'],how='left')
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_parent_CIG_cv,on=['TF'],how='left')
        
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_child_CIG_mean,on=['TF'],how='left')
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_child_CIG_median,on=['TF'],how='left')
        master_tf_cig_final=pd.merge(master_tf_cig_final,master_tf_child_CIG_cv,on=['TF'],how='left')
        master_tf_cig_final=master_tf_cig_final.fillna(0)
        del master_tf_cig_final['celltype']
        master_tf_cig_final.columns=['Geneid','TF','no_of_parent_edges', 'no_of_children_edges','CIG', 'EXP','CIG_distance','CIGscore_parent_mean','Exp_parent_mean','CIG_distance_parent_mean','CIGscore_parent_median','Exp_parent_median','CIG_distance_parent_median','CIGscore_parent_CV','Exp_parent_CV','CIG_parent_distance_CV',
                                     'CIGscore_child_mean','Exp_child_mean','CIG_distance_child_mean','CIGscore_child_median','Exp_child_median','CIG_distance_child_median','CIGscore_child_CV','Exp_child_CV','CIG_distance_child_CV']         

        
        ##cig_reg-pred trained model  human
        X_train = training_table_cigreg_human
        X_train = X_train.drop(columns=['label'])
        Y_train = training_table_cigreg_human['label']
        mastf_scaler = StandardScaler()
        X_train1 = mastf_scaler.fit_transform(X_train)
        X_train = pd.DataFrame(X_train1, index=X_train.index, columns=X_train.columns)
        mastf_loaded_model = LogisticRegression(penalty='l1',solver='liblinear')
        mastf_loaded_model.fit(X_train,Y_train)

        req_feat_saved_model=mastf_loaded_model.feature_names_in_
        master_tf_cig_final=master_tf_cig_final.set_index('TF')
        del master_tf_cig_final['Geneid']
        master_tf_cig_final=master_tf_cig_final[req_feat_saved_model]
        master_tf_cig_final.replace([np.inf, -np.inf], 0, inplace=True)

        pr = mastf_scaler.fit_transform(master_tf_cig_final)
        pr = pd.DataFrame(pr, index=master_tf_cig_final.index, columns=master_tf_cig_final.columns)
        pred_cols = list(pr.columns.values)
        pred = pd.Series(mastf_loaded_model.predict_proba(pr[pred_cols])[:,1])
        pred = pd.DataFrame(pred)
        pred.columns=[c_type]
        gene_ids = master_tf_cig_final.index
        gene_ids = pd.DataFrame(gene_ids)
        mf_pred_f=pd.concat([gene_ids, pred], axis=1)
        
        mastf_desc_fun_pred = pd.Series(mastf_loaded_model.decision_function(pr[pred_cols]))
        mastf_desc_fun_pred = pd.DataFrame(mastf_desc_fun_pred)    
        mastf_desc_fun_pred.columns=[c_type+'_distance']
        mf_pred_f = pd.concat([mf_pred_f, mastf_desc_fun_pred], axis=1)
        stats = importr('stats')
        from scipy.stats import norm
        for l_mean in np.arange(-2, 10, 0.25):
            mf_pred_f[c_type+'_p_value'] = 1 - norm.cdf(mf_pred_f[c_type+'_distance'],loc=l_mean)
            mf_pred_f[c_type+'_FDR'] = stats.p_adjust(FloatVector(mf_pred_f[c_type+'_p_value'].tolist()),method='BH')
            if mf_pred_f[mf_pred_f[c_type+'_FDR'] < 0.05].shape[0]*1 <50:
                break
        mf_pred_f=mf_pred_f.sort_values(by=[c_type+'_distance'], ascending=[False])
        mf_pred_f[c_type+'_rank']=range(1, mf_pred_f.shape[0] + 1)
        if n_flag==1:
            
            mas_tf_pred_final=pd.concat([mas_tf_pred_final,mf_pred_f], axis=1)
        else:
            
            del mf_pred_f['TF']
            mas_tf_pred_final=pd.concat([mas_tf_pred_final,mf_pred_f], axis=1)
        n_flag=n_flag+1
    mas_tf_pred_final.index=mas_tf_pred_final['TF']
    mas_tf_pred_final=mas_tf_pred_final.round(3)
    return mas_tf_pred_final
    
def cig_pred_singlecell_human(adata,features,all_seq_features_human,training_table_human,exp_filename_pred):   
    adata=adata
    cell_ranger_file_name=exp_filename_pred
    all_seq_features1=all_seq_features_human    
    ##PROTEIN_CODING GENES
    gene_length = pd.read_csv("./Generic_features_GTF_gene_with_features_human.txt",sep="\t") 
    gene_length =gene_length[['Geneid', 'length']]
    gene_length=gene_length.drop_duplicates(subset='Geneid', keep="first")
    cc=adata[:,[g in gene_length['Geneid'].tolist() for g in adata.var['Geneid']]]
    cc.var_names_make_unique
    ##BASIC FILTEERING BASED GENE EXP AND CELLS
    #sc.pp.filter_cells(cc, min_genes=200)
    
    req_prot_cod_genes=cc.to_df()
    req_prot_cod_genes=req_prot_cod_genes.T
    req_prot_cod_genes.index.names = ['Geneid']
    
    #TPM (Transcripts per million)
    gene_length.index=gene_length['Geneid']
    del gene_length['Geneid']
    req_prot_cod_genes = pd.merge(req_prot_cod_genes, gene_length, on=["Geneid"], how="inner")
    from bioinfokit.analys import norm, get_data
    nm = norm()
    nm.tpm(df=req_prot_cod_genes, gl='length')
    req_prot_cod_genes = nm.tpm_norm
    req_prot_cod_genes.head(2)
    
    ##REQ features
    combine_feat=all_seq_features1[['Geneid', 'hpa_256_tiss_ntpmtau', 'gtex_db_tissuegini',
       'Total_utr_3_len', 'tss_dis', 'TF_motif_total_exons',
       'rna_motif_total_gene_norm_len_cv', 'utr3_mirna_count',
       'exon_cs_mean', 'tss_500bp_cs_median', 'V_codon_bias_sum',
       'P_codon_bias_max', 'Y_codon_bias_max', 'A_codon_bias_median',
       'cds_len', 'AT_cds_pct', 'AT_intron_pct', 'AT_tss_3kb_PTG_pct',
       'GC_tss_4kb_PTG_pct']]
    
    ##cig_pred trained model  human
    X_train = training_table_human
    X_train = X_train.drop(columns=['label', 'Gene_name'])
    Y_train = training_table_human['label']
    cig_scaler = StandardScaler()
    X_train1 = cig_scaler.fit_transform(X_train)
    X_train = pd.DataFrame(X_train1, index=X_train.index, columns=X_train.columns)
    cig_loaded_model = LogisticRegression(penalty='l1',solver='liblinear')
    cig_loaded_model.fit(X_train,Y_train)
    
    pro_cds_gens_names_m=features[['Geneid','gene_symbols']]
    ###ADD EACH CELL GENE EXP TO PREDICT CIG SCORE
    cig_score_matrix1 = pd.DataFrame()
    cig_score_matrix = pd.DataFrame()
    desc_fun_pred_cig_list= pd.DataFrame()
    desc_fun_pred_cig=[]
    cells_list=list(req_prot_cod_genes.columns)
    cells_list=cells_list[-5:]
    n_c=1
    for i in cells_list:
        print("CIG-Pred: Identifying the CIGs in ",i,"; ", "Processed:", n_c, "out of ",len(cells_list), " cells")
        gene_exp_each_cell=req_prot_cod_genes[[i]]
        gene_exp_each_cell.index=req_prot_cod_genes.index
        gene_exp_each_cell = np.log2(gene_exp_each_cell)
        gene_exp_each_cell.columns=['TPM_exon_exp_mean']
        test_file_model=pd.merge(gene_exp_each_cell,combine_feat, on='Geneid')
        test_file_model=test_file_model[['Geneid','TPM_exon_exp_mean', 'hpa_256_tiss_ntpmtau', 'gtex_db_tissuegini',
            'Total_utr_3_len', 'tss_dis', 'TF_motif_total_exons',
            'rna_motif_total_gene_norm_len_cv', 'utr3_mirna_count',
            'exon_cs_mean', 'tss_500bp_cs_median', 'V_codon_bias_sum',
            'P_codon_bias_max', 'Y_codon_bias_max', 'A_codon_bias_median',
            'cds_len', 'AT_cds_pct', 'AT_intron_pct', 'AT_tss_3kb_PTG_pct',
            'GC_tss_4kb_PTG_pct']]
        
        test_file_model.replace([np.inf, -np.inf], 0, inplace=True)
        test_file_model=test_file_model.fillna(0)
        test_file_model=test_file_model.set_index('Geneid')
    
        pr = cig_scaler.fit_transform(test_file_model)
        pr = pd.DataFrame(pr, index=test_file_model.index, columns=test_file_model.columns)
        pr_q=pr
        pr_q['Geneid'] = pr_q.index
        
        pr_q.reset_index(drop=True, inplace=True)    
        pr_q=pr
        pr_q['Geneid'] = pr_q.index
        del pr_q['Geneid']
        
        pred_cols = list(pr.columns.values)
        # apply the whole pipeline to data
        pred = pd.Series(cig_loaded_model.predict_proba(pr[pred_cols])[:,1])
        desc_fun_pred = pd.Series(cig_loaded_model.decision_function(pr[pred_cols]))
        desc_fun_pred = pd.DataFrame(desc_fun_pred) 
        pred = pd.DataFrame(pred)
        pred.columns =[i]
        desc_fun_pred.columns =[i]
        desc_fun_pred=desc_fun_pred.round(2)
        #cig_score_matrix1=pd.concat([cig_score_matrix, pred], axis=1)#for cig prob score
        cig_score_matrix=pd.concat([cig_score_matrix, desc_fun_pred], axis=1)##for distance
        desc_fun_pred_cig_list=pd.concat([cig_score_matrix, desc_fun_pred], axis=1)
        desc_fun_pred.index=test_file_model.index
        desc_fun_pred=pd.merge( pro_cds_gens_names_m,desc_fun_pred, on=["Geneid"], how="inner")
        test_file_model1=pd.merge( desc_fun_pred,test_file_model, on=["Geneid"], how="inner")
        i11=str(i)+'_CIGscore'
        test_file_model1[i11] = np.where(test_file_model1['TPM_exon_exp_mean'] == 0,  test_file_model1['TPM_exon_exp_mean'],test_file_model1[i])
        desc_fun_pred=test_file_model1[[i11,'Geneid', 'gene_symbols']]
        desc_fun_pred.columns=[i,'Geneid', 'gene_symbols']
        desc_fun_pred1=test_file_model1[['Geneid',i11]]
        #cig_score_matrix1=pd.concat([cig_score_matrix1, desc_fun_pred1], axis=1)##for distance
        if n_c==1:
            cig_score_matrix1=pd.concat([cig_score_matrix1, desc_fun_pred1], axis=1)
        else:
            cig_score_matrix1=pd.merge(cig_score_matrix1, desc_fun_pred1, on=['Geneid'])
    
        print(n_c)    
        n_c+=1
        
        desc_fun_pred=desc_fun_pred.sort_values(by=[i],ascending=False)
        desc_fun_pred=desc_fun_pred.head(500)
        hgs_cig=desc_fun_pred['gene_symbols'].tolist()
        
        desc_fun_pred_cig.extend(hgs_cig)
    
    desc_fun_pred_cig_list.index=test_file_model.index
    #ADD CELL BARCODES AND GENES NAMES
    #pro_cds_gens_names_m=features[['Geneid','gene_symbols']]
    cig_score_matrix=cig_score_matrix1
    cig_score_matrix=pd.merge( pro_cds_gens_names_m,cig_score_matrix, on=["Geneid"], how="inner")
    del cig_score_matrix['Geneid']
    cig_score_matrix=cig_score_matrix.set_index('gene_symbols')
    #cig_score_matrix.to_csv(r'./cigscore_matrix_SRR10590805.csv',index = False,sep=",")
    v=pd.DataFrame(cig_score_matrix.columns)
    v.to_csv(str(cell_ranger_file_name)+'cig_barcodes.tsv',sep="\t",index=False,header=False)
    v=pd.DataFrame(cig_score_matrix.index)
    v.to_csv(str(cell_ranger_file_name)+'cig_genes.tsv',sep="\t",index=False,header=False)
    ##CONVERT INTO MATRIX FORM
    adatass = sc.AnnData(cig_score_matrix)
    return adatass
    #adatass.T.write_h5ad("cig_matrix.h5ad")
