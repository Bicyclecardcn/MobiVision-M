#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
import glob
import json
import shutil
import warnings
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from pandas import Series
from anndata import AnnData
from scipy.sparse import issparse
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
sys.path.append(os.path.dirname(__file__))
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor

warnings.filterwarnings("ignore")
os.environ["PYTHONHASHSEED"] = "0"
os.environ["OMP_NUM_THREADS"] = "1"
sc.settings.max_memory = 40


# the input is exp_raw and integrated adata.

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()



class BrowserExpData:
    def __init__(self, adata_bat, exp_data, gene_df, rebuild_form_X, out_prefix_ame, mobilogger:MobiLoggingSystem, sample_name=None):
        self.__exp_data = exp_data
        self.__adata = adata_bat
        self.__rebuild_form_X = rebuild_form_X
        self.__gene_df = pd.read_csv(gene_df, sep="\t", header=None)
        self.__outJsonDir = out_prefix_ame    # outJson + '.zip'
        self.__sample_name = sample_name  # 样本名参数
        self.mobilogger = mobilogger
       # if not os.path.isdir(self.__outJsonDir):
       #     os.makedirs(self.__outJsonDir)

    # def __get_Allgenes(self):
    #     geneName_list = list(self.__adata.var.index)
    #     geneID_list = list(self.__adata.var['gene_ids'])

    #     return geneName_list, geneID_list
        

    def __get_bdata(self):
        if self.__rebuild_form_X:
            adata01 = AnnData(self.__exp_data.X)
            #if os.path.isdir("batch_mtx"):
            #    df = pd.read_csv(os.path.join("batch_mtx", "features.tsv.gz"), header=None, sep='\t')
            #else:
            #    path_prefix = os.getcwd() + "/map_result"
            #    for root, dir, files in os.walk(path_prefix):
            #        if 'filtered' in dir:
            #            feature_path = root + "/filtered/genes.tsv"            
            #            df = pd.read_csv(feature_path, header=None, sep="\t")
            #            break       
            feature_names = self.__gene_df[self.__gene_df.iloc[:,2] == 'Gene Expression'].iloc[:,1].tolist()
            feature_ids = self.__gene_df[self.__gene_df.iloc[:,2] == 'Gene Expression'].iloc[:,0].tolist()
            adata01.var["gene_ids"] = feature_ids
            adata01.var_names = feature_names
            adata01.obs_names = self.__adata.obs_names
            # adata01.raw = self.__exp_data
            adata01.var['n_count'] = adata01.X.sum(axis=0).A1
        else:
            adata01 = self.__adata
        # adata01.var["gene_ids"] = self.__exp_data.var_names
        # gene_exp_list = adata01.var.loc[adata01.var['n_count'] != 0].index
        # bdata = adata01[:,gene_exp_list]
        return adata01


    def __get_NoneZero_genes(self):
        self.__adata.raw = self.__exp_data
        if issparse(self.__adata.X):
            self.__adata.var['n_count'] = self.__adata.X.sum(axis=0).A1
        else:
            self.__adata.var['n_count'] = self.__adata.X.sum(axis=0)
        gene_exp_list = self.__adata.var.loc[self.__adata.var['n_count'] != 0].index
        return gene_exp_list

    def __get_exp_each_dict(self, geneList, annData , splitGeneNum = 1500):
        sub_genes_list = []
        all_expr_genes = geneList.tolist()
        for i in range(0, len(all_expr_genes), splitGeneNum):
            sub_genes_list.append(all_expr_genes[i: i+splitGeneNum])

        for n, sub_list in enumerate(sub_genes_list):
            Json_Data = {}
            Json_Data["cluster_plots"] = {}
            Json_Data["cluster_plots"]["gene_feature_data"] = {}
            Json_Data["cluster_plots"]["gene_feature_data"]["Gene_UMI_Values"] = []

            for i in range(len(sub_list)):
                gene_dict = {}
                gene_dict['name'] = sub_list[i]
                # coo_data = self.__exp_data[:,int(i)+ splitGeneNum * n].tocoo()
                coo_data = annData.X[:,int(i)+ splitGeneNum * n].tocoo()
                #coo_data = annData.X[:,int(i)+ splitGeneNum * n]
                gene_dict["indices"] = coo_data.row.tolist()
                #gene_dict["indices"] = coo_data.tolist()
                #gene_dict["umi_values"] = coo_data.data.tolist()
                gene_dict["umi_values"] = [int(x) for x in coo_data.data.tolist()]
                Json_Data["cluster_plots"]["gene_feature_data"]["Gene_UMI_Values"] += [gene_dict]
            yield Json_Data

    def __ReductCluster(self, adata):
        adata.obs['n_count'] = adata.X.sum(axis=1).A1
        sc.pp.normalize_total(adata,target_sum=10000)
        sc.pp.log1p(adata)
        adata.raw=adata
        #sc.pp.highly_variable_genes(adata,n_top_genes=2000)
        #sc.pp.highly_variable_genes(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
        #sc.pp.highly_variable_genes(
        #    adata,
        #    flavor="seurat",
        #    min_disp = 0.05,
        #    min_mean = 0.000001, max_mean= 6)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.regress_out(adata,['n_count'])
        sc.pp.scale(adata)
        sc.tl.pca(adata,n_comps=30, use_highly_variable=True)
        sc.pp.neighbors(adata,n_neighbors=20)
        sc.tl.leiden(adata,resolution=0.4)
        sc.tl.louvain(adata,resolution=0.4)
        sc.tl.umap(adata)
        sc.tl.tsne(adata)
        return adata

    def __fdr_lfc(self, adata):
        df_all = sc.get.rank_genes_groups_df(adata, group = None)
        df_all = df_all.set_index("names")
        df_filter = df_all.loc[(df_all['pvals'] <= 1e-5) & ((df_all['logfoldchanges'] >= 5) | (df_all['logfoldchanges'] <= -5))]
        if df_filter.shape[0] == 0:
            if df_all.shape[0] > 1000:
                df_all = df_all.sort_values(by='pvals')
                df_filter = df_all.head(1000)
            else:
                df_filter = df_all
        df_temp = df_all.loc[df_filter.index.unique()]
        total_cluster = df_temp['group'].unique()
        dfs = []
        cycle = int(len(total_cluster)/2)  
        for i in range(cycle):
            left = df_temp.loc[df_temp['group'] == str(i)]
            right = df_temp.loc[df_temp['group'] == str(2*cycle - 1 - i)]
            suffix_left = "_" + str(i)
            suffix_right = "_" + str(2*cycle-1-i)
            df_merge = pd.merge(left, right,left_index=True,right_index=True,suffixes = (suffix_left, suffix_right),how = "outer")
            dfs += [df_merge]
            if i == 0:
                df_lfc = pd.DataFrame(df_merge['logfoldchanges_0'],columns=['logfoldchanges_0'])
                df_lfc['logfoldchanges'+ suffix_right] = pd.DataFrame(df_merge['logfoldchanges' + suffix_right])
            else:
                df_lfc['logfoldchanges' + suffix_left] = pd.DataFrame(df_merge['logfoldchanges' + suffix_left])
                df_lfc['logfoldchanges' + suffix_right] = pd.DataFrame(df_merge['logfoldchanges' + suffix_right])
        df_out = pd.concat(dfs, axis=1, sort=False)
        if len(total_cluster)%2 != 0:   
            df_temp = df_temp.loc[df_temp['group'] == str(len(total_cluster)-1)]
            group = "group_" + str(len(total_cluster)-1)
            score = "scores_" + str(len(total_cluster)-1)
            logfoldchanges = "logfoldchanges_" + str(len(total_cluster)-1)
            pvals = "pvals_" + str(len(total_cluster)-1)
            pvals_adj = "pvals_adj_" + str(len(total_cluster)-1)
            df_temp.columns = Series([group, score, logfoldchanges, pvals, pvals_adj])
            df_out = pd.merge(df_out,df_temp, right_index=True, left_index=True, how = "outer")
            df_lfc[logfoldchanges] = pd.DataFrame(df_temp[logfoldchanges])
        
        # Handle NaN values - Key new code
        #print(f"├── Number of genes before clustering: {len(df_lfc)}")
        #print(f"NaN check: {df_lfc.isnull().sum().sum()} NaN values")
        #self.mobilogger._mobilogrecorder(log_message=f"├── Number of genes before clustering: {len(df_lfc)}", log_level="INFO")
        #self.mobilogger._mobilogrecorder(log_message=f"├── NaN check: {df_lfc.isnull().sum().sum()} NaN values", log_level="INFO")
        # Method 1: Drop rows with NaN (recommended)
        df_lfc_clean = df_lfc.dropna()
        #print(f"Number of genes after cleaning NaN: {len(df_lfc_clean)}")
        #self.mobilogger._mobilogrecorder(log_message=f"├── Number of genes after cleaning NaN: {len(df_lfc_clean)}", log_level="INFO")
        # If the number of genes is too small after cleaning, use filling method
        if len(df_lfc_clean) < 10:
            #print("Warning: Too few genes after cleaning NaN, using filling method")
            self.mobilogger._mobilogrecorder(log_message=f"├── Too few genes after cleaning NaN, using filling method", log_level="WARNING")
            # Method 2: Fill NaN values with 0
            df_lfc_clean = df_lfc.fillna(0)
        
        # Check if there are still issues
        if len(df_lfc_clean) < 3:
            #print("Error: Insufficient number of genes for clustering")
            self.mobilogger._mobilogrecorder(log_message=f"├── Insufficient number of genes for clustering", log_level="WARNING")
            return df_out, df_lfc.index.tolist()
            
        try:
            matrix = linkage(df_lfc_clean, method='ward', metric='euclidean')
            labels = fcluster(matrix, 3, criterion='maxclust')
            df_lfc_clean['labels'] = labels
            df_lfc_clean.sort_values(by=['labels', 'logfoldchanges_0'], inplace=True, ascending=[True, False])
            return df_out, df_lfc_clean.index.tolist()
        except Exception as e:
            #print(f"Hierarchical clustering failed: {e}")
            #print("Returning unclustered gene list")
            self.mobilogger._mobilogrecorder(log_message=f"├── Hierarchical clustering failed.", log_level="WARNING")
            self.mobilogger._mobilogrecorder(log_message=f"├── Returning unclustered gene list.", log_level="WARNING")
            # If clustering fails, at least return sorted by the first column
            df_lfc_sorted = df_lfc.dropna().sort_values(by='logfoldchanges_0', ascending=False)
            return df_out, df_lfc_sorted.index.tolist()

    def __leiden_cluster(self, adata, bdata, Json_Data):
        total_cell = len(adata.obs['leiden'])
        count_cluster = len(set(adata.obs['leiden']))
        cell_cluster = pd.DataFrame(adata.obs['leiden'].values)
        leiden_data = {}
        leiden_data["data"] = []
        leiden_cluster = {}
        leiden_cluster["data"] = {}
        leiden_cluster["data"]["fc_data"] = []
        leiden_cluster["data"]["clusters"] = []
        
        # 检查每个聚类的细胞数量
        cluster_sizes = cell_cluster[0].value_counts()
        valid_clusters = cluster_sizes[cluster_sizes >= 2].index.tolist()  # 至少2个细胞才能做差异分析
        
        for j in range(count_cluster):
            cluster_name = "cluster" + str(int(j)+1).rjust(2,'0')
            leiden_cluster["data"]["clusters"] += [cluster_name]
            bc_indice = cell_cluster.loc[cell_cluster[0] == str(j)].index
            cluster_rate = str(round(int(len(bc_indice))/int(total_cell)*100,2)) + "%"
            text = cluster_name + ": " + cluster_rate
            cluster_data = {}
            cluster_data["indices"] = bc_indice.tolist()
            cluster_data["name"] = cluster_name
            cluster_data["text"] = text
            leiden_data["data"] += [cluster_data]
        leiden_data["key"] = "leidenclust"
        leiden_data["name"] = "Leiden-Graph-Based"
        Json_Data["cluster_plots"]["cell_clusterings"] += [leiden_data]

        # 只有当有足够的有效聚类时才进行差异基因分析
        if len(valid_clusters) >= 2:
            try:
                sc.tl.rank_genes_groups(adata, "leiden", use_raw=True, method = 't-test_overestim_var')
                df_out, gene_list = self.__fdr_lfc(adata)
                res_gene_list = list(filter(lambda x: x in bdata.var_names, gene_list))
                for gene in res_gene_list:
                    gene_dict = {}
                    gene_dict['genename'] = gene
                    gene_dict['id'] = bdata.var.loc[gene]['gene_ids']
                    gene_dict['b_index'] = '0'
                    gene_set = df_out.loc[gene]
                    for cluster in range(count_cluster):
                        fc_temp = str(cluster) + "_0"
                        try:
                            gene_dict[fc_temp] = float(gene_set['logfoldchanges_'+str(cluster)])
                        except (KeyError, ValueError):
                            gene_dict[fc_temp] = 0.0  # 默认值
                        p_temp = str(cluster) + "_1"
                        try:
                            gene_dict[p_temp] = float(gene_set['pvals_adj_'+str(cluster)])
                        except (KeyError, ValueError):
                            gene_dict[p_temp] = 1.0  # 默认p值
                    leiden_cluster["data"]["fc_data"] += [gene_dict]
            except (ValueError, Exception) as e:
                #print(f"Leiden聚类差异基因分析失败: {e}")
                #print("跳过差异基因分析，仅保留聚类信息")
                self.mobilogger._mobilogrecorder(log_message=f"├── DEG analysis of Leiden failed. Skiping... {e}", log_level="WARNING")
                # 不添加差异基因数据，仅保留聚类信息
        else:
            #print(f"Leiden聚类中有效聚类数量不足 ({len(valid_clusters)} < 2)，跳过差异基因分析")
            self.mobilogger._mobilogrecorder(log_message="├── No enough cluster of Leiden for DEG analysis. Skiping...", log_level="WARNING")

        leiden_cluster["key"] = "leidenclust"
        leiden_cluster["name"] = "Leiden-Graph-Based"
        Json_Data["differential_tables_FDR"]["clusterings"] += [leiden_cluster]
        return Json_Data

    def __louvain_cluster(self, adata, bdata, Json_Data):
        total_cell = len(adata.obs['louvain'])
        count_cluster = len(set(adata.obs['louvain']))
        cell_cluster = pd.DataFrame(adata.obs['louvain'].values)
        louvain_data = {}
        louvain_data["data"] = []
        louvain_cluster = {}
        louvain_cluster["data"] = {}
        louvain_cluster["data"]["fc_data"] = []
        louvain_cluster["data"]["clusters"] = []
        
        # 检查每个聚类的细胞数量
        cluster_sizes = cell_cluster[0].value_counts()
        valid_clusters = cluster_sizes[cluster_sizes >= 2].index.tolist()
        
        for j in range(count_cluster):
            cluster_name = "cluster"+str(int(j)+1).rjust(2,'0')
            louvain_cluster["data"]["clusters"] += [cluster_name]
            bc_indice = cell_cluster.loc[cell_cluster[0] == str(j)].index
            cluster_rate = str(round(int(len(bc_indice))/int(total_cell)*100,2)) + "%"
            text = cluster_name + ": " + cluster_rate
            cluster_data = {}
            cluster_data["indices"] = bc_indice.tolist()
            cluster_data["name"] = cluster_name
            cluster_data["text"] = text
            louvain_data["data"] += [cluster_data]
        louvain_data["key"] = "louvainclust"
        louvain_data["name"] = "Louvain-Graph-Based"
        Json_Data["cluster_plots"]["cell_clusterings"] += [louvain_data]

        # 只有当有足够的有效聚类时才进行差异基因分析
        if len(valid_clusters) >= 2:
            try:
                sc.tl.rank_genes_groups(adata, "louvain", use_raw=True, method = 't-test_overestim_var')
                df_out, gene_list = self.__fdr_lfc(adata)
                res_gene_list = list(filter(lambda x: x in bdata.var_names, gene_list))
                
                for gene in res_gene_list:
                    gene_dict = {}
                    gene_dict['genename'] = gene
                    gene_dict['id'] = bdata.var.loc[gene]['gene_ids']
                    gene_dict['b_index'] = '0'
                    gene_set = df_out.loc[gene]
                    for cluster in range(count_cluster):
                        fc_temp = str(cluster) + "_0"
                        try:
                            gene_dict[fc_temp] = float(gene_set['logfoldchanges_'+str(cluster)])
                        except (KeyError, ValueError):
                            gene_dict[fc_temp] = 0.0
                        p_temp = str(cluster) + "_1"
                        try:
                            gene_dict[p_temp] = float(gene_set['pvals_adj_'+str(cluster)])
                        except (KeyError, ValueError):
                            gene_dict[p_temp] = 1.0
                    louvain_cluster["data"]["fc_data"] += [gene_dict]
                    
            except (ValueError, Exception) as e:
                #print(f"Louvain聚类差异基因分析失败: {e}")
                #print("跳过差异基因分析，仅保留聚类信息")
                self.mobilogger._mobilogrecorder(log_message=f"├── DEG analysis of Louvain failed. Skiping...{e}", log_level="WARNING")
        else:
            #print(f"Louvain聚类中有效聚类数量不足 ({len(valid_clusters)} < 2)，跳过差异基因分析")
            self.mobilogger._mobilogrecorder(log_message="├── No enough cluster of Louvain for DEG analysis. Skiping...", log_level="WARNING")
            
        louvain_cluster["key"] = "louvainclust"
        louvain_cluster["name"] = "Louvain-Graph-Based"
        Json_Data["differential_tables_FDR"]["clusterings"] += [louvain_cluster]
        return Json_Data

    def __kmeans_cluster(self, adata, bdata, min_cls_num, max_cls_num):
        total_cell = len(adata.obs['leiden'])
        for i in range(min_cls_num, max_cls_num):        ## minimal 2, maximal 111
            kmeans_tsne = KMeans(n_clusters=int(i), \
                        init='k-means++',\
                        n_init=10, \
                        max_iter=300, \
                        tol=0.0001, \
                        verbose=0, \
                        random_state=None, \
                        copy_x=True, \
                        algorithm='auto').fit_predict(adata.obsm["X_tsne"])
            count_cluster = int(i)
            group_name = 'kmeans_' + str(i) + "_clusters"
            adata.obs[group_name] = kmeans_tsne.astype('str')
            cell_cluster = pd.DataFrame(kmeans_tsne)
            
            # 检查每个聚类的细胞数量
            cluster_sizes = pd.Series(kmeans_tsne).value_counts()
            valid_clusters = cluster_sizes[cluster_sizes >= 2].index.tolist()
            
            group_data = {}
            group_data["data"] = []
            group_cluster = {}
            group_cluster["data"] = {}
            group_cluster["data"]["fc_data"] = []
            group_cluster["data"]["clusters"] = []
            for j in range(int(i)):
                cluster_name = "cluster"+str(int(j)+1).rjust(2,'0')
                group_cluster["data"]["clusters"] += [cluster_name]
                bc_indice = cell_cluster.loc[cell_cluster[0] == int(j)].index
                cluster_rate = str(round(int(len(bc_indice))/int(total_cell)*100,2)) + "%"
                text = cluster_name + ": " + cluster_rate
                cluster_data = {}
                cluster_data["indices"] = bc_indice.tolist()
                cluster_data["name"] = cluster_name
                cluster_data["text"] = text
                group_data["data"] += [cluster_data]
            group_data["key"] = group_name
            group_data["name"] = "K-Means {}".format(str(i))
            Json_Data = self.__blankJson()
            json_file = str(group_name + ".json")
            Json_Data["cluster_plots"]["cell_clusterings"] = [group_data]
            
            # 只有当有足够的有效聚类时才进行差异基因分析
            if len(valid_clusters) >= 2:
                try:
                    # 优化scanpy计算：预先过滤单细胞聚类，创建临时的聚类标签
                    cluster_sizes = pd.Series(kmeans_tsne).value_counts()
                    single_cell_clusters = cluster_sizes[cluster_sizes == 1].index.tolist()
                    
                    if single_cell_clusters:
                        # 创建过滤后的聚类标签，将单细胞聚类标记为NaN（scanpy会忽略）
                        filtered_clusters = adata.obs[group_name].copy()
                        for single_cluster in single_cell_clusters:
                            filtered_clusters[filtered_clusters == str(single_cluster)] = np.nan
                        
                        # 临时添加过滤后的聚类标签
                        temp_group_name = group_name + "_filtered"
                        adata.obs[temp_group_name] = filtered_clusters
                        
                        #print(f"├── K-means(K={i}): 跳过单细胞聚类 {single_cell_clusters}，使用过滤后的聚类进行差异分析")
                        self.mobilogger._mobilogrecorder(log_message=f"├── K-means(K={i}): Skipping cluster {single_cell_clusters}...", log_level="WARNING")
                        sc.tl.rank_genes_groups(adata, temp_group_name, use_raw=True, method = 't-test_overestim_var')
                        
                        # 清理临时标签
                        adata.obs.drop(columns=[temp_group_name], inplace=True)
                    else:
                        # 没有单细胞聚类，正常计算
                        single_cell_clusters = []  # 确保变量存在
                        sc.tl.rank_genes_groups(adata, group_name, use_raw=True, method = 't-test_overestim_var')
                    
                    df_out, gene_list = self.__fdr_lfc(adata)
                    res_gene_list = list(filter(lambda x: x in bdata.var_names, gene_list))
                    
                    for gene in res_gene_list:
                        gene_dict = {}
                        gene_dict['genename'] = gene
                        gene_dict['id'] = bdata.var.loc[gene]['gene_ids']
                        gene_dict['b_index'] = '0'
                        target_gene = df_out.loc[gene]
                        for cluster in range(count_cluster):
                            # 跳过单细胞聚类的数据填充
                            if cluster in single_cell_clusters:
                                fc_temp = str(cluster) + "_0"
                                gene_dict[fc_temp] = 0  # 单细胞聚类无差异基因信息
                                p_temp = str(cluster) + "_1"  
                                gene_dict[p_temp] = 1   # 单细胞聚类无法计算p值
                                continue
                                
                            fc_temp = str(cluster) + "_0"
                            try:
                                gene_dict[fc_temp] = float(target_gene['logfoldchanges_'+str(cluster)])
                            except (KeyError, ValueError):
                                gene_dict[fc_temp] = 0.0
                            p_temp = str(cluster) + "_1"
                            try:
                                gene_dict[p_temp] = float(target_gene['pvals_adj_'+str(cluster)])
                            except (KeyError, ValueError):
                                gene_dict[p_temp] = 1.0
                        group_cluster["data"]["fc_data"] += [gene_dict]
                        
                except (ValueError, Exception) as e:
                    #print(f"K-means(K={i})聚类差异基因分析失败: {e}")
                    #print("跳过差异基因分析，仅保留聚类信息")
                    self.mobilogger._mobilogrecorder(log_message=f"├── DEG analysis of K-means(K={i}) failed. Skiping... {e}", log_level="WARNING")
            else:
                #print(f"K-means(K={i})聚类中有效聚类数量不足 ({len(valid_clusters)} < 2)，跳过差异基因分析")
                self.mobilogger._mobilogrecorder(log_message="├── No enough cluster of K-means(K={i}) for DEG analysis. Skiping...", log_level="WARNING")

            group_cluster["key"] = group_name
            group_cluster["name"] = "K-Means {}".format(str(i))
            Json_Data["differential_tables_FDR"]["clusterings"] = [group_cluster]
            with open(json_file, "w") as outfile:
                json.dump(Json_Data, outfile, indent = 4, cls=MyEncoder)
        return True
    
    def __blankJson(self):
        Json_Data = {}
        Json_Data["cluster_plots"] = {}
        Json_Data["cluster_plots"]["cell_clusterings"] = []
        Json_Data["differential_tables_FDR"] = {}
        Json_Data["differential_tables_FDR"]["clusterings"] = []
        return Json_Data

    def __create_sample_clustering(self):
        """创建基于样本的聚类方法"""
        try:
            total_cell = len(self.__adata.obs_names)
            
            # 添加调试信息
            if 'sample' in self.__adata.obs.columns:
                sample_ids = self.__adata.obs['sample'].unique()
                #print(f"调试信息 - 数据中的样本: {sample_ids.tolist()}")
                #print(f"调试信息 - 传入的样本名: {self.__sample_name}")
                #print(f"调试信息 - 样本数量: {len(sample_ids)}")
            
            # 修改逻辑：优先使用传入的sample_name（单样本场景）
            if self.__sample_name and 'sample' in self.__adata.obs.columns:
                # 单样本场景但有sample列：使用传入的样本名
                sample_ids = self.__adata.obs['sample'].unique()
                if len(sample_ids) == 1:
                    # 确认是单样本，使用传入的样本名覆盖
                    all_indices = list(range(total_cell))
                    cluster_name = f"sample_{str(self.__sample_name)}"
                    text = f"{cluster_name}: 100.0%"
                    
                    sample_data = [{
                        "indices": all_indices,
                        "name": cluster_name,
                        "text": text
                    }]
                    
                    sample_clustering = {
                        "data": sample_data,
                        "key": "sampleclust", 
                        "name": "Sample-based"
                    }
                    #print(f"使用传入样本名创建单样本聚类: {cluster_name}")
                    return sample_clustering
            
            # 检查是否存在样本信息 - 注意：multi_convert.py中使用的是'sample'而不是'sample_id'
            if 'sample' in self.__adata.obs.columns:
                # 多样本情况：每个样本是一个cluster
                sample_ids = self.__adata.obs['sample'].unique()
                sample_data = []
                
                for i, sample_id in enumerate(sample_ids):
                    sample_indices = self.__adata.obs[self.__adata.obs['sample'] == sample_id].index.tolist()
                    # 将样本名转换为数字索引
                    indices = [self.__adata.obs_names.get_loc(idx) for idx in sample_indices]
                    
                    cluster_name = f"sample_{str(sample_id)}"
                    cluster_rate = str(round(len(indices)/total_cell*100, 2)) + "%"
                    text = f"{cluster_name}: {cluster_rate}"
                    
                    cluster_data = {
                        "indices": indices,
                        "name": cluster_name,
                        "text": text
                    }
                    sample_data.append(cluster_data)
                
                sample_clustering = {
                    "data": sample_data,
                    "key": "sampleclust",
                    "name": "Sample-based"
                }
                
            else:
                # 单样本情况：所有细胞属于一个sample cluster
                all_indices = list(range(total_cell))
                # 使用传入的样本名，如果没有则使用默认名
                if self.__sample_name:
                    cluster_name = f"sample_{str(self.__sample_name)}"
                else:
                    cluster_name = "sample_01"
                text = f"{cluster_name}: 100.0%"
                
                sample_data = [{
                    "indices": all_indices,
                    "name": cluster_name,
                    "text": text
                }]
                
                sample_clustering = {
                    "data": sample_data,
                    "key": "sampleclust", 
                    "name": "Sample-based"
                }
            
            return sample_clustering
            
        except Exception as e:
            #print(f"创建样本聚类时出错: {e}")
            self.mobilogger._mobilogrecorder(log_message="├── Clustering samples failed.", log_level="WARNING")
            return None

    def __create_sample_differential_clustering(self):
        """创建基于样本的差异分析聚类结构"""
        try:
            # 修改逻辑：优先使用传入的sample_name（单样本场景）
            if self.__sample_name and 'sample' in self.__adata.obs.columns:
                sample_ids = self.__adata.obs['sample'].unique()
                if len(sample_ids) == 1:
                    # 单样本场景，使用传入的样本名
                    clusters = [f"sample_{str(self.__sample_name)}"]
                    sample_differential_clustering = {
                        "data": {
                            "fc_data": [],  # 空的差异基因数据
                            "clusters": clusters
                        },
                        "key": "sampleclust",
                        "name": "Sample-based"
                    }
                    return sample_differential_clustering
            
            # 检查是否存在样本信息
            if 'sample' in self.__adata.obs.columns:
                sample_ids = self.__adata.obs['sample'].unique()
                # 构建clusters列表
                clusters = [f"sample_{str(sample_id)}" for sample_id in sample_ids]
            else:
                # 单样本情况
                if self.__sample_name:
                    clusters = [f"sample_{str(self.__sample_name)}"]
                else:
                    clusters = ["sample_01"]
            
            # 创建差异分析结构，但fc_data为空（因为样本间无法做差异基因分析，或单样本无差异可比较）
            sample_differential_clustering = {
                "data": {
                    "fc_data": [],  # 空的差异基因数据
                    "clusters": clusters
                },
                "key": "sampleclust",
                "name": "Sample-based"
            }
            
            return sample_differential_clustering
            
        except Exception as e:
            #print(f"创建样本差异分析聚类时出错: {e}")
            self.mobilogger._mobilogrecorder(log_message="├── DEG analysis of sample failed.", log_level="WARNING")
            return None
        
    def run(self):
        try:
            # 获取基础数据
            bdata = self.__get_bdata()
            
            # 检查数据基本信息
            total_cells = len(self.__adata.obs_names)
            total_genes = len(self.__adata.var_names)
            #print(f"数据基本信息: {total_cells} 个细胞, {total_genes} 个基因")
            
            # 检查是否有高变基因信息
            if 'highly_variable' not in self.__adata.var.columns:
                #print("警告: 数据中没有高变基因信息，这可能导致后续分析失败")
                self.mobilogger._mobilogrecorder(log_message="├── No highly variable gene detected.", log_level="WARNING")
            else:
                highly_variable_count = self.__adata.var['highly_variable'].sum()
                #print(f"高变基因数量: {highly_variable_count}")
                if highly_variable_count < 100:
                    #print(f"警告: 高变基因数量过少 ({highly_variable_count} < 100)，可能影响聚类质量")
                    self.mobilogger._mobilogrecorder(log_message="├── Too few highly variable gene detected. (%s < 100)" %(str(highly_variable_count)), log_level="WARNING")
            
            # 检查聚类信息
            cluster_methods = ['leiden', 'louvain']
            for method in cluster_methods:
                if method in self.__adata.obs.columns:
                    cluster_count = len(set(self.__adata.obs[method]))
                    #print(f"{method}聚类数量: {cluster_count}")
                    if cluster_count > total_cells * 0.5:
                        #print(f"警告: {method}聚类数量过多，可能存在过度聚类")
                        self.mobilogger._mobilogrecorder(log_message="├── Too many cluster(%s) detected." %(str(cluster_count)), log_level="WARNING")
                else:
                    #print(f"警告: 缺少{method}聚类信息")
                    self.mobilogger._mobilogrecorder(log_message="├── No cluster detected.", log_level="WARNING")
            
            # 生成降维数据 - 将直接合并到info.json中，不再单独生成Reduct.json
            dim_reduction_data = {
                "umap_x": self.__adata.obsm['X_umap'][:,0].tolist(),
                "umap_y": self.__adata.obsm['X_umap'][:,1].tolist(),
                "tsne_x": self.__adata.obsm['X_tsne'][:,0].tolist(),
                "tsne_y": self.__adata.obsm['X_tsne'][:,1].tolist()
            }

            # 生成聚类数据
            Json_Data = self.__blankJson()
            Json_Data = self.__leiden_cluster(self.__adata, bdata, Json_Data)
            leiden_Json = "leiden.json"
            with open(leiden_Json, "w") as outfile:
                json.dump(Json_Data, outfile, indent = 4, cls=MyEncoder)
            
            Json_Data = self.__blankJson()
            Json_Data = self.__louvain_cluster(self.__adata, bdata, Json_Data)
            louvain_Json = "louvain.json"
            with open(louvain_Json, "w") as outfile:
                json.dump(Json_Data, outfile, indent = 4, cls=MyEncoder)

            status = self.__kmeans_cluster(self.__adata, bdata, 2, 11)
            kmean_List = []
            for i in range(2, 11):
                kmean_List.append('kmeans_' + str(i) +'_clusters.json')

            # 生成基因表达数据
            try:
                gene_exp_list = self.__get_NoneZero_genes()
                #print(f"非零表达基因数量: {len(gene_exp_list)}")
                if len(gene_exp_list) < 500:
                    #print(f"警告: 非零表达基因数量较少 ({len(gene_exp_list)} < 500)")
                    self.mobilogger._mobilogrecorder(log_message="├── Too few non-zero gene detected. (%s < 500)" %(str(len(gene_exp_list))), log_level="WARNING")
            except ValueError as e:
                #print(f"获取基因表达数据时出错: {e}")
                # 创建空的基因列表，继续处理
                gene_exp_list = self.__adata.var_names[:1000] if len(self.__adata.var_names) > 1000 else self.__adata.var_names
                #print(f"使用前 {len(gene_exp_list)} 个基因作为备选")
                self.mobilogger._mobilogrecorder(log_message="├── No gene expression data found. Using top 1000...", log_level="WARNING")
                
            if not os.path.exists("exp"):
                os.mkdir("exp")
            exp_n = 0
            
            # 限制基因表达数据的分割大小，避免内存问题
            split_size = min(3000, max(1000, len(gene_exp_list) // 10))
            for expJsonData in self.__get_exp_each_dict(gene_exp_list, bdata, split_size):
                exp_n += 1
                with open(os.path.join("exp", str(exp_n) + ".json"), 'w') as outfile:
                    json.dump(expJsonData, outfile, indent = 4, cls=MyEncoder)

            if status:
                # 创建clusters.json - 合并所有聚类数据（包含完整的聚类和差异基因信息）
                all_clusterings = []
                all_differential_clusterings = []

                # 处理leiden聚类
                if os.path.exists(leiden_Json):
                    with open(leiden_Json,'r') as load_leiden_f:
                        leiden_data = json.load(load_leiden_f)
                        leiden_clustering = leiden_data["cluster_plots"]["cell_clusterings"][0]
                        all_clusterings.append(leiden_clustering)
                        all_differential_clusterings.extend(leiden_data["differential_tables_FDR"]["clusterings"])

                # 处理louvain聚类
                if os.path.exists(louvain_Json):
                    with open(louvain_Json,'r') as load_louvain_f:
                        louvain_data = json.load(load_louvain_f)
                        louvain_clustering = louvain_data["cluster_plots"]["cell_clusterings"][0]
                        all_clusterings.append(louvain_clustering)
                        all_differential_clusterings.extend(louvain_data["differential_tables_FDR"]["clusterings"])

                # 处理kmeans聚类
                for kmeans_file in kmean_List:
                    if os.path.exists(kmeans_file):
                        with open(kmeans_file,'r') as load_kmeans_f:
                            kmeans_data = json.load(load_kmeans_f)
                            kmeans_clustering = kmeans_data["cluster_plots"]["cell_clusterings"][0]
                            all_clusterings.append(kmeans_clustering)
                            all_differential_clusterings.extend(kmeans_data["differential_tables_FDR"]["clusterings"])

                # 添加基于样本的聚类方法
                sample_clustering = self.__create_sample_clustering()
                if sample_clustering:
                    all_clusterings.append(sample_clustering)

                # 添加基于样本的差异分析聚类结构（包含完整结构但fc_data为空）
                sample_differential_clustering = self.__create_sample_differential_clustering()
                if sample_differential_clustering:
                    all_differential_clusterings.append(sample_differential_clustering)

                # 创建clusters.json
                clusters_json_data = {
                    "cluster_plots": {
                        "cell_clusterings": all_clusterings
                    },
                    "differential_tables_FDR": {
                        "clusterings": all_differential_clusterings
                    }
                }
                
                clusters_Json = "clusters.json"
                with open(clusters_Json, "w") as outfile:
                    json.dump(clusters_json_data, outfile, indent = 4, cls=MyEncoder)

                # 创建新的info.json结构 - 更新library_type字段
                info_json_data = {
                    "exp_split": split_size,
                    "exp_count": exp_n, 
                    "library_type": "single-cell-RNA-seq",  # 更新字段名和值
                    "cell_barcodes": self.__adata.obs_names.tolist(),  # 添加细胞条形码
                    "cluster_plots": {
                        "dim_reduction_data": dim_reduction_data  # 合并降维数据
                        # 按需求移除cell_clusterings，不添加任何替代字段
                    }
                }
                
                with open('info.json', "w") as outfile:
                    json.dump(info_json_data, outfile, indent = 4, cls=MyEncoder)

                # 更新打包文件列表 - 移除Reduct.json，添加clusters.json
                package_JsonList = ['exp/*', 'clusters.json', 'info.json']

                all_json_file = ' '.join(package_JsonList)
                cmd = "zip -m {} {}".format(self.__outJsonDir, all_json_file) 
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = p.communicate()

                # 清理临时文件
                if os.path.exists("exp"): 
                    shutil.rmtree("exp")
                # 清理单独的聚类JSON文件
                for temp_file in [leiden_Json, louvain_Json] + kmean_List:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)

                zipfile = self.__outJsonDir + ".zip"
                target_f = self.__outJsonDir + ".qlentille"   # .mobi
                shutil.move(zipfile,target_f)
                
                #print("qlentille文件生成成功!")
                self.mobilogger._mobilogrecorder(log_message="├── Generating qlentille file is done successfully.", log_level="INFO")
            return True
            
        except Exception as e:
            #print(f"生成qlentille文件时发生错误: {e}")
            #print("错误详情:")
            import traceback
            self.mobilogger._mobilogrecorder(log_message="├── Generating qlentille file failed.", log_level="ERROR")
            traceback.print_exc()
            # 清理可能的临时文件
            for temp_dir in ["exp"]:
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
            
            return False
