from typing import Any
import pandas as pd
import numpy as np
import os,sys
import json
import scanpy as sc
import textwrap
import math
import matplotlib.pyplot as plt
import gc
import warnings
import re
from datetime import datetime
import multiprocessing
import copy
sys.path.append(os.path.dirname(__file__))
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor
import shutil
from sklearn.metrics import silhouette_score, pairwise_distances
from sklearn.neighbors import KernelDensity
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse import issparse
from scipy.signal import find_peaks, peak_prominences
import networkx as nx
from collections import defaultdict
from xopen import xopen

sc.settings.max_memory = 40

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr

def Hex_to_RGB(hex:str):
    hex = hex.replace("#", "")
    r = int(hex[0:2], 16)
    g = int(hex[2:4], 16)
    b = int(hex[4:6], 16)
    rgb = [r,g,b]
    return rgb

def RGB_to_Hex(RGB:list):
    color = "#"
    for i in RGB:
        if i < 1:
            i = 255 * i
        num = int(i)
        color += str(hex(num))[-2:].replace('x', '0').upper()
    return color

def RGBA_to_RGB(RGBA:list):
    R = RGBA[0]
    G = RGBA[1]
    B = RGBA[2]
    A = RGBA[3]
    R = A*R
    G = A*G
    B = A*B
    return [R,G,B]

def continuous_color(start_color:str, end_color:str, rate:float):
    start_cor = Hex_to_RGB(start_color)
    end_cor = Hex_to_RGB(end_color)
    return_cor = []
    for i in range(len(start_cor)):
        if start_cor[i] > end_cor[i]:
            delta_value = (start_cor[i] - end_cor[i]) * (rate)
            tmp_value = round(start_cor[i] - delta_value)
            if tmp_value > start_cor[i]:
                tmp_value = start_cor[i]
            if tmp_value < end_cor[i]:
                tmp_value = end_cor[i]
            return_cor.append(tmp_value)
        else:
            delta_value = (end_cor[i] - start_cor[i]) * (rate)
            tmp_value = round(end_cor[i] - delta_value)
            if tmp_value > end_cor[i]:
                tmp_value = end_cor[i]
            if tmp_value < start_cor[i]:
                tmp_value = start_cor[i]
            return_cor.append(tmp_value)
    return RGB_to_Hex(return_cor)

def softmax_vectorized(apply_array, species_list, Temp=1):
    """
    Vectorized softmax calculation, replacing the original per-dictionary processing
    
    Parameters:
    -----------
    apply_array : np.ndarray, shape (n_cells, n_species)
        Alignment counts for each cell across all species
    species_list : list
        List of species names corresponding to columns of apply_array
    Temp : float
        Temperature parameter
    
    Returns:
    --------
    species_assignments : np.ndarray, shape (n_cells,)
        Assigned species index for each cell, -1 indicates uncertain assignment
    """
    # Handle zero-sum rows
    row_sums = apply_array.sum(axis=1)
    non_zero_mask = row_sums > 0

    # Initialize results with -1 (unknown)
    result = np.full(len(apply_array), -1, dtype=int)

    if not non_zero_mask.any():
        return result

    # Only process non-zero rows
    valid_data = apply_array[non_zero_mask]

    # Auto-scaling (simulating original logic)
    max_vals = valid_data.max(axis=1, keepdims=True)
    log10_max = np.ceil(np.log10(np.maximum(max_vals, 1)))
    scales = np.power(10, np.maximum(log10_max - 3, 0))

    # Scale and apply temperature
    scaled_data = valid_data / scales / Temp

    # Compute softmax
    exp_data = np.exp(scaled_data - scaled_data.max(axis=1, keepdims=True))  # Numerical stability
    softmax_probs = exp_data / exp_data.sum(axis=1, keepdims=True)

    # Determine maximum probability indices and values
    max_probs = softmax_probs.max(axis=1)
    max_indices = softmax_probs.argmax(axis=1)

    # Get the original UMI counts (before scaling) for ratio check
    # Sort each row to get top 2 UMI counts
    sorted_umis = np.sort(valid_data, axis=1)[:, ::-1]  # Descending order
    top_umi = sorted_umis[:, 0]      # Highest UMI count
    second_umi = sorted_umis[:, 1]   # Second highest UMI count

    # Check if second highest UMI is 0 (avoid division by zero)
    # If second highest is 0, any positive top_umi satisfies "more than 2x"
    umi_ratio_check = np.where(
        second_umi == 0,
        top_umi > 0,  # If second is 0, check if top is positive
        top_umi > 2 * second_umi  # Otherwise check 2x ratio
    )

    # Combined confidence check: both probability >= 0.9 AND UMI > 2x
    confident_mask = (max_probs >= 0.9) & umi_ratio_check

    result[non_zero_mask] = np.where(confident_mask, max_indices, -1)

    return result

def process_secondary_analysis(filter_data, species_info, cell_stat_data, output_path, threads, mobilogger, sim_method="pearson", gene_filter=None, h5ad_out_file=None):
    mobilogger._mobilogrecorder(log_message="├── Running secondary analysis...", log_level="INFO")
    os.environ.pop("OMP_DISPLAY_ENV", None)
    species_info.pop('all', None)
    species_info.pop('unknown', None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        warnings.filterwarnings("ignore")
        try_filter_data = filter_data.copy()
        ###filter cells with gene count less than 30, more than 10000
        #sc.pp.filter_cells(try_filter_data, min_genes=30)
        #sc.pp.filter_cells(try_filter_data, max_genes=10000)
        #mobilogger._mobilogrecorder(log_message="├──├── Filtering cell by gene count is done successfully.", log_level="INFO")
        ####filter cells with UMI count less than 50
        #sc.pp.filter_cells(try_filter_data, min_counts=50)
        #mobilogger._mobilogrecorder(log_message="├──├── Filtering cell by UMI count is done successfully.", log_level="INFO")
        ####filter genes with present less in 3 cells
        #sc.pp.filter_genes(try_filter_data, min_cells=3)
        #mobilogger._mobilogrecorder(log_message="├──├── Filtering gene by min expression is done successfully.", log_level="INFO")
        #if try_filter_data.n_obs >= 100 and try_filter_data.n_vars >= 100:
        #    filter_data = try_filter_data.copy()
        #sc.pp.normalize_total(filter_data,target_sum=1e6, exclude_highly_expressed=False, max_fraction=0.1)
        #mobilogger._mobilogrecorder(log_message="├──├── Processing of normalization is done successfully.", log_level="INFO")
        #sc.pp.log1p(filter_data)
        #mobilogger._mobilogrecorder(log_message="├──├── Processing of log1p is done successfully.", log_level="INFO")
        ##try_filter_data = filter_data.copy()
        #sc.pp.regress_out(filter_data,['n_count'])
        #mobilogger._mobilogrecorder(log_message="├──├── Processing of regress is done successfully.", log_level="INFO")
        #sc.pp.scale(filter_data)
        #mobilogger._mobilogrecorder(log_message="├──├── Processing of scale is done successfully.", log_level="INFO")
        if len(species_info) == 1:
            expect_bins = 4
            reso = 0.05
        else:
            expect_bins = len(species_info) * 2
            reso = 0.5
        apply_high_variable_gene_list = []
        if len(species_info) >= 2:
            sc.pp.normalize_total(filter_data,target_sum=1e6, exclude_highly_expressed=False, max_fraction=0.1)
            sc.pp.log1p(filter_data)
            variable_gene_dict = {}
            species_barcode = {}
            if cell_stat_data is not None:
                for i in cell_stat_data.index:
                    tmp_species = cell_stat_data.loc[i, "species"]
                    if tmp_species != "unknown":
                        if not tmp_species in species_barcode.keys():
                            species_barcode[tmp_species] = [cell_stat_data.loc[i, "barcode"]]
                        else:
                            species_barcode[tmp_species].append(cell_stat_data.loc[i, "barcode"])
            else:
                species_barcode["unknown"] = filter_data.obs_names.tolist()
            for i in species_barcode.keys():
                high_expr_genes = []
                highly_variable_genes = []
                if len(filter_data.obs_names.isin(species_barcode[i])) != 0:
                    try:
                        sub_filter_data = filter_data[filter_data.obs_names.isin(species_barcode[i])]
                        sub_filter_data = sub_filter_data[:, sub_filter_data.var_names.str.startswith(i)]
                        sc.pl.highest_expr_genes(sub_filter_data, n_top=30)
                        gene_totals = np.ravel(sub_filter_data.X.sum(axis=0))
                        top_idx = gene_totals.argsort()
                        if len(top_idx) >= 20:
                            apply_index = 20
                        else:
                            apply_index = len(top_idx)
                        high_expr_genes = sub_filter_data.var_names[top_idx].tolist()
                        high_expr_genes = [x for x in high_expr_genes if x.startswith(i)]
                        sc.pp.highly_variable_genes(sub_filter_data, n_top_genes=apply_index)
                        highly_variable_genes = sub_filter_data.var[sub_filter_data.var['highly_variable']].index.tolist()
                        highly_variable_genes = [x for x in highly_variable_genes if x.startswith(i)]
                        variable_gene_dict[i] = list(set(high_expr_genes + highly_variable_genes))
                    except ZeroDivisionError:
                        high_expr_genes = []
                        highly_variable_genes = []
                    except Exception as e:
                        mobilogger._mobilogrecorder(log_message="├──├── ERROR occured when detecting highly variable gene of %s.\n %s\n Skip" %(e, i), log_level="WARNING")
                        high_expr_genes = []
                        highly_variable_genes = []
                apply_high_variable_gene_list += high_expr_genes + highly_variable_genes
            with open(os.path.join(output_path, "species_high_variable_genes.json"), "w") as json_file:
                json.dump(variable_gene_dict, json_file, indent=4)
            apply_high_variable_gene_list = list(set(apply_high_variable_gene_list))
            filter_data.var['highly_variable'] = filter_data.var_names.isin(apply_high_variable_gene_list)
            sc.pp.regress_out(filter_data,['n_count'])
            sc.pp.scale(filter_data)
        if apply_high_variable_gene_list == []:
            filter_data = try_filter_data.copy()
            sc.pp.filter_cells(filter_data, min_genes=30)
            sc.pp.filter_cells(filter_data, max_genes=10000)
            sc.pp.filter_cells(filter_data, min_counts=50)
            min_cells = min(int(filter_data.n_obs * 0.2), 10)
            sc.pp.filter_genes(filter_data, min_cells=min_cells)
            if filter_data.n_obs < 100 or filter_data.n_vars < 100:
                filter_data = try_filter_data.copy()
                mobilogger._mobilogrecorder(log_message="├──├── No filter applied.", log_level="INFO")
            sc.pp.normalize_total(filter_data,target_sum=1e6, exclude_highly_expressed=False, max_fraction=0.1)
            sc.pp.log1p(filter_data)
            sc.pp.highly_variable_genes(filter_data,n_bins=expect_bins)
            highly_variable_genes = filter_data.var[filter_data.var['highly_variable']]
            if len(highly_variable_genes) < 10:
                sc.pp.highly_variable_genes(filter_data)
            sc.pp.regress_out(filter_data,['n_count'])
            sc.pp.scale(filter_data)
        mobilogger._mobilogrecorder(log_message="├──├── Detecting of highly variable gene is done successfully.", log_level="INFO")
        highly_variable_genes = filter_data.var[filter_data.var['highly_variable']]
        with open(os.path.join(output_path, "highly_variable_genes.tsv"), "w") as f:
            for i in highly_variable_genes.index.tolist():
                f.write(i + "\n")
        try:
            #if len(highly_variable_genes) > 5:
            sc.tl.pca(filter_data, use_highly_variable=True)
        except Exception as e:
            #npcs = 2
            sc.tl.pca(filter_data)
            npcs = 3
        else:
            if filter_data.obsm['X_pca'].shape[1] > 10:
                npcs = 10
            else:
                npcs = filter_data.obsm['X_pca'].shape[1]
        mobilogger._mobilogrecorder(log_message="├──├── Processing of PCA is done successfully.", log_level="INFO")
        sc.pp.neighbors(filter_data, n_pcs = npcs)
        mobilogger._mobilogrecorder(log_message="├──├── Processing of neighbors is done successfully.", log_level="INFO")
        sc.tl.leiden(filter_data, resolution=reso)
        sc.tl.louvain(filter_data, resolution=reso)
        mobilogger._mobilogrecorder(log_message="├──├── Processing of leiden and louvain is done successfully.", log_level="INFO")
        sc.tl.umap(filter_data)
        sc.tl.tsne(filter_data)
        mobilogger._mobilogrecorder(log_message="├──└── Processing of UMAP and t-SNE is done successfully.", log_level="INFO")
        ###add species info
        if cell_stat_data is not None:
            for i in filter_data.obs.index.tolist():
                tmp_barcode = filter_data.obs.loc[i, "barcode"]
                tmp_k = cell_stat_data.loc[:,"barcode"].tolist().index(tmp_barcode)
                filter_data.obs.loc[i, "species"] = cell_stat_data.loc[tmp_k, "species"]
        else:
            for i in filter_data.obs.index.tolist():
                filter_data.obs.loc[i, "species"] = "unknown"
        filter_data.write(h5ad_out_file)
    return filter_data

class NpEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Any:
        if isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        elif isinstance(o, np.ndarray):
            return o.tolist()
        else:
            return super(NpEncoder, self).default(o)

def sanitize_keys(obj):
    if isinstance(obj, dict):
        return {int(k) if isinstance(k, np.integer) else str(k): sanitize_keys(v)
                for k, v in obj.items()}
    if isinstance(obj, list):
        return [sanitize_keys(i) for i in obj]
    return obj

def process_gtf_line(line_list:list, by_anno:list, split_anno:list, secondary_split_anno:list, with_prefix:bool, default_species:str):
    split_dict = {}
    for line in line_list:
        temp_line = line.replace("\n","").split("\t")
        anno = temp_line[8]
        anno_list = anno.strip(";").split(";")
        temp_by_anno = "unknown"
        temp_split_anno = "NA"
        temp_split_anno_secondary = "NA"
        for i in range(len(anno_list)):
            if not pd.isna(anno_list[i]):
                j = anno_list[i].split('"')
                if len(j) >= 2:
                    j[0] = j[0].strip()
                    j[1] = j[1].replace('"','').strip()
                    if j[0] in by_anno:
                        temp_by_anno = j[1]
                    elif j[0] in split_anno:
                        temp_split_anno = j[1]
                        if with_prefix:
                            tmp_species = temp_split_anno.split("_")[0].strip("_")
                        else:
                            tmp_species = default_species
                    elif j[0] in secondary_split_anno:
                        temp_split_anno_secondary = j[1]
        if temp_split_anno == "NA" and temp_split_anno_secondary != "NA":
            temp_split_anno = temp_split_anno_secondary
            if with_prefix:
                tmp_species = temp_split_anno.split("_")[0].strip("_")
            else:
                tmp_species = default_species
        if temp_split_anno != "NA":
            if not temp_split_anno in split_dict.keys():
                split_dict[temp_split_anno] = {"type":temp_by_anno, "species":tmp_species}
            elif split_dict[temp_split_anno]["type"] == "unknown" and temp_by_anno != "unknown":
                split_dict[temp_split_anno]["type"] = temp_by_anno
    return split_dict

def readin_custom_mtx(mobilogger, i_dir:str, mtx_file:str, feature_file:str, barcode_file:str, ref:dict, input_type:str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        species_list = ref['genomes']
        gene_info = []
        gene_name = []
        with xopen(i_dir + "/" + feature_file, "rt") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                s_line = line.split("\t")
                if not s_line[1].split("_")[0] in species_list and len(species_list) > 1:
                    s_line[1] = "_".join([species_list[0], s_line[1]])
                gene_info.append(s_line[0])
                gene_name.append(s_line[1])
        if input_type == "raw":
            mtx_data = pd.DataFrame(columns=["gene_name", "UMI_count"], index=gene_name)
            mtx_data.loc[:,"gene_name"] = gene_name
            mtx_data.loc[:,"UMI_count"] = 0
            n = 0
            with xopen(i_dir + "/" + mtx_file, "rt") as f:
                while True:
                    line = f.readline()
                    if not line:
                        break
                    if line.startswith("%"):
                        continue
                    n += 1
                    if n == 1:
                        continue
                    line = line.replace("\n","").split(" ")
                    tmp_name = gene_name[int(line[0])-1]
                    mtx_data.loc[tmp_name, "UMI_count"] += float(line[2])
            mtx_data = mtx_data[mtx_data["UMI_count"] != 0]
            return_gene_name = mtx_data.index.tolist()
            detected_dict = {}
            if len(species_list) > 1:
                for i in return_gene_name:
                    tmp_species = i.split("_")[0]
                    if not tmp_species in detected_dict.keys():
                        detected_dict[tmp_species] = 1
                    else:
                        detected_dict[tmp_species] += 1
            else:
                detected_dict[species_list[0]] = len(return_gene_name)
            return mtx_data, return_gene_name, "sparse"
        else:
            mobilogger._mobilogrecorder(log_message="├──├── Reading filtered mtx...", log_level="INFO")
            mtx_data = sc.read(i_dir + "/" + mtx_file)
            cell_info = []
            with xopen(i_dir + "/" + barcode_file, "rt") as f:
                while True:
                    line = f.readline()
                    if not line:
                        break
                    cell_info.append(line.replace("\n",""))
            mobilogger._mobilogrecorder(log_message="├──├── Creating an annoData object...", log_level="INFO")
            X = mtx_data.X.T
            raw_data = sc.AnnData(X, obs=cell_info, var = gene_info)
            raw_data.obs_names = cell_info
            raw_data.var_names = gene_info
            raw_data.obs.columns = ["barcode"]
            raw_data.var.columns = ["gene"]
            raw_data.var.index = gene_name
            raw_data.var["gene_ids"] = gene_info
            raw_data = raw_data.copy()
            mobilogger._mobilogrecorder(log_message="├──├── Removing zero-expression genes in mtx...", log_level="INFO")
            qc_metrics = sc.pp.calculate_qc_metrics(raw_data, qc_vars=[], percent_top=None, log1p=False, inplace=False)
            gene_qc = qc_metrics[1]
            raw_data.var = gene_qc
            non_zero_genes = raw_data.var_names[raw_data.var["n_cells_by_counts"] != 0].tolist()
            sc.pp.filter_genes(raw_data, min_cells=1)
            mobilogger._mobilogrecorder(log_message="├──└── Creating an annoData object is done successfully.", log_level="INFO")
            detected_dict = {}
            if len(species_list) > 1:
                for i in non_zero_genes:
                    tmp_species = i.split("_")[0]
                    if not tmp_species in detected_dict.keys():
                        detected_dict[tmp_species] = 1
                    else:
                        detected_dict[tmp_species] += 1
            else:
                detected_dict[species_list[0]] = len(non_zero_genes)
            return raw_data, non_zero_genes, "dataframe"


class Data_InjectTool:
    def __init__(self, sample_ID, output_path, pre_stat_file, input_stat_file, filter_stats_file, saturation_file, 
                 raw_mtx_path, filter_mtx_path, cell_stat_file, kit, reference_json_path, run_cmd, threads, 
                 run_cluster=True, last_json='NA', gene_type_file=None, multiplet_method="auto", host_remove=False, 
                 mobilogger=None, dev_mod=False, Temperature = 2):
        self.output_path = output_path
        if mobilogger == None:
            self.mobilogger = MobiLoggingSystem(o_dir=self.output_path, dev_mode=dev_mod)
        else:
            self.mobilogger = mobilogger
        self.mobicommandlogger = MobiCommandLogSystem(o_dir=self.mobilogger.working_path, dev_mode=False)
        self.mobiexecutor = CommandExecutor(log_system=self.mobicommandlogger, console_output=False)
        self.sample_ID = sample_ID
        self.reference = reference_json_path
        self.pre_stat_file = pre_stat_file
        self.input_data_file = input_stat_file
        self.filter_stats_file = filter_stats_file
        self.saturation_file = saturation_file
        self.raw_mtx_path = raw_mtx_path
        self.filter_mtx_path = filter_mtx_path
        self.cell_stat_file = cell_stat_file
        self.__vers_kit = kit
        self.run_cmd = run_cmd
        self.version = "v1.3.2"
        self.dev_mod = dev_mod
        self.run_cluster = run_cluster
        self.threads = int(threads)
        if not gene_type_file is None:
            self.gene_type_file = gene_type_file
        else:
            self.gene_type_file = None
        self.multiplet_method = multiplet_method
        self.host_remove = host_remove
        os.environ["PYTHONHASHSEED"] = "0"
        os.environ['OMP_DISPLAY_ENV'] = '0'
        os.environ['OMP_VERBOSE'] = '0'
        os.environ["OMP_NUM_THREADS"] = "1"
        np.random.seed(0)
        if last_json != "NA" and os.path.exists(last_json):
            with open(last_json, "r", encoding="utf-8") as json_file:
                self.o_json = json.load(json_file)
        else:
            self.o_json = {}
        self.down_list = list(np.arange(0, 1, 0.1))
        self.Temperature = Temperature

    def run_export(self):
        out_jsonf = self.export_json_microbe()    
        return out_jsonf
    
    def split_gtf(self, gtf_file:str, split_anno:list, secondary_split_anno:list, by_anno:list, o_dir:str, with_prefix:bool, default_species:str):
        split_dict = {}
        if self.threads <= 4:
            tmp_threads = self.thread
        else:
            tmp_threads = 4
        pool = multiprocessing.Pool(processes = tmp_threads)
        all_task = []
        line_list = []
        with open(gtf_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if "#" == line[0]:
                    continue
                line_list.append(line)
                if len(line_list) == 1000:
                    tmp_args = [line_list, by_anno, split_anno, secondary_split_anno, with_prefix, default_species]
                    all_task.append(pool.apply_async(process_gtf_line, tmp_args))
                    line_list = []
        pool.close()
        out_dict = {}
        for res in all_task:
            stat = res.get()
            split_dict = stat
            for i in split_dict.keys():
                if not i in out_dict.keys():
                    out_dict[i] = split_dict[i]
                elif out_dict[i]["type"] == "unknown" and split_dict[i]["type"] != "unknown":
                    out_dict[i]["type"] = split_dict[i]["type"]
        pool.join()
        with open(o_dir + '/gene_type.json', "w") as json_file:
            json.dump(out_dict, json_file, indent=4)
        return o_dir + '/gene_type.json'

    def get_gene_type_and_species(self, gene_name, empty_anno, run_type, feature_type_list, feature_dict, species_list):
        """
        返回 (gene_type, species)
        失败则 raise KeyError
        """
        if run_type == "df":
            if empty_anno:
                gtype = "unknown"
            else:
                gtype = feature_type_list.loc[gene_name, "gene_type"]
            if len(species_list) == 1:
                species = species_list[0]
            else:
                species = gene_name.split("_")[0]
            #species = gene_name.split("_")[0]
            return gtype, species

        elif run_type == "dict":
            if empty_anno:
                return "unknown", "unknown"
            return feature_dict[gene_name]["type"], feature_dict[gene_name]["species"]

        raise ValueError("invalid run_type")

    def get_RNA_type(self, i_data, reference, type_plot, out_file, outer_type, ref, input_type, gene_names):
        if input_type == "sparse":
            p_data = pd.DataFrame(columns=["n_count", "n_per"], index = i_data.index)
            p_data.loc[:,"n_count"] = i_data.loc[:,"UMI_count"]
            p_data["p_per"] = p_data["n_count"] / sum(p_data["n_count"])
            p_data.loc[:,"gene_name"] = gene_names
        else:
            p_data = pd.DataFrame(columns=["n_count", "n_per"], index = gene_names)
            n_count = i_data.X.sum(axis=0)
            p_data.loc[:, "n_count"] = n_count.A1  # 将稀疏矩阵的行向量转换为一维
            p_data["p_per"] = p_data["n_count"] / sum(p_data["n_count"])
            p_data.loc[:,"gene_name"] = gene_names
        species_list = ref['genomes']
        if not self.gene_type_file is None:
            gene_type_file = self.gene_type_file
        else:
            gene_type_file = "/".join(reference.split("/")[:-1]) + '/genes/gene_type.tsv'
        if not os.path.exists(gene_type_file):
            gene_type_file = "/".join(reference.split("/")[:-1]) + '/genes/gene_info.json'
        if not os.path.exists(gene_type_file):
            gtf_file = "/".join(reference.split("/")[:-1]) + '/genes/genes.gtf'
            by_anno = "gene_biotype,transcript_biotype"
            split_anno = "gene_name"
            secondary_split_anno = "gene_id"
            if len(species_list) == 1:
                tmp_with_prefix = False
                default_species = species_list[0]
            else:
                tmp_with_prefix = True
                default_species = species_list[0]
            gene_type_file = self.split_gtf(gtf_file=gtf_file, 
                                            split_anno=split_anno.split(","), 
                                            secondary_split_anno = secondary_split_anno.split(","), 
                                            by_anno=by_anno.split(","), 
                                            o_dir=self.output_path, 
                                            with_prefix=tmp_with_prefix, 
                                            default_species=default_species)
        empty_anno = True
        if os.path.dirname(gene_type_file) != self.output_path:
            #cmd = "cp %s %s" %(gene_type_file, self.output_path)
            #os.system(cmd)
            shutil.copy2(gene_type_file, self.output_path)
        feature_type_list = None
        feature_dict = None
        if gene_type_file.endswith(".tsv"):
            feature_type_list = pd.read_csv(gene_type_file, sep="\t", index_col=0)
            feature_type_list.index = range(feature_type_list.shape[0])
            for i in feature_type_list.index:
                if not str(feature_type_list.loc[i, "gene_name"]).split("_")[0] in species_list:
                    feature_type_list.loc[i, "gene_name"] = species_list[0] + "_" + feature_type_list.loc[i, "gene_name"]
            feature_type_list.index = feature_type_list.loc[:,"gene_name"].tolist()
            feature_type_list = feature_type_list[~feature_type_list.index.duplicated(keep='first')]
            run_type = "df"
            if feature_type_list.shape[0] != 0:
                empty_anno = False
        elif gene_type_file.endswith(".json"):
            with open(gene_type_file, "r", encoding="utf-8") as json_file:
                feature_dict = json.load(json_file)
            run_type = "dict"
            if len(feature_dict) != 0:
                empty_anno = False
        p_data.loc[:,"gene_type"] = "unknown"
        p_data.loc[:,"species"] = "unknown"
        p_data.index = range(p_data.shape[0])
        species_gene_dict = {}
        for i in p_data.index.tolist():
            tmp_gene_name1 = p_data.loc[i, "gene_name"]
            if len(tmp_gene_name1.split("_", 1)) > 1:
                tmp_gene_name2 = tmp_gene_name1.split("_", 1)[1].strip("_")
            else:
                tmp_gene_name2 = tmp_gene_name1.strip("_")
            for try_name in (tmp_gene_name1, tmp_gene_name2):
                try:
                    gtype, sp = self.get_gene_type_and_species(
                        try_name, empty_anno, run_type,
                        feature_type_list, feature_dict, species_list)
                except KeyError:
                    continue          # skip when failed
                else:
                    # sucess update
                    sp = sp.strip("_")
                    p_data.loc[i, "gene_name"]   = try_name   
                    p_data.loc[i, "gene_type"]   = gtype
                    p_data.loc[i, "species"]     = sp
                    species_gene_dict[sp] = species_gene_dict.get(sp, 0) + 1
                    break
            else:
                # skip when key error
                pass
        p_data.to_csv(os.path.join(self.output_path, "processed_gene_type.tsv"), sep="\t", index=False)
        type_list = pd.DataFrame(columns=["species", "gene_type", "n_count", "n_percentage"])
        for i in p_data.index.tolist():
            temp_type = p_data.loc[i, "gene_type"]
            temp_species = p_data.loc[i, "species"]
            if temp_species != "unknown":
                unique_index = temp_species + "_" + temp_type
                if not pd.isna(temp_type):
                    if not unique_index in type_list.index.tolist():
                        add_shape = unique_index
                        type_list.loc[add_shape, "species"] = temp_species
                        type_list.loc[add_shape, "gene_type"] = temp_type
                        type_list.loc[add_shape, "n_count"] = p_data.loc[i, "n_count"]
                        type_list.loc[add_shape, "n_percentage"] = 0
                    else:
                        add_shape = unique_index
                        type_list.loc[add_shape, "n_count"] += p_data.loc[i, "n_count"]
        species_sum = {}
        for i in type_list.index:
            if not type_list.loc[i, "species"] in species_sum.keys():
                species_sum[type_list.loc[i, "species"]] = type_list.loc[i, "n_count"]
            else:
                species_sum[type_list.loc[i, "species"]] += type_list.loc[i, "n_count"]
        for i in type_list.index:
            type_list.loc[i, "n_percentage"] = type_list.loc[i, "n_count"] / species_sum[type_list.loc[i, "species"]]
        out_gene_type_result_file = self.output_path + "/" + out_file
        type_list = type_list.fillna(0)
        type_list.to_csv(out_gene_type_result_file, index=False)
        ###to type plot
        type_annote = True
        std_config = {
            "displaylogo": 'false',
            "staticPlot": 'false',
            "displayModeBar": 'true',
            "showAxisDragHandles": 'true',
            "toImageButtonOptions": {},
            "scrollZoom": 'true'
        }
        std_layout = {
                        "title": {
                            "text": "Gene Type Fraction",
                            "font":  {
                                "size": 16
                            }
                        },
                        "hovermode": "closest",
                        "showlegend": 'true',
                        "yaxis": {
                            "tickformat": ",.1%", "rangemode": 'tozero', 'autorange': 'true'
                        },
                        "barmode": "stack"
                    }
        std_data = {'x':[], 'y':[], "text":[], "type": "histogram", 
                    "name": "", 
                    "hovertemplate": "%{y}", 
                    "meanline" : {"visible": 'true'},
                    "histfunc": "sum",  
                    }
        if not "config" in type_plot.keys():
            type_plot["config"] = std_config
        if not "layout" in type_plot.keys():
            type_plot["layout"] = std_layout
        if not "data" in type_plot.keys():
            type_plot["data"] = []
            hier = False
        else:
            hier = True
        added_type = {}
        for i in type_list.index.tolist():
            if type_list.loc[i, "n_percentage"] < 0.01:
                tmp_type = "other"
            else:
                tmp_type = type_list.loc[i, "gene_type"]
            found_existed = False
            if hier:
                for j in type_plot["data"]:
                    if tmp_type == j["name"] or (tmp_type == "protein_coding" and j["name"] == "mRNA"):
                        found_existed = True
                        break
                if not found_existed:
                    tmp_type = "other"
                    for j in type_plot["data"]:
                        if tmp_type == j["name"]:
                            found_existed = True
                            break
            if tmp_type == "protein_coding":
                tmp_type = "mRNA"
            if not found_existed:
                if not tmp_type in added_type.keys():
                    tmp_data = std_data.copy()
                    tmp_data["x"] = [type_list.loc[i, "species"] + "_" + outer_type]
                    tmp_data["y"] = [type_list.loc[i, "n_percentage"]]
                    tmp_data["text-d"] = [tmp_type]
                    tmp_data["name"] = tmp_type
                    added_type[tmp_type] = tmp_data.copy()
                else:
                    if type_list.loc[i, "species"] + "_" + outer_type in added_type[tmp_type]["x"]:
                        tmp_index = added_type[tmp_type]["x"].index(type_list.loc[i, "species"] + "_" + outer_type)
                        added_type[tmp_type]["y"][tmp_index] += type_list.loc[i, "n_percentage"]
                    else:
                        added_type[tmp_type]["x"].append(type_list.loc[i, "species"] + "_" + outer_type)
                        added_type[tmp_type]["y"].append(type_list.loc[i, "n_percentage"])
                        added_type[tmp_type]["text-d"].append(tmp_type)
            else:
                if type_list.loc[i, "species"] + "_" + outer_type in j["x"]:
                    tmp_index = j["x"].index(type_list.loc[i, "species"] + "_" + outer_type)
                    j["y"][tmp_index] += type_list.loc[i, "n_percentage"]
                else:
                    j["x"].append(type_list.loc[i, "species"] + "_" + outer_type)
                    j["y"].append(type_list.loc[i, "n_percentage"])
                    j["text-d"].append(j["text-d"][0])
        for i in added_type.keys():
            type_plot["data"].append(added_type[i])
        return type_plot, type_annote, species_gene_dict
    
    def get_UMAP_plot_count(self, i_data):
        umap_plot_data=dict()
        umap_plot_data["config"] = {
            "displaylogo": 'false',
            "staticPlot": 'false',
            "displayModeBar": 'true',
            "showAxisDragHandles": 'true',
            "toImageButtonOptions": {},
            "scrollZoom": 'true'
        }
        umap_plot_data["layout"] = {
            "title": {
                "text": "UMAP Projection of Microbes Colored by UMI Count",
                "font":  {
                    "size": 16
                }, 
                "y": 0.94,
                "xanchor": 'center',
                "yanchor": 'top'
            },
            "xaxis": {
                "type": "linear",
                "title": "UMAP1",
                "showline": 'false',
                "zeroline": 'true',
                "fixedrange": 'false',
                "titlefont": {
                    "size": 14
                }
            },
            "yaxis": {
                "type": "linear",
                "title": "UMAP2",
                "showline": 'false',
                "zeroline": 'true',
                "fixedrange": 'false',
                "titlefont": {
                    "size": 14
                }
            },
            "margin": {
                "t" : 50
            }
        }
        umap_plot_data["data"] = []
        tmp_data = {}
        tmp_data["type"] = "scattergl"
        tmp_data["name"] = "Cells"
        tmp_data["mode"] = "markers"
        tmp_data["x"]=i_data.obsm['X_umap'][:,0].tolist()
        tmp_data["y"]=i_data.obsm['X_umap'][:,1].tolist()
        tmp_data["text"] = list(map(math.log10, i_data.obs['n_count'].tolist()))
        tmp_data["marker"] = {}
        tmp_data["marker"]["opacity"] = 0.9
        tmp_data["marker"]["size"] = 3.8
        tmp_data["marker"]["color"] = list(map(math.log10, i_data.obs['n_count'].tolist()))
        tmp_data["marker"]["colorscale"] = "Jet"
        tmp_data["marker"]["cmax"] = math.log10(max(i_data.obs['n_count'].tolist()))
        tmp_data["marker"]["cmin"] = math.log10(min(i_data.obs['n_count'].tolist()))
        tmp_data["marker"]["colorbar"] = {"title": "log10(UMI counts)"}
        umap_plot_data["data"].append(tmp_data)
        return umap_plot_data

    def get_UMAP_plot_cluster(self, i_data):
        i_data.obs['x'] = i_data.obsm['X_umap'][:,0]
        i_data.obs['y'] =i_data.obsm['X_umap'][:,1]
        i_data.obs = i_data.obs.sort_values("leiden")
        i_data.obs.reset_index(drop=True,inplace=True)
        flag = []
        for i in range(len(i_data.obs)):
            if i >= 1:
                if i_data.obs['leiden'][i] != i_data.obs['leiden'][i-1]:
                    flag += [i]
        flag += [len(i_data.obs)]
        p_dict = {}
        p_dict["config"] = {}
        p_dict["config"]["displaylogo"] = "false"
        p_dict["config"]["staticPlot"] = "false"
        p_dict["config"]["displayModeBar"] = "true"
        p_dict["config"]["showAxisDragHandles"] = "true"
        p_dict["config"]["toImageButtonOptions"] = {}
        p_dict["config"]["scrollZoom"] = "true"
        p_dict["layout"] = {}
        p_dict["layout"]["title"] = {}
        p_dict["layout"]["title"]["text"] = "UMAP Projection of Microbes Colored by Cluster"
        p_dict["layout"]["title"]["font"] = {"size":16}
        p_dict["layout"]["title"]["y"] = 0.94
        p_dict["layout"]["title"]["xanchor"] = "center"
        p_dict["layout"]["title"]["yanchor"] = "top"
        p_dict["layout"]["hovermode"] = "closest"
        p_dict["layout"]["xaxis"] = {
                                        "type": "linear",
                                        "title": "UMAP1",
                                        "showline": 'false',
                                        "zeroline": 'true',
                                        "fixedrange": 'false',
                                        "titlefont": {
                                            "size": 14
                                        }   
                                    }
        p_dict["layout"]["yaxis"] = {
                                        "type": "linear",
                                        "title": "UMAP2",
                                        "showline": 'false',
                                        "zeroline": 'true',
                                        "fixedrange": 'false',
                                        "titlefont": {
                                            "size": 14
                                        }
                                    }
        p_dict["layout"]["margin"] = {"t" : 50}
        p_dict["data"] = []
        for i in range(len(flag)):
            tmp_add = {}
            tmp_add["type"] = "scatter"
            tmp_add["mode"] = "markers"
            tmp_add["marker"] = {}
            tmp_add["marker"]["opacity"] = 0.9
            tmp_add["marker"]["size"] = 3.8
            if i == 0:
                tmp_add['name'] = "Cluster" + str(i_data.obs['leiden'][flag[i]-1])
                tmp_add["hovertemplate"] = "(%%{x:.3f}, %%{y:.3f}) <br> %s: %%{text}<extra></extra>" %(tmp_add["name"])
                tmp_add['x'] = i_data.obs['x'][0:flag[i]].tolist()
                tmp_add['y'] = i_data.obs['y'][0:flag[i]].tolist()
                tmp_add["text"] = i_data.obs['n_count'][0:flag[i]].tolist()
            else:
                tmp_add['name'] ="Cluster" + str(i_data.obs['leiden'][flag[i]-1])
                tmp_add["hovertemplate"] = "(%%{x:.3f}, %%{y:.3f}) <br> %s: %%{text}<extra></extra>" %(tmp_add["name"])
                tmp_add['x'] = i_data.obs['x'][flag[i-1]:flag[i]].tolist()
                tmp_add['y'] = i_data.obs['y'][flag[i-1]:flag[i]].tolist()
                tmp_add['text'] = i_data.obs['n_count'][flag[i-1]:flag[i]].tolist()
            p_dict["data"].append(tmp_add)
        return p_dict
        
    def prepare_raw_mtx(self, i_dir:str, mtx_file:str, feature_file:str, barcode_file:str, ref:dict, barcode_list:list):
        species_list = ref['genomes']
        gene_info = []
        gene_name = []
        with open(i_dir + "/" + feature_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                s_line = line.split("\t")
                if not s_line[1].split("_")[0] in species_list and len(species_list) > 1:
                    s_line[1] = "_".join([species_list[0], s_line[1]])
                gene_info.append(s_line[0])
                gene_name.append(s_line[1])
        mtx_data = sc.read(i_dir + "/" + mtx_file)
        X = mtx_data.X.T
        cell_info = []
        with open(i_dir + "/" + barcode_file, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                cell_info.append(line.replace("\n",""))
        #cell_info_df = pd.DataFrame(cell_info)
        #gene_name_df = pd.DataFrame(gene_name)
        raw_data = sc.AnnData(X, obs=cell_info, var = gene_name)
        #raw_data.obs = cell_info_df
        #raw_data.var = gene_name_df
        raw_data.obs_names = cell_info
        raw_data.var_names = gene_name
        raw_data.var_names_make_unique()
        raw_data.obs['n_count'] = raw_data.X.sum(axis=1).A1
        exp_raw = raw_data[barcode_list].copy()
        return exp_raw 

    def estimate_saturation(self, gene_count_array:np.array, down_list:list):
        return_dict = {}
        for i in down_list:
            tmp_array = 1-np.power(1-i, gene_count_array)
            return_dict[i] = sum(tmp_array)
        return return_dict

    def export_json_microbe(self):
        #self.mobilogger._mobilogrecorder(log_message="Exporting report json...", log_level="INFO")
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
        species_info = {}
        h5ad_out_file = os.path.join(self.output_path, self.sample_ID + ".h5ad")
        if not "input" in self.o_json.keys():
            input_data = pd.read_csv(self.input_data_file, sep="\t")
            pre_data = pd.read_csv(self.pre_stat_file, sep="\t")
            filter_data = pd.read_csv(self.filter_stats_file, sep="\t")
            input_info = {}
            input_info["raw_reads"] = filter_data.loc[0, "processed_reads"]
            if self.dev_mod:
                input_info["valid_reads"] = input_data.loc[0, "Raw_reads"]
            else:
                input_info["valid_reads"] = filter_data.loc[0, "passed_reads"]
            input_info["Fraction_of_reads"] = filter_data.loc[0, "passed_reads"] / filter_data.loc[0, "processed_reads"]
            input_info["star_input_reads"] = input_data.loc[0, "Host_unremoved_reads"]
            if self.host_remove:
                input_info["adaptor_free_reads"] = input_data.loc[0, "Fastp_reads"]
                input_info["fraction_of_host_reads"] = 1 - input_data.loc[0, "Host_unremoved_reads"] / input_data.loc[0, "Fastp_reads"]
            input_info["saturation_rate"] = 0
            input_info["Q30_in_all"] = pre_data.loc[0, "Q30"]
            if pd.isna(input_data.loc[0, "Q30"]):
                input_info["Q30_in_valid"] = pre_data.loc[0, "Q30"]
            else:
                input_info["Q30_in_valid"] = input_data.loc[0, "Q30"]
            self.o_json["input"] = input_info
        else:
            self.mobilogger._mobilogrecorder(log_message="├── Using info from last json...", log_level="INFO")
            input_info = self.o_json["input"]
        ###readin ref json
        if self.reference != "NA":
            with open(self.reference,'r') as load_f:
                ref_dict = json.load(load_f)
                self.o_json["ref"] = ref_dict
        else:
            ref_dict = self.o_json["ref"]
        self.mobilogger._mobilogrecorder(log_message="├── Processing gene type from raw mtx...", log_level="INFO")
        for tmp_f in os.listdir(self.raw_mtx_path):
            if tmp_f.endswith(".gz"):
                cmd = "gunzip %s" %(os.path.join(self.raw_mtx_path, tmp_f))
                #os.system(cmd)
                exit_code = self.mobiexecutor.execute(
                    command=cmd,
                    context={
                        "export json": "gunzip"
                    }, console_output=False
                )
                if exit_code != 0:
                    self.mobilogger._mobilogrecorder(log_message="Gunzip failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                        log_level="ERROR")
                    sys.exit()
        raw_data, gene_names, tmp_output_type = readin_custom_mtx(mobilogger=self.mobilogger, i_dir=self.raw_mtx_path, mtx_file="matrix.mtx", feature_file="features.tsv", barcode_file="barcodes.tsv", ref=ref_dict, input_type="raw")
        ###rRNA
        type_plot = {}
        type_plot, type_annote, _ = self.get_RNA_type(i_data=raw_data, reference=self.reference, type_plot=type_plot, out_file="raw_gene_type_stats.csv", outer_type="MI", ref=ref_dict, input_type=tmp_output_type, gene_names=gene_names)
        del raw_data
        gc.collect()
        self.mobilogger._mobilogrecorder(log_message="├── Processing gene type from raw mtx is done successfully.", log_level="INFO")
        ###readin filtered mtx
        self.mobilogger._mobilogrecorder(log_message="├── Processing gene type from filtered mtx...", log_level="INFO")
        for tmp_f in os.listdir(self.filter_mtx_path):
            if tmp_f.endswith(".gz"):
                cmd = "gunzip %s" %(os.path.join(self.filter_mtx_path, tmp_f))
                #os.system(cmd)
                exit_code = self.mobiexecutor.execute(
                    command=cmd,
                    context={
                        "export json": "gunzip"
                    }, console_output=False
                )
                if exit_code != 0:
                    self.mobilogger._mobilogrecorder(log_message="Gunzip failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                        log_level="ERROR")
                    sys.exit()
        filter_data, gene_names, tmp_output_type = readin_custom_mtx(mobilogger=self.mobilogger, i_dir=self.filter_mtx_path, mtx_file="matrix.mtx", feature_file="features.tsv", barcode_file="barcodes.tsv", ref=ref_dict, input_type="filtered")
        filter_data.var_names_make_unique()
        if type_annote:
            type_plot, type_annote, species_genes = self.get_RNA_type(i_data=filter_data, reference=self.reference, type_plot=type_plot, out_file="filtered_gene_type_stats.csv", outer_type="ME", ref=ref_dict, input_type=tmp_output_type, gene_names=gene_names)
        else:
            type_plot["type"] = "NA"
        self.mobilogger._mobilogrecorder(log_message="├── Processing gene type from filtered mtx is done successfully.", log_level="INFO")
        ###species info
        self.mobilogger._mobilogrecorder(log_message="├── Processing species info...", log_level="INFO")
        species_info = {}
        species_info["all"] = {}
        species_info["all"]["cell_number"] = 0
        species_info["all"]["found"] = True
        species_info["all"]["gene_detected"] = 0
        for i in ref_dict["genomes"]:
            if not i in species_info.keys() and i in species_genes.keys():
                tmp_species = i
                species_info[tmp_species] = {}
                species_info[tmp_species]["UMI_count"] = []
                species_info[tmp_species]["gene_count"] = []
                species_info[tmp_species]["reads_count"] = []
                species_info[tmp_species]["cell_number"] = 0
                species_info[tmp_species]["map_rate"] = 0
                species_info[tmp_species]["confidently_map_rate"] = 0
                species_info[tmp_species]["found"] = False
                species_info[tmp_species]["gene_detected"] = species_genes[i]
                species_info[tmp_species]["down_sample_gene_count"] = {}
                for j in self.down_list:
                    species_info[tmp_species]["down_sample_gene_count"][j] = []
                species_info["all"]["gene_detected"] += species_genes[i]
        self.mobilogger._mobilogrecorder(log_message="├── Processing species info is done successfully.", log_level="INFO")
        filter_df = filter_data.to_df().fillna(0)
        filter_df.index = filter_data.obs.loc[:,"barcode"].tolist()
        filter_df.columns = filter_data.var.index.tolist()
        cell_stat_data = pd.read_csv(self.cell_stat_file, sep="\t")
        if "gene_UMI_count" in cell_stat_data.columns:
            UMI_col = "gene_UMI_count"
        else:
            UMI_col = "UMI_count"
        cell_stat_data = cell_stat_data.sort_values(by=UMI_col, ascending=False)
        cell_stat_data.index = range(cell_stat_data.shape[0])
        if (len(ref_dict["genomes"]) >= 10 and self.multiplet_method == "auto") or \
            self.multiplet_method == "scaled_softmax":
            self.mobilogger._mobilogrecorder(log_message="├── Re-assigning barcode-species by softmax...", log_level="INFO")
            # 向量化
            apply_array = cell_stat_data[ref_dict["genomes"]].values.astype(float)
            species_indices = softmax_vectorized(apply_array, ref_dict["genomes"], self.Temperature)
            # 直接赋值
            cell_stat_data["species"] = [
                ref_dict["genomes"][idx] if idx >= 0 else "unknown" 
                for idx in species_indices
            ]
            self.mobilogger._mobilogrecorder(log_message="├── Re-assigning barcode-species by softmax is done successfully.", log_level="INFO")
        cell_stat_data.to_csv(os.path.join(self.output_path, "barcode_info.tsv.gz"), sep="\t", index=False, compression='gzip')
        passed_barcode = filter_df.index.tolist()
        barcode_rank = 1
        ###rank plot
        cell_color = "#004799"
        end_color = "#cbe1f0"
        background_color = "#dddddd"
        log_step = 20
        rank_plot = {}
        rank_plot["config"] = {
                                "displaylogo": 'false',
                                "staticPlot": 'false',
                                "displayModeBar": 'true',
                                "showAxisDragHandles": 'true',
                                "toImageButtonOptions": {},
                                "scrollZoom": 'true'
                            }
        rank_plot["layout"] = {
                            "legend": {
                                "x":0.49,
                                "y":1.02,
                                "orientation": "h",
                                "legend" : {"traceorder": "grouped+reversed"} ,
                                "xanchor": "center",
                                "yanchor": "bottom"
                            },
                            "title": {
                                "text": "Barcode Rank Plot",
                                "font":  {
                                    "size": 16
                                }
                            },
                            "xaxis": {
                                "title": "Barcodes",
                                "type": "log",
                                "showline": 'true',
                                "zeroline": 'false',
                                "fixedrange": 'false',
                                "titlefont": {
                                    "size": 14
                                }
                            },
                            "yaxis": {
                                "title": "UMI Counts",
                                "type": "log",
                                "showline": 'true',
                                "zeroline": 'false',
                                "fixedrange": 'false',
                                "nticks": 6,
                                "tickformat": ".0f",
                                "titlefont": {
                                    "size": 14
                                }
                            },
                            "hovermode": "closest",
                            "margin": {
                                "t" : 80
                            }
                        }
        rank_plot["data"] = []
        primer_dot = {}
        primer_dot['x'] = []
        primer_dot['y'] = []
        primer_dot['type'] = "scattergl"
        primer_dot["mode"] = "lines"
        dot_template = primer_dot.copy()
        background_dot = dot_template.copy()
        current_dot = dot_template.copy()
        current_dot["name"] = "Microbes"
        current_dot["line"] = {"color":cell_color,"width":3}
        current_dot["markers"] = {"color":cell_color,"size":3}
        current_dot["showlegend"] = 'true'
        current_dot["xaxis"] = "x"
        current_dot["yaxis"] = "y"
        current_dot["hovertemplate"] = "Microbes<br>Rank: %{x}<br>UMI: %{y}<extra></extra>"
        background_dot["line"] = {"color":background_color,"width":3}
        background_dot["showlegend"] = 'true'
        background_dot["hoverinfo"] = "text"
        background_dot["text"] = ["Background"]
        background_dot["name"] = "Background"
        found_bk = False
        end_cell = False
        current_start = -1
        current_cell = -1
        current_bk = -1
        self.mobilogger._mobilogrecorder(log_message="├── Processing barcode-rank plot...", log_level="INFO")
        found_bc = 0
        drop_index = []
        for i in cell_stat_data.index:
            tmp_barcode = cell_stat_data.loc[i, "barcode"]
            if not tmp_barcode in passed_barcode:
                if not found_bk:
                    found_bk = True
                    rank_plot["data"].append(current_dot)
                    if len(current_dot['x']) > 0:
                        last_point = True
                        last_point_x = current_dot['x'][-1]
                        last_point_y = current_dot['y'][-1]
                    else:
                        last_point = False
                    current_dot = dot_template.copy()
                    current_start = cell_stat_data.loc[i, UMI_col]
                    current_dot["type"] = "scattergl"
                    current_dot["mode"] = "lines"
                    current_dot["showlegend"] = 'false'
                    current_dot["name"] = "Mix"
                    current_dot["hoverinfo"] = "text"
                    current_dot['x'] = []
                    current_dot['y'] = []
                    if last_point:
                        current_dot['x'].append(last_point_x)
                        current_dot['y'].append(last_point_y)
                    current_cell = 0
                    current_bk = 0
                else:
                    if cell_stat_data.loc[i, UMI_col] < current_start - log_step and not end_cell :
                        current_rate = current_cell / (current_bk + current_cell)
                        if current_rate > 0:
                            current_color = continuous_color(start_color=cell_color, end_color=end_color, rate=current_rate)
                            current_dot["line"] = {"color":current_color,"width":3}
                            current_dot["marker"] = {"color":current_color,"size":3}
                            show_rate = "%.2f%%" %(current_rate * 100)
                            current_dot["hovertemplate"] = "%s Microbes(%s/%s)<br>Rank:" %(show_rate, str(current_cell), str(current_cell + current_bk)) + "%{x}<br>UMI: %{y}<extra></extra>" 
                            if len(current_dot['x']) > 0:
                                last_point = True
                                last_point_x = current_dot['x'][-1]
                                last_point_y = current_dot['y'][-1]
                            else:
                                last_point = False
                            rank_plot["data"].append(current_dot)
                            current_dot = dot_template.copy()
                            current_start = cell_stat_data.loc[i, UMI_col]
                            current_dot["type"] = "scattergl"
                            current_dot["mode"] = "lines"
                            current_dot["showlegend"] = 'false'
                            current_dot["name"] = "Mix"
                            current_dot["hoverinfo"] = "text"
                            current_dot['x'] = []
                            current_dot['y'] = []
                            if last_point:
                                current_dot['x'].append(last_point_x)
                                current_dot['y'].append(last_point_y)
                            current_cell = 0
                            current_bk = 0
                        else:
                            if found_bc == len(passed_barcode):
                                current_color = background_color
                                current_dot["line"] = {"color":current_color,"width":3}
                                current_dot["showlegend"] = 'true'
                                current_dot["hoverinfo"] = "text"
                                current_dot["hovertemplate"] = "Background<br>Rank: %{x}<br>UMI: %{y}<extra></extra>"
                                current_dot["name"] = "Background"
                                end_cell = True
                current_dot['x'].append(barcode_rank)
                current_dot['y'].append(cell_stat_data.loc[i, UMI_col])
                current_bk += 1
                drop_index.append(i)
            else:
                found_bc += 1
                current_dot['x'].append(barcode_rank)
                current_dot['y'].append(cell_stat_data.loc[i, UMI_col])
                current_cell += 1
                ###species_info
                species_info["all"]["cell_number"] += 1
                tmp_species = cell_stat_data.loc[i, "species"]
                filter_df.loc[tmp_barcode, "species"] = tmp_species
                filter_df.loc[tmp_barcode, "reads_count"] = cell_stat_data.loc[i, "read_count"]
                filter_df.loc[tmp_barcode, "barcode_rank"] = barcode_rank
                if not tmp_species in species_info.keys():
                    species_info[tmp_species] = {}
                    species_info[tmp_species]["UMI_count"] = []
                    species_info[tmp_species]["gene_count"] = []
                    species_info[tmp_species]["reads_count"] = []
                    species_info[tmp_species]["cell_number"] = 0
                    species_info[tmp_species]["map_rate"] = 0
                    species_info[tmp_species]["confidently_map_rate"] = 0
                    species_info[tmp_species]["down_sample_gene_count"] = {}
                    species_info[tmp_species]["gene_detected"] = 0
                    for j in self.down_list:
                        species_info[tmp_species]["down_sample_gene_count"][j] = []
                species_info[tmp_species]["found"] = True
            barcode_rank += 1
        cell_stat_data.drop(drop_index, inplace=True)
        sum_valid_reads = cell_stat_data["read_count"].sum()
        input_info["Fraction_of_reads_in_cell"] = sum_valid_reads / input_info["raw_reads"]
        ###update color
        nn = 0
        if len(species_info) >= 10:
            rm_unfound = True
        else:
            rm_unfound = False
        drop_index = []
        for i in species_info.keys():
            if rm_unfound and not species_info[i]["found"]:
                drop_index.append(i)
            else:
                species_info[i]["line_color"] = nn
                species_info[i]["fill_color"] = nn
                nn += 1
        if len(drop_index) > 0 and len(drop_index) < len(species_info) - 2:
            for i in drop_index:
                if i in species_info.keys():
                    del species_info[i]
        else:
            for i in drop_index:
                species_info[i]["line_color"] = nn
                species_info[i]["fill_color"] = nn
                nn += 1
        tmp_array = range(len(species_info) + 1)
        tmp_array = [i/(len(species_info)+1) for i in tmp_array]
        apply_colors = plt.cm.rainbow(tmp_array)
        for i in species_info.keys():
            tmp_n = species_info[i]["line_color"]
            tmp_color = apply_colors[tmp_n]
            if len(tmp_color) == 3:
                tmp_color = RGB_to_Hex(list(tmp_color))
            else:
                tmp_color = RGB_to_Hex(RGBA_to_RGB(list(tmp_color)))
            species_info[i]["line_color"] = tmp_color
            fill_color = apply_colors[tmp_n]
            if len(fill_color) == 3:
                fill_color = list(fill_color)
                fill_color.append(0.5)
                fill_color = tuple(fill_color)
            else:
                fill_color[3] = 0.5
            fill_color = RGB_to_Hex(list(fill_color))
            species_info[i]["fill_color"] = fill_color
        current_rate = current_cell / (current_bk + current_cell)
        if current_rate > 0:
            current_color = continuous_color(start_color=cell_color, end_color=end_color, rate=current_rate)
            current_dot["line"] = {"color":current_color,"width":3}
            current_dot["marker"] = {"color":current_color,"size":3}
            show_rate = "%.2f%%" %(current_rate * 100)
            current_dot["text"] = "%s Cells<br>(%s/%s)" %(show_rate, str(current_cell), str(current_cell + current_bk))
        else:
            current_color = background_color
            current_dot["line"] = {"color":current_color,"width":3}
            current_dot["showlegend"] = 'true'
            current_dot["hoverinfo"] = "text"
            current_dot["text"] = "Background"
            current_dot["name"] = "Background"
        rank_plot["data"].append(current_dot)
        for i in filter_df.index.tolist():
            if pd.isna(filter_df.loc[i, "species"]):
                filter_df.loc[i, "species"] = "unknown"
            tmp_species = filter_df.loc[i, "species"]
            if tmp_species in species_info.keys():
                tmp_list = filter_df.loc[i,:].tolist()[:-5]
                tmp_UMI_sum = sum(tmp_list)
                tmp_gene_count = np.count_nonzero(np.array(tmp_list))
                species_info[tmp_species]["UMI_count"].append(tmp_UMI_sum)
                species_info[tmp_species]["gene_count"].append(tmp_gene_count)
                species_info[tmp_species]["reads_count"].append(int(filter_df.loc[i, "reads_count"]))
                down_gene_dict = self.estimate_saturation(gene_count_array = np.array(tmp_list), down_list = self.down_list)
                for j in down_gene_dict.keys():
                    species_info[tmp_species]["down_sample_gene_count"][j].append(down_gene_dict[j])
        self.mobilogger._mobilogrecorder(log_message="├── Processing barcode-rank plot is done successfully.", log_level="INFO")
        ###saturation plot
        self.mobilogger._mobilogrecorder(log_message="├── Processing saturation plot...", log_level="INFO")
        sat_data = pd.read_csv(self.saturation_file, sep="\t")
        sat_data.index = range(sat_data.shape[0])
        try:
            sat_plot_UMI = self.o_json["all_plot"]["saturation_plot_UMI"]
            sat_plot_gene = self.o_json["all_plot"]["saturation_plot_GENE"]
        except KeyError:
            sat_plot = {}
            sat_plot["config"] = {
                                    "displaylogo": 'false',
                                    
                                    "staticPlot": 'false',
                                    "displayModeBar": 'true',
                                    "showAxisDragHandles": 'true',
                                    "toImageButtonOptions": {},
                                    "scrollZoom": 'true'
                                }
            sat_plot["layout"] = {
                                    "title": {
                                        "text": "Sequencing Saturation Plot based on UMI",
                                        "font":  {
                                            "size": 16
                                        }
                                    },
                                    "showlegend": 'false',
                                    "hovermode": "closest",
                                    "xaxis": {
                                        "title": "Downsample ratio",
                                        "fixedrange": 'false',
                                        "titlefont": {
                                            "size": 14
                                        }
                                    },
                                    "yaxis": {
                                        "title": "Sequencing Saturation",
                                        "range": [
                                            0,
                                            1
                                        ],
                                        "fixedrange": 'false',
                                        "titlefont": {
                                            "size": 14
                                        }
                                    },
                                    "shapes": [
                                        {
                                            "type": "line",
                                            "x0": 0,
                                            "y0": 0.9,
                                            "x1": 1,
                                            "y1": 0.9,
                                            "line": {
                                                "color": "#dddddd",
                                                "width": 4,
                                                "dash": "dot"
                                            }
                                        }
                                    ],
                                    "margin": {
                                        "t" : 80
                                    }
                                }
            sat_plot["data"] = []
            sat_plot_gene = copy.deepcopy(sat_plot)
            sat_plot_UMI = copy.deepcopy(sat_plot)
            tmp_sat_plot = {}
            found_sat = False
            td_sat_plot = {}
            td_sat_plot['x'] = [0]
            td_sat_plot['y'] = [0]
            gene_sat_plot = {}
            gene_sat_plot['x'] = [0]
            gene_sat_plot['y'] = [0]
            for j in range(len(sat_data.columns.tolist())):
                if found_sat:
                    tmp_sat = 0
                    tmp_sum = 0
                    for i in sat_data.index.tolist():
                        if i != "All":
                            tmp_sat += sat_data.loc[i, "mapped_reads_count"] * sat_data.loc[i, sat_data.columns.tolist()[j]]
                            tmp_sum += sat_data.loc[i, "mapped_reads_count"]
                            tmp_species = sat_data.loc[i, "species"]
                            if not tmp_species in tmp_sat_plot.keys():
                                tmp_sat_plot[tmp_species] = {}
                                tmp_sat_plot[tmp_species]['x'] = []
                                tmp_sat_plot[tmp_species]['y'] = []
                            try: 
                                tmp_sat_plot[tmp_species]['x'].append(float(sat_data.columns.tolist()[j]))
                                tmp_sat_plot[tmp_species]['y'].append(round(float(sat_data.loc[i, sat_data.columns.tolist()[j]]),4))
                            except ValueError:
                                pass
                    sat_data.loc["All", sat_data.columns.tolist()[j]] = tmp_sat / tmp_sum
                    input_info["saturation_rate"] = tmp_sat / tmp_sum
                    try:
                        td_sat_plot['x'].append(float(sat_data.columns.tolist()[j]))
                        td_sat_plot['y'].append(round(float(sat_data.loc["All", sat_data.columns.tolist()[j]]),4))
                    except ValueError:
                        pass
                if sat_data.columns.tolist()[j] == "uniquely_mapped_UMI_count":
                    found_sat = True
            #gene saturation
            for i in cell_stat_data.columns.tolist():
                if "gene_saturation_" in i:
                    tmp_level = float(i.replace("gene_saturation_",""))
                    tmp_value = np.median(cell_stat_data.loc[:,i])
                    gene_sat_plot['x'].append(tmp_level)
                    gene_sat_plot['y'].append(tmp_value)
            if len(gene_sat_plot["x"]) > 1:
                ###sort
                sorted_indices = sorted(range(len(gene_sat_plot['x'])), key=lambda i: gene_sat_plot['x'][i])
                gene_sat_plot['x'] = [gene_sat_plot['x'][i] for i in sorted_indices]
                gene_sat_plot['y'] = [gene_sat_plot['y'][i] for i in sorted_indices]
                sat_plot_gene["layout"]["yaxis"]={
                                        "title": "Median Gene Count",
                                        "fixedrange": 'false',
                                        "titlefont": {
                                            "size": 14
                                        }
                                    }
                sat_plot_gene["layout"]["title"]={
                                        "text": "Sequencing Saturation Plot based on Median Gene Count",
                                        "font":  {
                                            "size": 16
                                        }
                                    }
                sat_plot_gene["data"].append(gene_sat_plot)
            sat_plot_UMI["data"].append(td_sat_plot)
        self.mobilogger._mobilogrecorder(log_message="├── Processing saturation plot is done successfully.", log_level="INFO")
        ###plots relate to species
        std_config = {
                        "displaylogo": 'false',
                        "staticPlot": 'false',
                        "displayModeBar": 'true',
                        "showAxisDragHandles": 'true',
                        "toImageButtonOptions": {},
                        "scrollZoom": 'true'
                    }
        std_layout = {
                        "title": {
                            "text": "Mitochondrial Fraction",
                            "font":  {
                                "size": 16
                            }
                        },
                        "hovermode": "closest",
                        "showlegend": 'false',
                        "yaxis": {
                            "rangemode": 'tozero', 'autorange': 'true'
                        }
                    }
        std_data = {'y':[], "type": "violin", "color":[], 
                    "name": "Ecoli", "hoverinfo": "y", 
                    "points": 'false', "spanmode": "manual", 
                    "hoveron": "violins",
                    "meanline" : {"visible": 'true'}, "span" : ["0", "y.max()"], 
                    }
        abundance_plot = {}
        abundance_plot["config"] = std_config.copy()
        abundance_plot["layout"] = std_layout.copy()
        abundance_plot["layout"]["showlegend"] = "true"
        abundance_plot["layout"]["title"] = "Abundance Distribution"
        abundance_plot["data"] = []
        abundance_data_dict = {}
        UMI_counts_plot = {}
        UMI_counts_plot["config"] = std_config.copy()
        UMI_counts_plot["layout"] = std_layout.copy()
        UMI_counts_plot["layout"]["title"] = "UMI Count Distribution"
        UMI_counts_plot["data"] = []
        UMI_counts_data_dict = {}
        gene_counts_plot = {}
        gene_counts_plot["config"] = std_config.copy()
        gene_counts_plot["layout"] = std_layout.copy()
        gene_counts_plot["layout"]["title"] = "Gene Count Distribution"
        gene_counts_plot["data"] = []
        gene_counts_data_dict = {}
        all_UMI_list = []
        all_gene_list = []
        all_reads_count = 0
        abundance_tmp_data = {'values':[], 'labels':[], "type": "pie", 'marker':{"colors":[]}}
        self.mobilogger._mobilogrecorder(log_message="├── Processing mapping info...", log_level="INFO")
        for i in species_info.keys():
            if i != "all":
                try:
                    sat_index = sat_data.loc[:,"species"].tolist().index(i)
                except ValueError:
                    pass
                else:
                    if i != "unknown":
                        species_info[i]["map_rate"] = sat_data.loc[sat_index, "mapped_reads_count"] / input_info["star_input_reads"]
                        species_info[i]["confidently_map_rate"] = sat_data.loc[sat_index, "uniquely_mapped_reads_count"] / input_info["star_input_reads"]
                if len(species_info[i]["gene_count"]) == 0:
                    species_info[i]["median_gene_count"] = 0
                else:
                    species_info[i]["median_gene_count"] = int(np.median(species_info[i]["gene_count"]))
                species_info[i]["cell_number"] = len(species_info[i]["UMI_count"])
                if len(species_info[i]["UMI_count"]) == 0:
                    species_info[i]["median_UMI_count"] = 0
                else:
                    species_info[i]["median_UMI_count"] = int(np.median(species_info[i]["UMI_count"]))
                if not i in abundance_data_dict.keys() and species_info[i]["cell_number"] != 0:
                    abundance_tmp_data["values"] += [species_info[i]["cell_number"]]
                    if not species_info[i]["line_color"] in abundance_tmp_data["marker"]["colors"]:
                        abundance_tmp_data["marker"]["colors"].append(species_info[i]["line_color"])
                    if i != "unknown":
                        abundance_tmp_data["labels"] += [i]
                    else:
                        abundance_tmp_data["labels"] += ["Multiplet"]
                    abundance_data_dict[i] = species_info[i]["cell_number"]
                if i != "unknown":
                    if not i in UMI_counts_data_dict.keys():
                        tmp_data = std_data.copy()
                        tmp_data["name"] = i
                        tmp_data["fillcolor"] = species_info[i]["fill_color"]
                        tmp_data["line"] = {"color":species_info[i]["line_color"]}
                        tmp_data["y"] = list(map(int, species_info[i]["UMI_count"]))
                        UMI_counts_data_dict[i] = len(UMI_counts_plot["data"])
                        UMI_counts_plot["data"].append(tmp_data)
                    if not i in gene_counts_data_dict.keys():
                        tmp_data = std_data.copy()
                        tmp_data["name"] = i
                        tmp_data["fillcolor"] = species_info[i]["fill_color"]
                        tmp_data["line"] = {"color":species_info[i]["line_color"]}
                        tmp_data["y"] = list(map(int, species_info[i]["gene_count"]))
                        gene_counts_data_dict[i] = len(gene_counts_plot["data"])
                        gene_counts_plot["data"].append(tmp_data)
                all_UMI_list += species_info[i]["UMI_count"]
                all_gene_list += species_info[i]["gene_count"]
                all_reads_count += sum(species_info[i]["reads_count"])
        for i in species_info.keys():
            if i != "all":
                if species_info[i]["cell_number"] != 0:
                    rated_reads = input_info["valid_reads"] * sum(species_info[i]["reads_count"]) / all_reads_count
                    species_info[i]["mean_reads_count"] = rated_reads / species_info[i]["cell_number"]
                    species_info[i]["mean_keep_reads_count"] = sum(species_info[i]["reads_count"]) / species_info[i]["cell_number"]
                else:    
                    species_info[i]["mean_reads_count"] = 0
                    species_info[i]["mean_keep_reads_count"] = 0
        species_info["all"]["median_gene_count"] = int(np.median(all_gene_list))
        species_info["all"]["median_UMI_count"] = int(np.median(all_UMI_list))
        species_info["all"]["mean_keep_reads_count"] = all_reads_count / species_info["all"]["cell_number"]
        species_info["all"]["mean_reads_count"] = input_info["valid_reads"] / species_info["all"]["cell_number"]
        if "species_info" in self.o_json.keys():
            species_info["all"]["map_rate"] = self.o_json["species_info"]["all"]["map_rate"]
            species_info["all"]["confidently_map_rate"] = self.o_json["species_info"]["all"]["confidently_map_rate"]
        else:
            all_mapped_reads_count = 0
            all_mapped_unique_reads_count = 0
            for i in sat_data.index.tolist():
                if not pd.isna(sat_data.loc[i, "mapped_reads_count"]):
                    all_mapped_reads_count += sat_data.loc[i, "mapped_reads_count"]
                if not pd.isna(sat_data.loc[i, "uniquely_mapped_reads_count"]):
                    all_mapped_unique_reads_count += sat_data.loc[i, "uniquely_mapped_reads_count"]
            species_info["all"]["map_rate"] = all_mapped_reads_count / input_info["star_input_reads"]
            species_info["all"]["confidently_map_rate"] = all_mapped_unique_reads_count / input_info["star_input_reads"]
        abundance_plot["data"].append(abundance_tmp_data)
        ###change low abundance species to other
        for i in range(len(abundance_plot["data"][0]["values"])):
            if abundance_plot["data"][0]["values"][i] / species_info["all"]["cell_number"] < 0.05 and \
                abundance_plot["data"][0]["labels"][i] != "Multiplet":
                abundance_plot["data"][0]["labels"][i] = "Other_Species"
        filter_data.obs['n_count'] = filter_data.X.sum(axis=1)
        cell_stat_data.index = range(cell_stat_data.shape[0])
        self.mobilogger._mobilogrecorder(log_message="├── Processing mapping info is done successfully.", log_level="INFO")
        if self.run_cluster:
            with warnings.catch_warnings():
                filter_data = process_secondary_analysis(filter_data=filter_data.copy(), 
                                                 species_info=species_info.copy(), 
                                                 cell_stat_data=cell_stat_data, 
                                                 output_path=self.output_path, 
                                                 threads=self.threads, 
                                                 mobilogger=self.mobilogger, 
                                                 h5ad_out_file=os.path.join(self.output_path, "combined.h5ad"))
                UMAP_count_plot = self.get_UMAP_plot_count(filter_data)
                UMAP_cluster_plot = self.get_UMAP_plot_cluster(filter_data)
                filter_data.obs_names = filter_data.obs.loc[:, "barcode"]
                exp_barcodes = filter_data.obs_names.tolist()
                #exp_raw = filter_data[exp_barcodes].copy()
                if os.path.exists(os.path.join(self.filter_mtx_path, "features.tsv")):
                    tmp_gene_df = os.path.join(self.filter_mtx_path, "features.tsv")
                elif os.path.exists(os.path.join(self.filter_mtx_path, "features.tsv.gz")):
                    tmp_gene_df = os.path.join(self.filter_mtx_path, "features.tsv.gz")
                else:
                    self.mobilogger._mobilogrecorder(log_message="├── Feature.tsv not found in %s." %(self.filter_mtx_path), log_level="ERROR")
                    sys.exit()
                self.mobilogger._mobilogrecorder(log_message="└── Secondary analysis is done successfully.", log_level="INFO")
        else:
            UMAP_count_plot = dict()
            UMAP_cluster_plot = dict()
        ###ID
        new_textwrapID = '\n'.join([i.center(24) for i in textwrap.wrap(self.sample_ID, 24)])
        ###output json
        o_json = self.o_json
        ###sample info
        if not "sample" in o_json.keys():
            o_json["sample"] = {}
            o_json["sample"]["id"] = new_textwrapID
            o_json["sample"]["id_ori"] = self.sample_ID
            o_json["sample"]["description"] = 'Single-microbe RNA Sequencing Analysis Report'
            o_json["sample"]["reference"] = self.reference.split("/")[-2]
        o_json["sample"]["run_cmd"] = self.run_cmd
        o_json["sample"]["run_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        o_json["sample"]["kit"] = self.__vers_kit
        o_json["sample"]["version"] = self.version.split("_")[0]
        o_json["sample"]["version_detail"] = self.version
        ###species info
        o_json["species_info"] = species_info
        ###plot info
        o_json["all_plot"] = {}
        o_json["all_plot"]["barcode_rank_plot"] = rank_plot
        o_json["all_plot"]["saturation_plot_UMI"] = sat_plot_UMI
        o_json["all_plot"]["saturation_plot_GENE"] = sat_plot_gene
        o_json["all_plot"]["UMI_counts"] = UMI_counts_plot
        o_json["all_plot"]["gene_counts"] = gene_counts_plot
        o_json["all_plot"]["abundance_plot"] = abundance_plot
        o_json["all_plot"]["RNA_type_plot"] = type_plot
        o_json["all_plot"]["UMAP_count"] = UMAP_count_plot
        o_json["all_plot"]["UMAP_cluster"] = UMAP_cluster_plot
        json_file = self.output_path + "/report.json"
        o_json = sanitize_keys(o_json)
        with open(json_file, "w") as outfile:
            json.dump(o_json,outfile, indent=4, cls=NpEncoder)
        return json_file

if __name__ == "__main__":
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('-o', '--output_path', type=str, default="summary_data", help='The output path.')
    parse.add_argument('-ID', '--sample_ID', type=str, help='The ID of the sample')
    parse.add_argument('-r', '--reference_json_file', type=str, help='The path of json file of reference.')
    parse.add_argument('--pre_stat_file', type=str, help="The file of raw stat.")
    parse.add_argument('--input_stat_file', type=str, help="The file of pre-process stat.")
    parse.add_argument('--filter_stats_file', type=str, help="The file of filter stats.")
    parse.add_argument('--saturation_file', type=str, help='The path of seq saturation file.')
    parse.add_argument('--raw_mtx_path', type=str, help="The paht of raw mtx.")
    parse.add_argument('--filter_mtx_path', type=str, help="The path of filtered mtx.")
    parse.add_argument('--cell_stat_file', type=str, help="The path of cell stat file")
    parse.add_argument('--kit', type=str, help="The version of the kit.")
    parse.add_argument('--threads', type=int, help="The run threads.")
    parse.add_argument('--run_cmd', type=str, help="The run cmd of all analysis.")
    args = parse.parse_args()
    Exp_test = Data_InjectTool(sample_ID=args.sample_ID, 
                               output_path=args.output_path,
                               pre_stat_file=args.pre_stat_file,
                               input_stat_file=args.input_stat_file,
                               filter_stats_file = args.filter_stats_file,
                               saturation_file=args.saturation_file, 
                               raw_mtx_path=args.raw_mtx_path, 
                               filter_mtx_path=args.filter_mtx_path,
                               cell_stat_file=args.cell_stat_file, 
                               reference_json_path=args.reference_json_file, 
                               kit=args.kit, 
                               run_cmd=args.run_cmd, 
                               threads=args.threads)
    Exp_test.export_json_microbe()
    
