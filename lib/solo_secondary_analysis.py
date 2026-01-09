import multiprocessing as mp
mp.set_start_method('spawn', force=True)
import pandas as pd
import json
import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import colormaps
import sys
sys.path.append(os.path.dirname(__file__))
from dataexp import process_secondary_analysis, readin_custom_mtx
from mobivisionlogging import MobiLoggingSystem

colormaps.get_cmap = plt.get_cmap

if __name__ == "__main__":
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', '--input_path', type=str, help='The input path.')
    parse.add_argument('-o', '--output_path', type=str, help='The output path')
    parse.add_argument('-t', '--threads', type=int, help="The run threads", default=8)
    parse.add_argument('-r', '--ref_json', type=str, help="The path of reference json")
    parse.add_argument('-m', '--sim_method', type=str, help="The method of gene expression similarity. Valid values are pearson, spearman and cosine.", default="pearson")
    parse.add_argument('-p', '--pattern', type=str, help="The regular expression to filter gene", default=None)
    parse.add_argument('--barcode_info', type=str, help="Optional barcode info file", default=None)
    parse.add_argument('--species_count', type=int, help="Optional number of species", default=1)
    args = parse.parse_args()
    print("All input argments are:")
    for k,v in vars(args).items():
        print(k,'=',v)
    if os.path.exists(args.output_path):
        print("Won't overwrite the existen path. %s" %(args.output_path))
    else:
        os.makedirs(args.output_path)
    mobilogger = MobiLoggingSystem(o_dir=args.output_path, dev_mode=False)
    if not os.path.exists(args.input_path):
        mobilogger._mobilogrecorder(log_message="The input path not found. %s" %(args.input_path), log_level="ERROR")
        sys.exit()
    if not os.path.exists(args.ref_json):
        mobilogger._mobilogrecorder(log_message="The reference not found. %s" %(args.input_path), log_level="ERROR")
        sys.exit()
    if args.barcode_info is not None:
        if not os.path.exists(args.barcode_info):
            mobilogger._mobilogrecorder(log_message="The barcode info file not found. %s" %(args.barcode_info), log_level="ERROR")
            sys.exit()
        else:
            barcode_df = pd.read_csv(args.barcode_info, sep="\t")
    else:
        barcode_df = None
    species_info = {}
    for i in range(args.species_count):
        species_info[i] = {}
    with open(args.ref_json,'r') as load_f:
        ref_dict = json.load(load_f)
    adatas = []
    for i in os.listdir(args.input_path):
        try_path = os.path.join(args.input_path, i, "filtered_feature_bc_matrix")
        if os.path.exists(try_path):
            tmp_mtx_file = None
            tmp_feature_file = None
            tmp_barcode_file = None
            for j in os.listdir(try_path):
                if "matrix.mtx" in j:
                    tmp_mtx_file = j
                elif "features.tsv" in j:
                    tmp_feature_file = j
                elif "barcodes.tsv" in j:
                    tmp_barcode_file = j
            tmp_filter_data, gene_names, tmp_output_type = readin_custom_mtx(mobilogger=mobilogger, 
                                                                            i_dir=try_path, 
                                                                            mtx_file=tmp_mtx_file, 
                                                                            feature_file=tmp_feature_file, 
                                                                            barcode_file=tmp_barcode_file, 
                                                                            ref=ref_dict, 
                                                                            input_type="filtered")
            tmp_filter_data.obs['sample'] = i
            tmp_filter_data.obs['n_count'] = tmp_filter_data.X.sum(axis=1)
            adatas.append(tmp_filter_data)
    adata = sc.concat(adatas, axis=0,
                  join='outer',
                  merge='same',
                  label='sample',
                  keys=[a.obs['sample'][0] for a in adatas])
    adata.var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    processed_adata = process_secondary_analysis(filter_data=adata, 
                                                 species_info=species_info, 
                                                 cell_stat_data=barcode_df, 
                                                 output_path=args.output_path, 
                                                 threads=args.threads, 
                                                 mobilogger=mobilogger, 
                                                 sim_method=args.sim_method, 
                                                 gene_filter=args.pattern, 
                                                 h5ad_out_file=os.path.join(args.output_path, "combined.h5ad"))
    fig, ax = plt.subplots(figsize=(5,4), dpi=300)
    sc.pl.umap(processed_adata, color=['leiden'], ax=ax, show=False, frameon=False)
    fig.savefig(os.path.join(args.output_path, "combine_umap_cluster.png"), bbox_inches='tight', dpi=300)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5,4), dpi=300)
    sc.pl.umap(processed_adata, color=['sample'], ax=ax, show=False, frameon=False)
    fig.savefig(os.path.join(args.output_path, "combine_umap_sample.png"), bbox_inches='tight', dpi=300)
    plt.close(fig)
    hvgs = processed_adata.var_names[processed_adata.var['highly_variable']]
    feature_dir = os.path.join(args.output_path, "feature_plots")
    os.makedirs(feature_dir)
    for g in hvgs:
        fig = sc.pl.umap(processed_adata,
                        color=g,
                        frameon=False,
                        s=20,
                        title=g,
                        show=False,
                        return_fig=True) 
        fig.savefig(os.path.join(feature_dir, "%s_feature_plot.png" %(g)), dpi=200, bbox_inches='tight')
        plt.close(fig)
    sc.tl.rank_genes_groups(
        processed_adata,
        groupby='leiden',
        method='wilcoxon',
        key_added='rank_genes_groups',
        use_raw=False
    )
    gene_df = sc.get.rank_genes_groups_df(
        processed_adata,
        group=None,
        key='rank_genes_groups',
        pval_cutoff=0.05,
        log2fc_min=0
    )
    gene_df.to_csv(os.path.join(args.output_path, "cluster_feature_gene.tsv"), sep="\t")
    sc.settings.figdir = args.output_path
    sc.pl.rank_genes_groups_heatmap(
        processed_adata,
        n_genes=10,
        groups='all',
        key='rank_genes_groups',
        use_raw=False,
        show_gene_labels=True,
        swap_axes=False,
        dendrogram=False,
        figsize=(10,6),
        save='feature_heatmap'
    )
    
