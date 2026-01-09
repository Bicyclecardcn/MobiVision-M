import pandas as pd
import os 
import numpy as np
import math
import multiprocessing
import json
import sys
import gzip 
import shutil
sys.path.append(os.path.dirname(__file__))
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor


def softmax(num_list: list, base=2):
    for i in range(len(num_list)):
        if num_list[i] > 0:
            num_list[i] = math.log(num_list[i], base)
        else:
            num_list[i] = -10
        num_list[i] = np.exp(num_list[i])
    temp_sum = sum(num_list)
    for i in range(len(num_list)):
        try:
            num_list[i] = num_list[i] / temp_sum
        except ZeroDivisionError:
            num_list[i] = 0
    return num_list

def get_alignBAM_param(map_result_dir):
    for root, dir, files in os.walk(map_result_dir):
        for fname in files:
            if '_Aligned.sortedByCoord.out.bam' in fname and not ".bai" in fname:
                bam_file_path = os.path.join(root, fname)
                return bam_file_path

def readin_coo_mtx(matrix_file:str, UMI_Count_dict:dict):
    output_dict = {}
    for i in UMI_Count_dict.keys():
        output_dict[i] = {}
    for i in output_dict.keys():
        output_dict[i]["Count"] = {}
        output_dict[i]["gene_Count"] = {}
        output_dict[i]["Total_Count"] = 0
    first_line = True
    with open(matrix_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not "%" in line:
                if not first_line:
                    line = line.replace("\n", "").split(" ")
                    tmp_row = int(line[0])
                    tmp_col = int(line[1])
                    for j in UMI_Count_dict.keys():
                        if tmp_row >= UMI_Count_dict[j]["start"] and tmp_row <= UMI_Count_dict[j]["start"] + UMI_Count_dict[j]["len"]:
                            break
                    if not tmp_col in output_dict[j]["Count"].keys():
                        output_dict[j]["Count"][tmp_col] = float(line[2])
                        output_dict[j]["gene_Count"][tmp_col] = 1
                    else:
                        output_dict[j]["Count"][tmp_col] += float(line[2])
                        output_dict[j]["gene_Count"][tmp_col] += 1
                    output_dict[j]["Total_Count"] += float(line[2])
                else:
                    first_line = False
    return output_dict

def count_UMI(matrix_info, UMI_Count_dict, barcode_list, i_list):
    return_dict = {}
    for i in i_list:
        ###UMI
        UMI_temp_Count = np.sum(matrix_info[UMI_Count_dict[i]["start"]:UMI_Count_dict[i]["start"] + UMI_Count_dict[i]["len"]],axis=0)
        UMI_temp_Count = UMI_temp_Count.reshape(-1,1)
        UMI_temp_Count =pd.DataFrame(UMI_temp_Count)
        #UMI_Count_dict[i]["Count"] = UMI_temp_Count
        ###gene Count
        Gene_temp_Count = []
        for j in range(len(barcode_list)):
            temp_Gene_Count = np.count_nonzero(matrix_info[UMI_Count_dict[i]["start"]:UMI_Count_dict[i]["start"] + UMI_Count_dict[i]["len"],j])
            Gene_temp_Count.append(temp_Gene_Count)
        UMI_Count_dict[i]["gene_Count"] = Gene_temp_Count
        return_dict[i] = [UMI_temp_Count, Gene_temp_Count]
    return return_dict, i_list

def find_best(UMI_Count_dict, Total_Data_info, j_list):
    return_dict = {}
    for j in j_list:
        max_UMI_count = -1
        best_species = ""
        for i in UMI_Count_dict.keys():
            if Total_Data_info.loc[j, "UMI_" + str(i) + "_Count"] > max_UMI_count:
                max_UMI_count = Total_Data_info.loc[j, "UMI_" + str(i) + "_Count"]
                best_species = str(i)
        return_dict[j] = [max_UMI_count, best_species]
    return return_dict, j_list

def apply_df(Temp_Data_info:pd.DataFrame, UMI_Count_dict:dict, i_list:list):
    for i in i_list:
        Temp_Data_info.loc[i, "UMI_Total_Count"] = 0
        for j in UMI_Count_dict.keys():
            temp_col = "UMI_" + str(j) + "_Count"
            try:
                Temp_Data_info.loc[i, temp_col] = UMI_Count_dict[j]["Count"][i+1]
                Temp_Data_info.loc[i, "UMI_Total_Count"] += UMI_Count_dict[j]["Count"][i+1]
            except KeyError:
                Temp_Data_info.loc[i, temp_col] = 0
            temp_col = "Gene_" + str(j) + "_Count"
            try:
                Temp_Data_info.loc[i, temp_col] = UMI_Count_dict[j]["gene_Count"][i+1]
            except KeyError:
                Temp_Data_info.loc[i, temp_col] = 0
    return Temp_Data_info

def fetch_bc_best(i_dir:str, o_dir:str, threads:int, species_list:list, method:str):
    feature_file = i_dir + "/features.tsv"
    barcode_file = i_dir + "/barcodes.tsv"
    matrix_file = i_dir + "/UniqueAndMult-Uniform.mtx"
    features_list = pd.read_csv(feature_file,names = ["features1", "features2", "features3"], sep='\t')
    barcode_list = pd.read_csv(barcode_file,names = ["barcodes"])
    barcode_list = list(barcode_list['barcodes'])
    features_list = list(features_list['features2'])
    UMI_Count_dict = {}
    for i in range(len(features_list)):
        temp_spe_name = features_list[i].split("_")[0]
        if not temp_spe_name in species_list:
            temp_spe_name = species_list[0]
        if not temp_spe_name in UMI_Count_dict.keys():
            UMI_Count_dict[temp_spe_name] = {}
            UMI_Count_dict[temp_spe_name]["start"] = i
            UMI_Count_dict[temp_spe_name]["len"] = 1
        else:
            UMI_Count_dict[temp_spe_name]["len"] += 1
    UMI_Count_dict = readin_coo_mtx(matrix_file=matrix_file, UMI_Count_dict=UMI_Count_dict)
    apply_columns = ["UMI_Total_Count"]
    for j in UMI_Count_dict.keys():
        temp_col = "UMI_" + str(j) + "_Count"
        apply_columns.append(temp_col)
        temp_col = "Gene_" + str(j) + "_Count"
        apply_columns.append(temp_col)
    Total_Data_info = pd.DataFrame(columns=['barcode_list'] + apply_columns)
    Total_Data_info.loc[:,'barcode_list'] = barcode_list
    Total_Data_info.index = range(Total_Data_info.shape[0])
    if method == "Uniform":
        all_task = []
        pool = multiprocessing.Pool(processes = threads)
        i_list = []
        for i in Total_Data_info.index.tolist():
            i_list.append(i)
            if len(i_list) >= 1e2:
                tmp_args = [Total_Data_info, UMI_Count_dict, i_list]
                all_task.append(pool.apply_async(apply_df, tmp_args))
        pool.close()
        for res in all_task:
            stat = res.get()
            for i in stat.index:
                Total_Data_info.loc[i, ] = stat.loc[i, ]
        pool.join()
        try:
            Total_Data_info.loc[:, 'MAX_UMI_Count'] = 0
        except ValueError:
            print("ERROR", Total_Data_info)
            print(barcode_list)
            Total_Data_info.loc[:, 'MAX_UMI_Count'] = 0
        pool = multiprocessing.Pool(processes = threads)
        all_task = []
        j_list = []
        for j in Total_Data_info.index.tolist():
            j_list.append(j)
            if len(j_list) > 1e2:
                tmp_args = [UMI_Count_dict, Total_Data_info, j_list]
                all_task.append(pool.apply_async(find_best, args=tmp_args))
                j_list = []
        if len(j_list) > 0:
            tmp_args = [UMI_Count_dict, Total_Data_info, j_list]
            all_task.append(pool.apply_async(find_best, args=tmp_args))
            j_list = []
        pool.close()
        for res in all_task:
            stat = res.get()
            tmp_j_list = stat[1]
            for j in tmp_j_list:
                Total_Data_info.loc[j, "MAX_UMI_Count"] = stat[0][j][0]
                Total_Data_info.loc[j, "best_species"] = stat[0][j][1]
        pool.join()
        #for i in Total_Data_info.index.tolist():
        #    species_dict[Total_Data_info.loc[i, "barcode_list"]] = Total_Data_info.loc[i, "best_species"]
        #try:
        #    Total_Data_info['Cell_Ratio'] = Total_Data_info['MAX_UMI_Count'] / Total_Data_info['UMI_Total_Count']
        #except ZeroDivisionError:
        #    Total_Data_info['Cell_Ratio'] = 0
        Total_Data_info.loc[:, "sum_UMI_Count"] = 0
        for i in Total_Data_info.index:
            softmax_list = []
            for j in UMI_Count_dict.keys():
                softmax_list.append(Total_Data_info.loc[i, "UMI_" + str(j) + "_Count"])
                Total_Data_info.loc[i, "sum_UMI_Count"] += Total_Data_info.loc[i, "UMI_" + str(j) + "_Count"]
            Total_Data_info.loc[i, "softmax"] = max(softmax(softmax_list))
    out_Total_Data_result_file = o_dir + "/pre_Total_cell_stats.csv"
    Total_Data_info.to_csv(out_Total_Data_result_file, index=False)
    return Total_Data_info.shape[0]


def call_cell_by_STAR(mtx_dir:str, tmp_dir:str, filter_args:str, filter_mtx_dir:str, mobiexecutor: CommandExecutor):
    if os.path.exists(tmp_dir):
        cmd = "rm -r %s " %(tmp_dir)
        os.system(cmd)
    for i in os.listdir(mtx_dir):
        if i.endswith(".gz"):
            cmd = "gunzip %s" %(os.path.join(mtx_dir, i))
            exit_code = mobiexecutor.execute(command=cmd, context={"call_cell": "gunzip"},console_output=False)
    tmp_output_path = os.path.dirname(filter_mtx_dir) + "/"
    cmd = "STAR --outTmpDir %s --runMode soloCellFiltering %s %s --soloCellFilter %s --outFileNamePrefix %s " %(tmp_dir, mtx_dir, filter_mtx_dir, filter_args, tmp_output_path)
    exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "STAR-call-cell"},console_output=False)
    if exit_code == 0:
        return_flag = "INFO"
        return_info = "Cell Calling by STAR is done successfully."
    else:
        return_flag = "ERROR"
        return_info = "Cell Calling by STAR failed."
    if os.path.exists(tmp_dir):
        cmd = "rm -r %s " %(tmp_dir)
        exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "STAR-call-cell"},console_output=False)
    for i in os.listdir(filter_mtx_dir):
        if not i.endswith(".gz"):
            if not i.endswith(".tsv") and not i.endswith(".mtx"):
                os.remove(os.path.join(filter_mtx_dir, i))
            else:
                cmd = "gzip %s " %(os.path.join(filter_mtx_dir, i))
                exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "gzip"},console_output=False)
    return filter_mtx_dir, return_flag, return_info

def check_coo_line(check_line:list, valid_index:list):
    return_list = []
    n = 0
    for split_line in check_line:
        tmp_index = split_line[1]
        if int(tmp_index) in valid_index:
            write_index = valid_index.index(int(tmp_index)) +1
            write_line = " ".join([split_line[0], str(write_index), split_line[2]])
            return_list.append(write_line.encode('utf-8'))
            n += 1
    return return_list, n

def call_cell_by_threshold(mtx_dir:str, cell_stat_file:str, filter_args:str, filter_mtx_dir:str, threads:int, mobiexecutor: CommandExecutor):
    if not os.path.exists(filter_mtx_dir):
        os.makedirs(filter_mtx_dir)
    ###filter by UMI or genes
    filter_flag = filter_args.split(":")[0]
    filter_value = float(filter_args.split(":")[1])
    if cell_stat_file.endswith(".gz"):
        cell_info_df = pd.read_csv(cell_stat_file, sep="\t", compression='gzip')
    else:
        cell_info_df = pd.read_csv(cell_stat_file, sep="\t")
    if filter_flag == "min_UMI":
        cell_info_df = cell_info_df[cell_info_df["gene_UMI_count"] >= filter_value]
    elif filter_flag == "min_reads":
        cell_info_df = cell_info_df[cell_info_df["gene_read_count"] >= filter_value]
    else:
        print("No filter found!")
        print(filter_args)
        sys.exit()
    valid_barcode = cell_info_df.loc[:,"barcode"].tolist()
    ###export filtered barcode info
    valid_index = []
    if os.path.exists(os.path.join(mtx_dir, "barcodes.tsv.gz")):
        with gzip.open(os.path.join(mtx_dir, "barcodes.tsv.gz"), "r") as f_i, \
        gzip.open(os.path.join(filter_mtx_dir, "barcodes.tsv.gz"), "w") as f_o:
            n = 1
            for line in f_i:
                tmp_barcode = line.decode('utf-8').replace("\n","")
                if tmp_barcode in valid_barcode:
                    valid_index.append(n)
                    f_o.write(tmp_barcode.encode('utf-8') + b"\n")
                n += 1
    else:
        with open(os.path.join(mtx_dir, "barcodes.tsv"), "r") as f_i, \
        gzip.open(os.path.join(filter_mtx_dir, "barcodes.tsv.gz"), "w") as f_o:
            n = 1
            for line in f_i:
                tmp_barcode = line.replace("\n","")
                if tmp_barcode in valid_barcode:
                    valid_index.append(n)
                    f_o.write(tmp_barcode.encode('utf-8') + b"\n")
                n += 1
    ###cp gene info
    if os.path.exists(os.path.join(mtx_dir, "features.tsv.gz")):
        shutil.copyfile(os.path.join(mtx_dir, "features.tsv.gz"), os.path.join(filter_mtx_dir, "features.tsv.gz"))
    else:
        shutil.copyfile(os.path.join(mtx_dir, "features.tsv"), os.path.join(filter_mtx_dir, "features.tsv"))
        cmd = "gzip %s " %(os.path.join(filter_mtx_dir, "features.tsv"))
        exit_code = mobiexecutor.execute(command=cmd, context={"call_cell": "gzip"},console_output=False)
    n = 0
    head_line = ""
    ###export filtered mtx
    if os.path.exists(os.path.join(mtx_dir, "matrix.mtx.gz")):
        f_i = gzip.open(os.path.join(mtx_dir, "matrix.mtx.gz"), "rt")
    else:
        f_i = open(os.path.join(mtx_dir, "matrix.mtx"), "r")
    f_o = gzip.open(os.path.join(filter_mtx_dir, "matrix_tmp.mtx.gz"), "wb")
    if not isinstance(threads, int):
        threads = int(threads)
    pool = multiprocessing.Pool(processes = threads)
    all_task = []
    apply_line = []
    for line in f_i:
        decode_line = line
        if not decode_line.startswith("%"):
            n += 1
            split_line = decode_line.split(" ")
            if n == 1:
                head_line = " ".join([split_line[0], str(len(valid_index)), split_line[2]])
                f_o.write(head_line.encode('utf-8'))
            else:
                apply_line.append(split_line)
                if len(apply_line) >= 500:
                    tmp_args = [apply_line, valid_index]
                    all_task.append(pool.apply_async(check_coo_line, tmp_args))
                    apply_line = []
        else:
            if isinstance(line, str):
                line = line.encode('utf-8')
            f_o.write(line)
    if len(apply_line) > 0:
        tmp_args = [apply_line, valid_index]
        all_task.append(pool.apply_async(check_coo_line, tmp_args))
        apply_line = []
    pool.close()
    n_line = 0
    for res in all_task:
        stat = res.get()
        for line in stat[0]:
            if isinstance(line, str):
                line = line.encode('utf-8')
            f_o.write(line)
            n_line += 1
    pool.join()
    f_i.close()
    f_o.close()
    n = 0
    with gzip.open(os.path.join(filter_mtx_dir, "matrix_tmp.mtx.gz"), "rt") as f_i, \
    gzip.open(os.path.join(filter_mtx_dir, "matrix.mtx.gz"), "wb") as f_o:
        for line in f_i:
            decode_line = line
            if not decode_line.startswith("%"):
                n += 1
                split_line = decode_line.split(" ")
                if n == 1:
                    head_line = " ".join([split_line[0], split_line[1], str(n_line)])
                    f_o.write(head_line.encode('utf-8') + b"\n")
                else:
                    if isinstance(line, str):
                        line = line.encode('utf-8')
                    f_o.write(line)
            else:
                if isinstance(line, str):
                    line = line.encode('utf-8')
                f_o.write(line)
    #os.remove(os.path.join(filter_mtx_dir, "matrix_tmp.mtx.gz"))
    return filter_mtx_dir

def re_assign_bam(map_result_dir:str, raw_mtx_dir:str, sample_id:str, 
                  split_bam_program:str, filter_args:str, threads:int, 
                  ref_json_path:str, read_number:int, dev_mod=True, 
                  allow_multi_target_UMI=False, reclaim_UMI=False, mobilogger=None, 
                  method="unique"):
    mobicommandlogger = MobiCommandLogSystem(o_dir=mobilogger.working_path, dev_mode=mobilogger.dev_mode)
    mobiexecutor = CommandExecutor(log_system=mobicommandlogger, console_output=False)
    with open(ref_json_path,'r') as load_f:
        ref = json.load(load_f)
    if sample_id == None:
        sample_id = ref["sample"]["id_ori"]
    species_list = ref['genomes']
    mtx_dir = raw_mtx_dir
    estimate_cell_number = fetch_bc_best(i_dir=mtx_dir , o_dir=map_result_dir, threads=threads, species_list=species_list, method="Species_Only")
    bam_file = get_alignBAM_param(map_result_dir)
    species_stat_file = map_result_dir + "/pre_Total_cell_stats.csv"
    #fetched_stat_file = map_result_dir + '/fetched_reads.tsv'
    alignment_stat_file = map_result_dir + '/map_stat.tsv'
    fetched_mtx_dir = map_result_dir + "/fetched_mtx"
    raw_feature_file = raw_mtx_dir + "/features.tsv"
    if not os.path.exists(fetched_mtx_dir):
        os.makedirs(fetched_mtx_dir)
    sort_name_bam = bam_file.replace("sortedByCoord.out", "sortedByName.out")
    cmd = "samtools sort -n -t UB -@ %s -o %s %s " %(str(threads), sort_name_bam, bam_file)
    #subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "samtools-sort"},console_output=False)
    if exit_code != 0:
        mobilogger._mobilogrecorder(log_message="Samtools sort failed.",
            log_level="ERROR")
        sys.exit()
    try_times = 0
    mobilogger._mobilogrecorder(log_message="Re-assigning multi-loci mapped reads...",
        log_level="INFO")
    while try_times < 5 and not os.path.exists(map_result_dir + "/map_stat.tsv"):
        #tmp_out = map_result_dir + "/tmp_go_%s.out" %(str(try_times))
        if os.path.exists(bam_file.replace(".bam","") + ".unique.unsort.bam"):
            cmd = "rm %s " %(bam_file.replace(".bam","") + ".unique.unsort.bam")
            os.system(cmd)
            exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "rm"},console_output=False)
            if exit_code != 0:
                mobilogger._mobilogrecorder(log_message="rm failed.",
                    log_level="ERROR")
                sys.exit()
        extra_args = ""
        if allow_multi_target_UMI:
            extra_args += "-Um "
        if reclaim_UMI:
            extra_args += "-re "
        with open(os.path.join(map_result_dir, "species_list.txt"), "w") as f:
            for i in species_list:
                f.write(i + "\n")
        cmd = " %s -CB=%s -NH=%s -s=%s -i=%s -f=%s -o=%s -os=%s -om=%s -t=%s -dg=%s -sp=%s -rn=%s -cn=%s -m=%s %s" %(split_bam_program, "CB", "NH", species_stat_file, 
                                                                        sort_name_bam, raw_feature_file, bam_file.replace(".bam","") + ".unique.unsort.bam", map_result_dir, 
                                                                        fetched_mtx_dir, str(threads), "0.1", os.path.join(map_result_dir, "species_list.txt"), str(read_number), str(estimate_cell_number), "Species_Only", extra_args)
        #subprocess.call(cmd, stdout=f_o, stderr=f_o, shell=True)
        if dev_mod:
            exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "fetch-reads"},console_output=True)
        else:
            exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "fetch-reads"},console_output=False)
        try_times += 1
    if exit_code != 0:
        mobilogger._mobilogrecorder(log_message="Re-assigning multi-loci mapped reads failed.",
            log_level="ERROR")
        sys.exit()
    else:
        mobilogger._mobilogrecorder(log_message="Re-assigning multi-loci mapped reads is done successfully.",
            log_level="INFO")
    if dev_mod and os.path.exists(map_result_dir + "/map_stat.tsv"):
        for i in range(try_times):
            cmd = "rm %s " %(map_result_dir + "/tmp_go_%s.out" %(str(i)))
    cmd = "samtools sort -@ %s -o %s %s " %(str(threads), bam_file.replace(".bam","") + ".unique.bam", bam_file.replace(".bam","") + ".unique.unsort.bam")
    exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "samtools-sort"},console_output=False)
    cmd = "samtools index %s" %(bam_file.replace(".bam","") + ".unique.bam")
    exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "samtools-index"},console_output=False)
    cmd = "rm -r %s " %(bam_file.replace(".bam","") + ".unique.unsort.bam")
    exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "samtools-rm"},console_output=False)
    filter_mtx_dir = os.path.join(map_result_dir, sample_id + "_Solo.unqie_only.out/GeneFull/filtered/")
    if method != "Unique":
        fetched_mtx_dir = os.path.join(map_result_dir, sample_id+"_Solo.out","GeneFull","raw")
        if os.path.exists(os.path.join(fetched_mtx_dir, "matrix.mtx")):
            os.remove(os.path.join(fetched_mtx_dir, "matrix.mtx"))
        expect_mtx_name = "UniqueAndMult-%s.mtx" %(method)
        cmd = "mv %s %s " %(os.path.join(fetched_mtx_dir, expect_mtx_name), os.path.join(fetched_mtx_dir, "matrix.mtx"))
        os.system(cmd)
    if not "min_UMI" in filter_args and not "min_reads" in filter_args:
        if os.path.exists("_STARtmp"):
            cmd = "rm -r %s " %("_STARtmp")
            exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "rm"},console_output=False)
        if os.path.exists(os.path.abspath('..') + "/_STARtmp"):
            cmd = "rm -r %s " %(os.path.abspath('..') + "/_STARtmp")
            exit_code = mobiexecutor.execute(command=cmd, context={"fetch-multi": "rm"},console_output=False)
        tmp_path = map_result_dir + "/_STARtmp2"
        mobilogger._mobilogrecorder(log_message="Cell Calling by STAR...", log_level="INFO")
        filter_mtx_dir, info_flag, info = call_cell_by_STAR(mtx_dir=fetched_mtx_dir, tmp_dir=tmp_path, filter_args=filter_args, filter_mtx_dir=filter_mtx_dir, mobiexecutor=mobiexecutor)
        mobilogger._mobilogrecorder(log_message=info, log_level=info_flag)
        if info_flag == "ERROR":
            sys.exit()
    elif "min_UMI" in filter_args or "min_reads" in filter_args:
        cell_stat_file = os.path.join(map_result_dir, "cell_stat.tsv")
        mobilogger._mobilogrecorder(log_message="Cell Calling by hard filter...", log_level="INFO")
        filter_mtx_dir = call_cell_by_threshold(mtx_dir=fetched_mtx_dir, cell_stat_file=cell_stat_file, filter_args=filter_args, filter_mtx_dir=filter_mtx_dir, threads=threads, mobiexecutor=mobiexecutor)
        mobilogger._mobilogrecorder(log_message="Cell Calling by hard filter is done successfully.", log_level="INFO")
    else:
        mobilogger._mobilogrecorder(log_message="Invalid filter arguments %s " %(filter_args), log_level="ERROR")
        sys.exit()
    return bam_file + ".fetched.bam", alignment_stat_file, fetched_mtx_dir, filter_mtx_dir

class fetch_multi_reads:
    def __init__(self, map_result_dir, sample_id, split_bam_program, threads, filter_args, o_dir, ref_json_path, read_number, mobilogger=None):
        self.map_result_dir = map_result_dir
        self.sample_id = sample_id
        self.split_bam_program = split_bam_program
        self.threads = threads
        self.filter_args = filter_args
        self.o_dir = o_dir
        self.ref_json_path = ref_json_path
        self.read_number = read_number
        if mobilogger != None:
            self.mobilogger = mobilogger
        else:
            self.mobilogger = MobiLoggingSystem(o_dir=self.o_dir, dev_mode=False)

    #def re_call_mtx(self):
    #    tmp_path = os.path.exists(self.o_dir, "_STARtmp2")
    #    i_list = os.listdir(self.map_result_dir)
    #    found = False
    #    for i in i_list:
    #        expect_file = os.path.join(self.map_result_dir, i, "raw_cell_gene_matrix", "matrix.mtx.gz")
    #        if os.path.exists(expect_file):
    #            found = True
    #            break
    #    if not found:
    #        #prog_runlog(time.strftime("%Y-%m-%d %H:%M:%S\t", time.localtime()) + "Mtx file not found in %s. Plz recheck." %(self.map_result_dir))
    #        self.mobilogger._mobilogrecorder(log_message="Mtx file not found in %s. Plz recheck." %(self.map_result_dir), 
    #            log_level="ERROR")
    #        sys.exit()
    #    fetched_mtx_dir = os.path.join(self.map_result_dir, i, "raw_cell_gene_matrix")
    #    filter_mtx_dir = os.path.join(self.o_dir, "filtered_cell_gene_matrix")
    #    cmd = "STAR --outTmpDir %s --runMode soloCellFiltering %s %s --soloCellFilter %s " %(tmp_path, fetched_mtx_dir, filter_mtx_dir, self.filter_args)
    #    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    #    if os.path.exists(tmp_path):
    #        cmd = "rm -r %s " %(tmp_path)
    #        os.system(cmd)
    def process(self):
        fetched_bam_file, alignment_stat_file, fetched_mtx, filtered_mtx = re_assign_bam(map_result_dir=self.map_result_dir, 
                                                                           sample_id=self.sample_id,
                                                                           split_bam_program=self.split_bam_program,
                                                                           filter_args=self.filter_args, 
                                                                           threads=self.threads, 
                                                                           ref_json_path=self.ref_json_path, 
                                                                           read_number=self.read_number, 
                                                                           mobilogger=self.mobilogger)
        return 0
