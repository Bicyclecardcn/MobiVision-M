import subprocess
import os
import time
from configparser import ConfigParser
from glob import glob
import re 
import shutil
import sys 
import copy
sys.path.append(os.path.dirname(__file__))
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor


def CheckFastqParam(fastq_dir:str, process_dir:str, mobiexecutor:CommandExecutor, mobilogger:MobiLoggingSystem, id_reset:str):
    dirfile_list = os.listdir(fastq_dir)
    valid_type = {"1": "R1", 
                  "R1": "R1", 
                  "2": "R2",
                  "R2": "R2", 
                  "0": "R0", 
                  "R0": "R0", 
                  "3": "R3", 
                  "R3": "R3"
                } 
    filepath_list = []
    sort_files = {}
    found_type = []
    if id_reset != None and id_reset != "":
        found_name = id_reset
    else:
        found_name = None
    for fastq_file in dirfile_list:
        x_find, x_name, x_type, x_suffx, x_file_suffx = check_Fq_NamePat(fastq_file)
        if x_find:
            if (found_name == None or found_name == "") and x_name != "" and x_name != "None":
                found_name = x_name
            if not x_file_suffx.endswith(".gz"):
                shutil.copy(os.path.join(fastq_dir, fastq_file), 
                            os.path.join(process_dir, fastq_file))
                cmd = "gzip %s" %(os.path.join(process_dir, fastq_file))
                exit_code = mobiexecutor.execute(
                    command=cmd,
                    context={
                        "quantify": "gzip"
                    }, console_output=False
                )
                if exit_code != 0:
                    mobilogger._mobilogrecorder(log_message="Gzip %s failed." %(os.path.join(fastq_dir, fastq_file)),
                        log_level="ERROR")
                else:
                    fastq_file += ".gz"
                    processed_fastq_dir = process_dir
            else:
                processed_fastq_dir = fastq_dir
            filepath_list.append(os.path.join(processed_fastq_dir, fastq_file))
            non_type_name = "_".join([x_name, x_suffx, x_file_suffx])
            try:
                tmp_type = valid_type[x_type]
            except KeyError:
                mobilogger._mobilogrecorder(log_message="Invalid flag found in %s. Found flag is %s." %(fastq_file, x_type), 
                    log_level="ERROR")
                sys.exit()
            if not non_type_name in sort_files.keys():
                sort_files[non_type_name] = {tmp_type:os.path.join(processed_fastq_dir, fastq_file)}
            else:
                sort_files[non_type_name][tmp_type] = os.path.join(processed_fastq_dir, fastq_file)
            if not tmp_type in found_type:
                found_type.append(tmp_type)
    final_fastq = {}
    if len(sort_files) == 0:
        mobilogger._mobilogrecorder(log_message="No fastq files found in  %s." %(fastq_dir), 
            log_level="ERROR")
        sys.exit()
    elif len(sort_files) > 1:
        cat_dict = {}
        for i in sort_files.keys():
            apply_type = copy.copy(found_type)
            for j in sort_files[i].keys():
                if not j in cat_dict.keys():
                    cat_dict[j] = [sort_files[i][j]]
                else:
                    cat_dict[j].append(sort_files[i][j])
                apply_type.remove(j)
            if len(apply_type) > 0:
                mobilogger._mobilogrecorder(log_message="The following file type(s) not found in %s : %s " %(",".join(apply_type), i), 
                    log_level="ERROR")
                sys.exit()
        for i in cat_dict.keys():
            tmp_out_file = os.path.join(process_dir, found_name + "_%s.fastq.gz" %(i))
            cmd = "cat %s > %s " %(" ".join(cat_dict[i]), tmp_out_file)
            exit_code = mobiexecutor.execute(
                command=cmd,
                context={
                    "quantify": "cat "
                }, console_output=False
            )
            final_fastq[i] = tmp_out_file
    else:
        key, value = next(iter(sort_files.items()))
        for i in sort_files[key].keys():
            final_fastq[i] = sort_files[key][i]
    if not "R1" in final_fastq.keys():
        mobilogger._mobilogrecorder(log_message="R1 file are not found in %s" %(fastq_dir), 
            log_level="ERROR")
        sys.exit()
    elif not "R2" in final_fastq.keys():
        mobilogger._mobilogrecorder(log_message="R2 file are not found in %s" %(fastq_dir), 
            log_level="ERROR")
        sys.exit()
    return final_fastq, found_name, os.path.dirname(final_fastq["R1"])

def check_Fq_NamePat(fqname:str):
    #pattern = re.compile("(.+)[-_\\.](R1|R2|1|2)(_0\d+)*.(fastq.gz|fq.gz|fastq|fq)$")
    pattern = re.compile("(.+)[_\\.](R1|R2|1|2)[_\\.](.*)(fastq.gz|fq.gz|fastq|fq)$")
    pat_match = re.match(pattern, fqname)
    if pat_match:
        return True, pat_match.group(1), pat_match.group(2), pat_match.group(3).strip("."), pat_match.group(4)
    else:
        return False, "", "", "", ""

def str2bool(x):
    return x.lower() in ('true')

def type_overwrite(data1, data2, head):
    prog_runlog("Update %s: %s to %s." %(head, str(data1), str(data2)))
    if isinstance(data1, bool):
        return str2bool(data2)
    elif isinstance(data1, int):
        return int(data2)
    elif isinstance(data1, float):
        return float(data2)
    else:
        return str(data2)

def run_command_safely(cmd, args, cwd=None, env=None):
    """
    call the shell commands safely, args is input parameters of cmd, use the start_new_session argument
    """

    arg_join_str = ' '.join([arg if not ' ' in arg else '"%s"' % arg for arg in args] ) ## ''.join(args)

    stderr_data = ''
    assert isinstance(cmd, str) and isinstance(args, list)
    if ">" in arg_join_str:
        # check_call command:
        # example: subprocess.check_call('gzip -dc %s > %s'%(f, new_uncomp_file), shell=True)

        check_call_args = {}
        if cwd is not None:
            check_call_args['cwd'] = cwd
        if env is not None:
            check_call_args['env'] = env

        if "1>" in arg_join_str:
            logf = ""
            for i in args:
                if i.startswith("1>"):
                    logf = i.replace("1>", "")
                    args.remove(i)
                    break

            fd_stdout=open(logf,"w+")
            p = subprocess.Popen([cmd] + args, stderr=subprocess.PIPE, stdout=fd_stdout, start_new_session=True)
            _, stderr_data = p.communicate()
            fd_stdout.close()
            if p.returncode != 0:
                raise Exception("%s: %s returned error code %d: %s" % (cmd, p, p.returncode, stderr_data))
            
        else:

            args = [cmd] + args
            try:
                #print('Running command \'%s\': %s' % (name, subprocess.list2cmdline(args)))
                subprocess.check_call(args, **check_call_args)
                #print('Command \'%s\' completed successfully' % name)
            except subprocess.CalledProcessError as e:
                stderr_data = cmd + " : call error."
                raise Exception('\'%s\' exited with error code: %s' % (cmd, e.returncode))

    else:
        p = subprocess.Popen([cmd] + args, stderr=subprocess.PIPE, stdout=subprocess.PIPE, start_new_session=True)
        _, stderr_data = p.communicate()
        if p.returncode != 0:
            raise Exception("%s: %s returned error code %d: %s" % (cmd, p, p.returncode, stderr_data))
    
    return stderr_data
    
## check all supplymentary data for sequencing

def check_fq_supp(in_dir):

    """  check the fastq files in input dir , if multiple R1 and R2, then merge those files in a merged directory . """
    
    all_f_names = glob(os.path.join(in_dir, "*"))
    
    ## data names store all fastq or fastq.gz file names in inputdir
    data_names = []
    
    for f in all_f_names:
        base_name = os.path.basename(f)
        ## check filename contain fastq or fq , if yes then store this filename
        if os.path.isfile(f) and (base_name.endswith('fq.gz') or base_name.endswith('fastq.gz') or base_name.endswith('fq') or base_name.endswith('fastq')):
            data_names.append(base_name)
    if len(data_names) >=2:
        if "_" in data_names[0]:
            tmp_name = data_names[0].split("_")[0]
        else:
            tmp_name = data_names[0]
        if len(tmp_name) > 20:
            tmp_name = tmp_name[:20]
        if len(data_names) <=2:
            return os.path.abspath(in_dir), tmp_name, True
        ## sort filename in data_names
        sort_names = sorted(data_names)
        
        merged_names = sort_names[:2]
        R1_files = [i for n,i in enumerate(sort_names) if n%2==0]
        R2_files = [i for n,i in enumerate(sort_names) if n%2==1]
        
        mergedir = os.path.join(os.path.abspath(in_dir), "merged")
        if not os.path.exists(mergedir):
            os.makedirs( mergedir,  mode=0o755 )
        
        mergeR1 = os.path.join(mergedir, merged_names[0])
        mergeR2 = os.path.join(mergedir, merged_names[1])
        
        cmd_R1 = " ".join([os.path.join(in_dir, f) for f in R1_files])
        cmd_R2 = " ".join([os.path.join(in_dir, f) for f in R2_files])
        if merged_names[0].endswith("gz"):
            subprocess.check_call("cat %s > %s"%(cmd_R1, mergeR1), shell=True)
            subprocess.check_call("cat %s > %s"%(cmd_R2, mergeR2), shell=True)
        else:
            cmd_R1 = "cat "  + cmd_R1
            cmd_R2 = "cat "  + cmd_R2
            mergeR1 = mergeR1+".gz"
            mergeR2 = mergeR2+".gz"
            subprocess.check_call("%s > %s"%(cmd_R1, mergeR1), shell=True)
            subprocess.check_call("%s > %s"%(cmd_R2, mergeR2), shell=True)
            subprocess.check_call("gzip %s"%(mergeR1), shell=True)
            subprocess.check_call("gzip %s"%(mergeR2), shell=True)
        found = True
    else:
        mergedir = in_dir
        tmp_name = "Less than 2 fastq file found in %s. The name if fastq file should endwith 'fastq', 'fq'. 'fastq.gz' or 'fq.gz'."
        found = False
    return mergedir, tmp_name, found
    
def prog_runlog(*args, **kwargs):
    # *kargs 为了通用, 可不传
    with open("./_log", "a") as f: # 打开文件 把print函数输出的数据写入到文件
        print(*args, file=f, **kwargs)

def readin_config(config_file:str, with_CB = False, UMI_adjust = "step_1"):
    found_config = False
    if os.path.exists(config_file):
        found_config = True
    else:
        config_file = os.path.split(os.path.realpath(__file__))[0] + "STAR.ini"
        if os.path.exists(config_file):
            found_config = True
    if with_CB:
        filter_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'V2.tsb')
    else:
        filter_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'V2_without_CB.tsb')
    default_config = {"Mobivision-M":{}, "STAR":{}}
    default_config["STAR"]["nExpectedCells"] = '3000'
    default_config["STAR"]["maxPercentile"] = '0.99'
    default_config["STAR"]["maxMinRatio"] = '10'
    default_config["STAR"]["indMin"] = '45000'
    default_config["STAR"]["indMax"] = '90000'
    default_config["STAR"]["umiMin"] = '150'
    default_config["STAR"]["umiMinFracMedian"] = '0.01'
    default_config["STAR"]["candMaxN"] = '20000'
    default_config["STAR"]["FDR"] = '0.01'
    default_config["STAR"]["simN"] = '10000'
    default_config["STAR"]["soloUMIfiltering"] = "MultiGeneUMI_CR"
    default_config["STAR"]["soloUMIdedup"] = "1MM_CR"
    default_config["STAR"]["allow_multi_target_UMI"] = str2bool("False")
    default_config["STAR"]["reclaim_UMI"] = str2bool("True")
    default_config["STAR"]["white_list_file"] = "NA"
    default_config["Mobivision-M"]["filter_pattern_file"] = filter_file_path
    default_config["Mobivision-M"]["output_pattern_file"] = "NA"
    default_config["Mobivision-M"]["output_type"] = "BacDrop"
    default_config["Mobivision-M"]["encrypt"] = "DES"
    default_config["Mobivision-M"]["encrypt_key"] = "mobidrop"
    default_config["Mobivision-M"]["seq_barcode_tag"] = "NA"
    default_config["Mobivision-M"]["barcode_len"] = 20
    default_config["Mobivision-M"]["seq_UMI_tag"] = "UMI"
    default_config["Mobivision-M"]["UMI_start"] = 21
    default_config["Mobivision-M"]["UMI_len"] = 10
    default_config["Mobivision-M"]["split_by"] = "NA"
    default_config["Mobivision-M"]["split_func"] = "GO"
    default_config["Mobivision-M"]["go_script"] = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))+"/sclib/scrna/main_pre-process"
    default_config["Mobivision-M"]["adpator_list_path"] = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'adaptor_list_BacDrop_v0.1.txt')
    default_config["Mobivision-M"]["fastp_adapter_path"] = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'fastp_adapter.fasta')
    UMI_method_dict = {"no_adjust":{"soloUMIfiltering":"MultiGeneUMI_CR", 
                                    "soloUMIdedup":"1MM_CR", 
                                    "allow_multi_target_UMI":False, 
                                    "reclaim_UMI":True}, 
                        "step_1":{"soloUMIfiltering":"-", 
                                    "soloUMIdedup":"1MM_CR", 
                                    "allow_multi_target_UMI":True, 
                                    "reclaim_UMI":False}, 
                        "step_1_and_2":{"soloUMIfiltering":"MultiGeneUMI_CR", 
                                    "soloUMIdedup":"1MM_CR", 
                                    "allow_multi_target_UMI":False, 
                                    "reclaim_UMI":False}}
    tmp_UMI_method = UMI_adjust
    if not found_config:
        pass
    else:
        prog_runlog(time.strftime("%Y-%m-%d %H:%M:%S\t", time.localtime())+"Updating config...")
        cfg = ConfigParser()
        cfg.optionxform = str
        cfg.read(config_file)
        star_config = dict(cfg.items("STAR"))
        for x in star_config.keys():
            if x == "UMI_adjust_method":
                tmp_UMI_method = star_config[x]
            default_config["STAR"][x] = type_overwrite(default_config["STAR"][x], star_config[x], x)
        run_config = dict(cfg.items("Mobivision-M"))
        for x in run_config.keys():
            default_config["Mobivision-M"][x] = type_overwrite(default_config["Mobivision-M"][x], run_config[x], x)
    if tmp_UMI_method in UMI_method_dict.keys():
        for x in UMI_method_dict[tmp_UMI_method].keys():
            prog_runlog("Update %s: %s to %s." %(x, str(default_config["STAR"][x]), str(UMI_method_dict[tmp_UMI_method][x])))
            default_config["STAR"][x] = UMI_method_dict[tmp_UMI_method][x]
    ###relative to abs
    relative_path = os.path.join(os.path.abspath(os.path.dirname(__file__)))
    check_key = ["filter_pattern_file", "output_pattern_file"]
    for tmp_key in check_key:
        if not os.path.exists(default_config["Mobivision-M"][tmp_key]):
            try_path = os.path.join(relative_path, default_config["Mobivision-M"][tmp_key])
            if os.path.exists(try_path):
                default_config["Mobivision-M"][tmp_key] = try_path
    return default_config