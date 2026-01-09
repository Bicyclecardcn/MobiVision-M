#! -*- utf-8 -*-

import os
import sys
import shutil
import pandas as pd
import datetime
import json
import traceback
import time
sys.path.append(os.path.dirname(__file__))
from dataexp import Data_InjectTool
from mofhandle import file_rearrange
import configparser
from configparser import ConfigParser
from mako_report import ExportReport
resource_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources", "report")
from cmdutil import CheckFastqParam
from BacDrop_preprocess import BacPreProcess
from re_assign_multi import re_assign_bam
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor


def str2bool(x):
    return x.lower() in ('true')

def update_config(star_config: str, with_CB: bool, species_number: int, UMI_adjust:str, mobilogger:MobiLoggingSystem):
    ###default config
    config_file = star_config
    found_config = False
    if os.path.exists(config_file):
        found_config = True
    else:
        config_file = os.path.split(os.path.realpath(__file__))[0] + "STAR.ini"
        if os.path.exists(config_file):
            found_config = True
    if with_CB:
        filter_file_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources", "barcodes", "V2.tsb")
    else:
        filter_file_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources", "barcodes", "V2_without_CB.tsb")
    default_config = {"Mobivision-M":{}, "STAR":{}}
    default_config["STAR"]["soloMultiMappers"] = "Unique"
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
    default_config["STAR"]["soloUMIfiltering"] = "-"
    default_config["STAR"]["soloUMIdedup"] = "Exact"
    default_config["STAR"]["allow_multi_target_UMI"] = str2bool("True")
    default_config["STAR"]["reclaim_UMI"] = str2bool("True")
    default_config["STAR"]["white_list_file"] = "NA"
    default_config["STAR"]["outFilterMismatchNmax"] = 10
    default_config["STAR"]["outFilterMismatchNoverLmax"] = 0.3
    default_config["STAR"]["outFilterMismatchNoverReadLmax"] = 1
    default_config["STAR"]["outFilterMatchNmin"] = 0
    default_config["STAR"]["outFilterMatchNminOverLread"] = 0.66
    default_config["STAR"]["outFilterScoreMinOverLread"] = 0.66
    default_config["STAR"]["outFilterScoreMin"] = 30
    default_config["STAR"]["outFilterMultimapNmax"] = 20
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
    default_config["Mobivision-M"]["go_script"] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin", "main_pre-process")
    default_config["Mobivision-M"]["process_cutadapt"] = True
    default_config["Mobivision-M"]["process_fastp"] = True
    default_config["Mobivision-M"]["adpator_list_path"] = None
    default_config["Mobivision-M"]["fastp_adapter_path"] = None
    default_config["Mobivision-M"]["logo_path"] = resource_dir + "/image/logo.png"
    default_config["Mobivision-M"]["multi_mappers"] = "Unique"
    UMI_method_dict = {"no_adjust":{"soloUMIfiltering":"-", 
                                    "soloUMIdedup":"Exact", 
                                    "allow_multi_target_UMI":True, 
                                    "reclaim_UMI":True}, 
                        "step_1":{"soloUMIfiltering":"-", 
                                    "soloUMIdedup":"1MM_CR", 
                                    "allow_multi_target_UMI":True, 
                                    "reclaim_UMI":False}, 
                        "step_1_and_2":{"soloUMIfiltering":"-", 
                                    "soloUMIdedup":"1MM_CR", 
                                    "allow_multi_target_UMI":False, 
                                    "reclaim_UMI":False}}
    tmp_UMI_method = UMI_adjust
    if not found_config:
        pass
    else:
        mobilogger._mobilogrecorder(log_message="Updating run config by arguments.",
            log_level="INFO")
        cfg = ConfigParser()
        cfg.optionxform = str
        cfg.read(config_file)
        update_star = True
        try:
            stararg_config = dict(cfg.items("STAR"))
        except configparser.NoSectionError:
            mobilogger._mobilogrecorder(log_message="Section 'STAR' not found in config file %s." %(config_file),
                log_level="WARNING")
            update_star = False
        except TypeError:
            stararg_config = cfg.items("STAR")
        if update_star:
            for x in stararg_config.keys():
                if x == "UMI_adjust_method":
                    tmp_UMI_method = stararg_config[x]
                if x == "multi_mappers":
                    if not stararg_config[x] in ["Unique", "Uniform", "Rescue", "PropUnique", "EM"]:
                        mobilogger._mobilogrecorder(log_message="The value %s of %s in STAR section is not a valid argument." \
                        %(stararg_config[x], x), log_level="WARNING")
                    else:
                        default_config["STAR"][x] = type_overwrite(default_config["STAR"][x], stararg_config[x], x, mobilogger)
                elif not x in default_config["STAR"].keys():
                    mobilogger._mobilogrecorder(log_message="The %s in STAR section is not a valid argument." %(x),
                        log_level="WARNING")
                else:
                    default_config["STAR"][x] = type_overwrite(default_config["STAR"][x], stararg_config[x], x, mobilogger)
        update_run = True
        try:
            run_config = dict(cfg.items("MobiVision-M"))
        except configparser.NoSectionError:
            mobilogger._mobilogrecorder(log_message="Section 'MobiVision-M' not found in config file %s." %(config_file),
                log_level="WARNING")
            update_run = False
        except TypeError:
            run_config = cfg.items("MobiVision-M")
        if update_run:
            for x in run_config.keys():
                if not x in default_config["Mobivision-M"].keys():
                    mobilogger._mobilogrecorder(log_message="The %s in Mobivision-M section is not a valid argument." %(x),
                        log_level="WARNING")
                else:
                    default_config["Mobivision-M"][x] = type_overwrite(default_config["Mobivision-M"][x], run_config[x], x, mobilogger)
    if not tmp_UMI_method == "no_adjust":
        if tmp_UMI_method in UMI_method_dict.keys():
            for x in UMI_method_dict[tmp_UMI_method].keys():
                mobilogger._mobilogrecorder(log_message="Update %s: %s to %s." %(x, str(default_config["STAR"][x]), str(UMI_method_dict[tmp_UMI_method][x])),
                    log_level="INFO")
                default_config["STAR"][x] = UMI_method_dict[tmp_UMI_method][x]
    mobilogger._mobilogrecorder(log_message="Config Updating is done successfully.",
        log_level="INFO")
    return default_config

def process_filter_arguments(top_cell: int, cr_flag: bool, hard_filter: str, config: dict):
    return_info = ""
    return_flag ="INFO"
    cell_filter_val = ""
    if top_cell != None:
        if not cr_flag and hard_filter == None:
            cell_filter_val = "TopCells " + str(top_cell)
        else:
            return_flag = "ERROR"
            return_info = "Detect conflicts. --cellnumber assigned and one of --cr2.2 or --hard_filter assigened. %s %s " %(cr_flag, hard_filter)
    if hard_filter != None:
        if not cr_flag and top_cell == None:
            tmp_split = hard_filter.split(":")
            if tmp_split[0] in ["min_UMI", "min_reads"]:
                try:
                    tmp_value = float(tmp_split[1])
                except Exception as e:
                    return_flag = "ERROR"
                    return_info = "The value of --hard_filter end with an positive integer. The input value was %s" %(hard_filter)
                else:
                    if tmp_value < 0 :
                        return_flag = "ERROR"
                        return_info = "The value of --hard_filter must end with an positive integer. The input value was %s" %(hard_filter)
                    else:
                        cell_filter_val = hard_filter
            else:
                return_flag = "ERROR"
                return_info = "The value of --hard_filter must start with 'min_UMI' or 'min_reads' and end with an positive integer. The input value was %s" %(hard_filter)
        else:
            return_flag = "ERROR"
            return_info = "Detect conflicts. --hard_filter assigned and one of --cellnumber or --cr2.2 assigned. %s %s " %(cr_flag, top_cell)
    if not cr_flag and cell_filter_val == "":
        if top_cell == None and hard_filter == None:
            cell_filter_val = "EmptyDrops_CR %s %s %s %s %s %s %s %s %s %s" %(config["STAR"]["nExpectedCells"], 
                                                                            config["STAR"]["maxPercentile"], 
                                                                            config["STAR"]["maxMinRatio"], 
                                                                            config["STAR"]["indMin"], 
                                                                            config["STAR"]["indMax"], 
                                                                            config["STAR"]["umiMin"], 
                                                                            config["STAR"]["umiMinFracMedian"], 
                                                                            config["STAR"]["candMaxN"], 
                                                                            config["STAR"]["FDR"], 
                                                                            config["STAR"]["simN"])
    if cr_flag and cell_filter_val == "":
        if top_cell == None and hard_filter == None:
            cell_filter_val = "CellRanger2.2 %s %s %s " %(config["STAR"]["nExpectedCells"], 
                                                        config["STAR"]["maxPercentile"], 
                                                        config["STAR"]["maxMinRatio"],) 
        else:
            return_flag = "ERROR"
            return_info = "Detect conflicts. --cr2.2 assigned and one of --cellnumber or --hard_filter assigened. %s %s " %(top_cell, hard_filter)
    if return_info == "":
        return_info = "Using '%s' as microbe filter algorithm." %(cell_filter_val)
    return cell_filter_val, return_flag, return_info

def type_overwrite(data1, data2, head, mobilogger):
    if data2 == "None":
        return_value = None
    elif isinstance(data1, bool):
        try:
            return_value = str2bool(data2)
        except Exception as e:
            mobilogger._mobilogrecorder(log_message="Update %s failed. Can't Transfer %s to Bool type. Error %s happened." %(head, str(data2), e),
                log_level="ERROR")
            sys.exit()
    elif isinstance(data1, int):
        try:
            return_value = int(data2)
        except Exception as e:
            mobilogger._mobilogrecorder(log_message="Update %s failed. Can't Transfer %s to Int type. Error %s happened." %(head, str(data2), e),
                log_level="ERROR")
            sys.exit()
    elif isinstance(data1, float):
        try:
            return_value = float(data2)
        except Exception as e:
            mobilogger._mobilogrecorder(log_message="Update %s failed. Can't Transfer %s to Float type. Error %s happened." %(head, str(data2), e),
                log_level="ERROR")
            sys.exit()
        return float(data2)
    else:
        return_value = str(data2)
    mobilogger._mobilogrecorder(log_message="Update %s: %s to %s." %(head, str(data1), str(data2)),
        log_level="INFO")
    return return_value
            
class MobivisionProcessM:
    def __init__(self, fastq_path: str, core_num: int, ref_path: str, out_path: str, vers_n: str, intron_incOpt: str, 
                 topcell_n:int, emptydrop_opt:bool, harder_filter: str, vers_kit: str, \
                 run_cmd: str, mobilogger: MobiLoggingSystem, new_sampleid=None, with_CB=False, star_config="NA", dev_mod=False, UMI_adjust="no_adjust", 
                 nosecondary=False, keep_bam=False, keep_unmapped=False, multiplet_method="auto", qc_only=False, 
                 host_remove=False, host_reference="NA"):
        self.mobilogger = mobilogger
        self.mobicommandlogger = MobiCommandLogSystem(o_dir=self.mobilogger.working_path, dev_mode=False)
        self.mobiexecutor = CommandExecutor(log_system=self.mobicommandlogger, console_output=False)
        self.fastq_path = fastq_path
        if core_num< 8:
            self.core_num = 8
        elif core_num > 64:
            self.core_num = 64
        else:
            self.core_num = core_num
        self.out_path = {}
        self.out_path["ALL"] = out_path
        self.ref_path = ref_path.rstrip('/')
        with open(os.path.join(ref_path, 'reference.json'), 'r', encoding='utf-8') as file:
            ref_dcit = json.load(file)
        self.species_num = len(ref_dcit["genomes"])
        self.id_reset = new_sampleid
        self.summary_data_dir = {}
        self.summary_data_dir["ALL"] = os.path.join(self.out_path["ALL"], 'summary_data')
        ###pre process dir 
        self.pre_process_dir = os.path.join(self.out_path["ALL"], 'pre_process')
        if not os.path.exists(self.pre_process_dir):
            os.makedirs(self.pre_process_dir)
        ###CB
        self.with_CB = with_CB
        self.map_result_dir = {}
        self.map_result_dir["ALL"] = os.path.join(self.out_path["ALL"], 'map_result')
        self.intron_iopt = intron_incOpt
        self.topcell_num = topcell_n
        self.empty_paras = emptydrop_opt
        self.hard_filter = harder_filter
        self.vers_kit = vers_kit
        self.config = star_config
        d_attrs = {"excluded":"Gene", "included":"GeneFull"}
        self.intron_option = d_attrs.get(self.intron_iopt, "Gene")
        self.run_cmd = run_cmd 
        self.rawfq_path, self.sample_id, self.fastq_path = CheckFastqParam(fastq_dir=fastq_path, 
                                                                              process_dir=self.pre_process_dir, 
                                                                              mobilogger=self.mobilogger, 
                                                                              mobiexecutor=self.mobiexecutor, 
                                                                              id_reset=self.id_reset)
        mobilogger._mobilogrecorder(log_message="Using %s as ID." %(self.sample_id), log_level="INFO")
        self.out_sample_id = self.sample_id
        self.id_reset = self.sample_id
        self.R1_file = self.rawfq_path["R1"].split("/")[-1]
        self.R2_file = self.rawfq_path["R2"].split("/")[-1]
        self.dev_mod = dev_mod
        self.split_bam_program = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin", "fetch_multihit")
        self.UMI_adjust = UMI_adjust
        self.nosecondary = nosecondary
        self.keep_bam = keep_bam
        self.keep_unmapped = keep_unmapped
        self.multiplet_method = multiplet_method
        self.qc_only = qc_only
        self.host_remove = host_remove
        self.host_reference = host_reference

    def map_result(self):
        if not os.path.isdir(self.map_result_dir["ALL"]):
            os.makedirs(self.map_result_dir["ALL"])
        r1_path = None
        r2_path = None
        white_list_path = None
        for root, dirs, files in os.walk(self.fastq_path):
            for fs in files:
                f = os.path.join(root, fs)
                if '_S0_L001_R1_001' in f:
                    r1_path = f
                elif '_S0_L001_R2_001' in f:
                    r2_path = f
                elif "white_list.tsv" in f:
                    white_list_path = f
        if r1_path == None:
            self.mobilogger._mobilogrecorder(log_message="R1 file not found in %s." %(self.fastq_path),
                log_level="ERROR")
            sys.exit()
        else:
            self.mobilogger._mobilogrecorder(log_message="Using %s as R1 file." %(r1_path),
                log_level="INFO")
        if r2_path == None:
            self.mobilogger._mobilogrecorder(log_message="R2 file not found in %s." %(self.fastq_path),
                log_level="ERROR")
            sys.exit()
        else:
            self.mobilogger._mobilogrecorder(log_message="Using %s as R2 file." %(r2_path),
                log_level="INFO")
        if white_list_path == None:
            self.mobilogger._mobilogrecorder(log_message="white list file not found in %s." %(self.fastq_path),
                log_level="ERROR")
            sys.exit()
        else:
            self.mobilogger._mobilogrecorder(log_message="Using %s as white list file." %(white_list_path),
                log_level="INFO")
        map_result_path = os.path.join(self.map_result_dir["ALL"], self.sample_id + '_')
        strand_option = "Unstranded"
        if self.__run_config["Mobivision-M"]["output_type"] == "Bacdrop":
            if self.with_CB:
                CBlen = 30
                UMIstart = 31
                UMIlen = 8
            else:
                CBlen = 20
                UMIstart = 21
                UMIlen = 10 
        else:
            CBlen = self.__run_config["Mobivision-M"]["barcode_len"]
            UMIstart = self.__run_config["Mobivision-M"]["UMI_start"]
            UMIlen = self.__run_config["Mobivision-M"]["UMI_len"]
        if self.keep_unmapped:
            out_unmapped = "Fastx"
        else:
            out_unmapped = "None"
        tmp_path = os.path.dirname(map_result_path)
        if not os.path.exists(tmp_path):
            os.makedirs(tmp_path)
        tmp_path = os.path.join(tmp_path, "_STARtmp")
        if self.__run_config["STAR"]["soloMultiMappers"] == "Unique":
            run_mapper = "Uniform"
        else:
            run_mapper = self.__run_config["STAR"]["soloMultiMappers"]
        cmd = 'STAR --runThreadN {} --genomeDir {} --outFilterMultimapNmax {} --quantMode GeneCounts ' \
              '--soloType Droplet --soloCBwhitelist {} --soloCBlen {} --soloUMIstart {} --soloUMIlen {} ' \
              '--soloStrand {} --soloCBmatchWLtype Exact --readFilesCommand zcat ' \
              '--soloUMIfiltering {} --soloUMIdedup {} --outFilterScoreMin {} ' \
              '--soloCellFilter {} --soloFeatures {} --soloMultiMappers {} --soloBarcodeReadLength 0 --outSAMattributes GX GN CB ' \
              'UB NH HI --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 65719476736 --readFilesIn {} {} ' \
              '--outFileNamePrefix {} --outTmpDir {} --outReadsUnmapped {} --outBAMsortingThreadN {}' \
              '--outFilterMismatchNmax {} --outFilterMismatchNoverLmax {} --outFilterMismatchNoverReadLmax {}' \
              '--outFilterMatchNmin {} --outFilterMatchNminOverLread {} --outFilterScoreMinOverLread {}'.format(self.core_num , self.ref_path+"/star" , 
                                              self.__run_config["STAR"]["outFilterMultimapNmax"], white_list_path , 
                                              str(CBlen), str(UMIstart), str(UMIlen), strand_option ,
                                              self.__run_config["STAR"]['soloUMIfiltering'], self.__run_config["STAR"]['soloUMIdedup'], 
                                              self.__run_config["STAR"]["outFilterScoreMin"], 
                                              "None", self.intron_option, run_mapper, 
                                              r2_path, r1_path, map_result_path, tmp_path, out_unmapped, self.core_num, 
                                              self.__run_config["STAR"]["outFilterMismatchNmax"], self.__run_config["STAR"]["outFilterMismatchNoverLmax"], 
                                              self.__run_config["STAR"]["outFilterMismatchNoverReadLmax"], self.__run_config["STAR"]["outFilterMatchNmin"], 
                                              self.__run_config["STAR"]["outFilterMatchNminOverLread"], self.__run_config["STAR"]["outFilterScoreMinOverLread"])
        self.mobilogger._mobilogrecorder(log_message="Mapping to reference...",
            log_level="INFO")
        if self.dev_mod:
            exit_code = self.mobiexecutor.execute(
                command=cmd,
                context={
                    "quantify": "STAR"
                }, console_output=True
            )
        else:
            exit_code = self.mobiexecutor.execute(
                command=cmd,
                context={
                    "quantify": "STAR"
                }, console_output=False
            )
        if exit_code != 0:
            self.mobilogger._mobilogrecorder(log_message="Mapping failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                log_level="ERROR")
        else:
            self.mobilogger._mobilogrecorder(log_message="Mapping is done successfully.",
                log_level="INFO")
        raw_mtx_dir = map_result_path + "Solo.out/" + self.intron_option + "/raw"
        summary_dir = map_result_path + "Solo.out/" + self.intron_option + "/Summary.csv"
        return True, raw_mtx_dir, summary_dir, self.__run_config["STAR"]["allow_multi_target_UMI"], self.__run_config["STAR"]["reclaim_UMI"]
    
    def export_json(self, final_map_stats_flie, raw_mtx_dir, final_mtx):
        pre_file = os.path.join(self.pre_process_dir, self.sample_id + "_filter_stat.tsv")
        Exp_test = Data_InjectTool(sample_ID=self.sample_id, 
                            output_path=os.path.join(self.out_path["ALL"], self.sample_id + "_outs"),
                            pre_stat_file=pre_file,
                            input_stat_file=os.path.join(self.pre_process_dir, self.sample_id + "_qc_stat.tsv"),
                            filter_stats_file=os.path.join(self.pre_process_dir, self.sample_id + "_filter_stat.tsv"),
                            saturation_file=final_map_stats_flie, 
                            raw_mtx_path=raw_mtx_dir,
                            filter_mtx_path=final_mtx,
                            cell_stat_file=os.path.join(self.map_result_dir["ALL"], "cell_stat.tsv"), 
                            reference_json_path=os.path.join(self.ref_path, "reference.json"), 
                            kit=self.vers_kit, 
                            run_cmd=self.run_cmd, 
                            threads=self.core_num,
                            run_cluster=(not self.nosecondary),
                            multiplet_method=self.multiplet_method, 
                            host_remove=self.host_remove, 
                            mobilogger=self.mobilogger, 
                            dev_mod=self.dev_mod)
        report_json = Exp_test.export_json_microbe()
        return report_json
    
    def process(self):
        self.mobilogger._mobilogrecorder(log_message="Mobivision-M Analysis Start.", 
            log_level="INFO")
        sub_process_data = pd.DataFrame(columns=["start", "end", "process"])
        ###note start pre-process time
        now = datetime.datetime.now()
        formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
        current_add = sub_process_data.shape[0]
        sub_process_data.loc[current_add, "start"] = formatted_date
        sub_process_data.loc[current_add, "process"] = "pre-process" 
        ###prepare config
        default_config = update_config(star_config=self.config,
                                       with_CB=self.with_CB, 
                                       species_number=self.species_num,  
                                       UMI_adjust=self.UMI_adjust, 
                                       mobilogger=self.mobilogger)
        ###check and make filter argument
        self.cell_filter_val, check_flag, check_info = process_filter_arguments(top_cell=self.topcell_num,
                                                                           cr_flag=self.empty_paras, 
                                                                           hard_filter=self.hard_filter, 
                                                                           config=default_config)
        self.mobilogger._mobilogrecorder(log_message=check_info, 
            log_level=check_flag)
        if check_flag == "ERROR":
            sys.exit()
        ###relative to abs
        relative_path = os.path.join(os.path.abspath(os.path.dirname(__file__)))
        check_key = ["filter_pattern_file", "output_pattern_file"]
        for tmp_key in check_key:
            if not os.path.exists(default_config["Mobivision-M"][tmp_key]):
                try_path = os.path.join(relative_path, default_config["Mobivision-M"][tmp_key])
                if os.path.exists(try_path):
                    default_config["Mobivision-M"][tmp_key] = try_path
        self.__run_config = default_config
        if not os.path.exists(os.path.join(self.pre_process_dir, "Job_done.flag")):
            mp = BacPreProcess(
                fastq_path=self.fastq_path,
                R1_file=self.R1_file,
                R2_file=self.R2_file, 
                output_path=self.pre_process_dir, 
                threads=self.core_num, 
                sample_name=self.sample_id, 
                pre_process="TRUE", 
                process_cutadapt=default_config["Mobivision-M"]["process_cutadapt"], 
                process_fastp=default_config["Mobivision-M"]["process_fastp"],
                adaptor_list=default_config["Mobivision-M"]["adpator_list_path"], 
                host_remove=self.host_remove, 
                host_ref=self.host_reference, 
                download_path=os.path.join(self.pre_process_dir + "/download_path"), 
                fastp_adapter=default_config["Mobivision-M"]["fastp_adapter_path"], 
                qualified_quality_phred=15, 
                unqualified_percent_limit=40, 
                n_base_limit=5, 
                length_required=30, 
                lib_type="Illumina", 
                with_CB = self.with_CB, 
                star_config = self.__run_config,
                dev_mod=self.dev_mod, 
                mobilogger=self.mobilogger)
            mp.process()
        if not self.qc_only:
            ori_sample_id = self.sample_id
            ori_outpath = self.out_path["ALL"]
            ready_sample = pd.read_csv(os.path.join(self.pre_process_dir, self.id_reset + "_sample_barcode_stat.tsv"), sep="\t")
            ###note the end pre-process time
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            sub_process_data.loc[current_add, "end"] = formatted_date
            ###note start map time
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            current_add = sub_process_data.shape[0]
            sub_process_data.loc[current_add, "start"] = formatted_date
            sub_process_data.loc[current_add, "process"] = "mapping" 
            temp_CB = ready_sample.loc[0, "sample_barcode"]
            self.fastq_path = os.path.join(self.pre_process_dir, "clean_data", temp_CB)
            self.sample_id = ori_sample_id
            self.out_path["ALL"] = os.path.join(ori_outpath, temp_CB)
            self.map_result_dir["ALL"] = os.path.join(self.out_path["ALL"], 'map_result')
            self.summary_data_dir["ALL"] = os.path.join(self.out_path["ALL"], 'summary_data')
            if not os.path.exists(self.map_result_dir["ALL"]):
                os.makedirs(self.map_result_dir["ALL"])
            if not os.path.exists(self.out_path["ALL"]):
                os.makedirs(self.out_path["ALL"])
            if not os.path.isdir(self.summary_data_dir["ALL"]):
                os.makedirs(self.summary_data_dir["ALL"])
            temp_name = self.id_reset + "_filter_stat.tsv"
            temp_filter_stat = pd.read_csv(os.path.join(self.pre_process_dir, temp_name), sep="\t")
            estimate_read_number = temp_filter_stat.loc[0, "passed_reads"]
            ###mapping
            map_status, raw_mtx_dir, summary_dir, allow_multi_target_UMI, reclaim_UMI = self.map_result()
            ###note end mapping time
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            sub_process_data.loc[current_add, "end"] = formatted_date
            ###note start fetch multi time
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            current_add = sub_process_data.shape[0]
            sub_process_data.loc[current_add, "start"] = formatted_date
            sub_process_data.loc[current_add, "process"] = "re-assign_mulit-alignments" 
            ###fetch mulihit
            fetched_bam_file, alignment_stat_file, fetched_mtx, filter_fetched_mtx = re_assign_bam(map_result_dir=self.map_result_dir["ALL"], 
                                                                    raw_mtx_dir=raw_mtx_dir,
                                                                    sample_id=self.sample_id,
                                                                    split_bam_program=self.split_bam_program,
                                                                    threads=self.core_num, 
                                                                    filter_args=self.cell_filter_val, 
                                                                    ref_json_path=os.path.join(self.ref_path, "reference.json"), 
                                                                    read_number=estimate_read_number, 
                                                                    dev_mod=self.dev_mod, 
                                                                    allow_multi_target_UMI=allow_multi_target_UMI, 
                                                                    reclaim_UMI=reclaim_UMI, 
                                                                    mobilogger=self.mobilogger, 
                                                                    method = default_config["STAR"]["soloMultiMappers"])
            ###re call cell 
            final_mtx = filter_fetched_mtx
            final_map_stats_flie = alignment_stat_file
            ###note end fetch multi time
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            sub_process_data.loc[current_add, "end"] = formatted_date
            ###note start report generating time
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            current_add = sub_process_data.shape[0]
            sub_process_data.loc[current_add, "start"] = formatted_date
            sub_process_data.loc[current_add, "process"] = "generating_report" 
            self.mobilogger._mobilogrecorder(log_message="Summarizing analysis results...", 
                log_level="INFO")
            ###export json
            try_times = 0
            while try_times < 2:
                #p = multiprocessing.Process(target=self.export_json, args=(final_map_stats_flie, 
                #                                                            raw_mtx_dir, 
                                                                            #fetched_mtx,
                #                                                            final_mtx))
                try_times += 1
                #p.start()
                #p.join()
                try:
                    tmp_json = self.export_json(final_map_stats_flie, raw_mtx_dir, final_mtx)
                except Exception as e:
                    self.mobilogger._mobilogrecorder(log_message="Summarizing analysis results failed. May because of low menory, retry in 1 minute.", 
                        log_level="WARNING")
                    traceback.print_exc()
                    time.sleep(60)
                else:
                    self.mobilogger._mobilogrecorder(log_message="Summarizing analysis results is done successfully.", 
                        log_level="INFO")
                    break
            report_json = os.path.join(self.out_path["ALL"], self.sample_id + "_outs", "report.json")
            if not os.path.exists(report_json):
                self.mobilogger._mobilogrecorder(log_message="Summarizing analysis results failed.", 
                        log_level="ERROR")
                sys.exit()
            ###create report
            self.mobilogger._mobilogrecorder(log_message="Making report...", 
                log_level="INFO")
            report_out_path = os.path.join(self.out_path["ALL"], self.sample_id + "_outs")
            p = ExportReport(template_file=resource_dir + "/report_template.html", 
                        json_file=report_json, 
                        output_file=os.path.join(report_out_path, self.sample_id + "_report.html"), 
                        jquery=os.path.join(resource_dir, "js", "jquery-latest.min.js"), 
                        plotly=os.path.join(resource_dir, "js", "plotly-latest.min.js"),
                        favicon_file=os.path.join(resource_dir, "image", "favicon.ico"), 
                        web_logo=default_config["Mobivision-M"]["logo_path"], 
                        web_back=os.path.join(resource_dir, "image","back.png"))
            p.process()
            tmp_files = os.listdir(report_out_path)
            if self.dev_mod:
                dev_path = os.path.join(self.out_path["ALL"], "dev_files")
                if not os.path.exists(dev_path):
                    os.makedirs(dev_path)
            for i in tmp_files:
                if i.endswith("stats.csv") or i.endswith('stat.tsv') or i.endswith('genes.tsv') or i.endswith('.out') or i.endswith('.log'):
                    if self.dev_mod:
                        shutil.move(os.path.join(report_out_path, i), dev_path)
                    else:
                        os.remove(os.path.join(report_out_path, i))
            self.mobilogger._mobilogrecorder(log_message="Making report is done successfully.", 
                log_level="INFO")
            #re-arrange output
            final_check = file_rearrange(raw_fdir=self.out_path["ALL"], 
                                        filter_mtx=final_mtx, 
                                        raw_mtx=raw_mtx_dir, 
                                        fetched_mtx=fetched_mtx, 
                                        summary_dir=summary_dir, 
                                        results_id=self.sample_id, 
                                        geneFeature_Output=self.intron_option, 
                                        output_path=report_out_path, 
                                        keep_bam=self.keep_bam, 
                                        dev_mod=self.dev_mod, 
                                        keep_qc=self.qc_only)
            final_check.fileMove()
        ## 分析正常结束， 输出结束时间
        self.mobilogger._mobilogrecorder(log_message="Whole process of MobiVision-M analysis is done successfully", 
            log_level="INFO")
        ###note end generating report time
        now = datetime.datetime.now()
        formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
        sub_process_data.loc[current_add, "end"] = formatted_date
        if self.dev_mod:
            temp_name = os.path.join(self.out_path["ALL"], "sub_process_annnote.tsv")
            sub_process_data.to_csv(temp_name, sep="\t", index=False)
