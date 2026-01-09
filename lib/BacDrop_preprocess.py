# -*- coding: utf-8 -*-
"""
Created on Wed May 18 09:04:32 2022

@author: BCc
"""

import gzip 
import pandas as pd
import os
import argparse
import numpy as np
import json
import datetime
import subprocess
import sys
import re 
sys.path.append(os.path.dirname(__file__))
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor

# returns reverse complement of a sequence
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

class CustomAdaptor:
    def __init__(self, adaptor_file, mobilogger):
        self.default_dict = {
            "M-RT-primer-1": {"seq":"GTGAGTGATGGTTGAGGTAGTGTGGAGNNNNNGGG", "cmd":"-A"},
            "M-RT-primer-2": {"seq":"GTGAGTGATGGTTGAGGTAGTGTGGAGNNNNNTTT", "cmd":"-A"}, 
            "M-RT-primer-rc": {"seq":"CTCCACACTACCTCAACCATCACTCAC", "cmd":"-A"}, 
            "Amplificatino-primer": {"seq":"AAGCAGTGGTATCAACGCAGAG", "cmd":"-A"}, 
            "Enrichment-primer-forward" : {"seq":"AGATCTGGAAGAGCGTCGTG", "cmd":"-a"}, 
            "P7": {"seq": "CAAGCAGAAGACGGCATACGAGAT", "cmd":"-a"}, 
            "P5": {"seq": "GTGTTAGATAGAGCCAGGCGGCCATGTT", "cmd":"-A"}, 
            "Truseq": {"seq": "AGATCGGAAGAG", "cmd":"-A"}, 
            "Nexetera": {"seq":"CTGTCTCTTATA", "cmd":"-A"}, 
            "Capture-sequence": {"seq":"AAAAAAAAAAAAAAACTCTGCGTTGATACCACTGCTT", "cmd":"-A"}, 
            "ploy-G1": {"seq": "G{10}", "cmd":"-A"}, 
            "ploy-T1": {"seq": "T{10}", "cmd":"-A"}, 
            "ploy-G1": {"seq": "G{10}", "cmd":"-a"}, 
            "ploy-T1": {"seq": "T{10}", "cmd":"-a"}
        }
        if adaptor_file != None:
            if not os.path.exists(adaptor_file):
                nd = mobilogger._mobilogrecorder(log_message="The custom adaptor file is not existed. %s" %(adaptor_file),
                log_level="ERROR")
                sys.exit()
            else:
                tmp_df = pd.read_csv(adaptor_file, sep="\t")
                expect_col = ["cmd", "seq", "remark"]
                for i in expect_col:
                    if not i in tmp_df.columns.tolist():
                        nd = mobilogger._mobilogrecorder(log_message="The required column %s is not existed in adaptor file %s." %(i, adaptor_file),
                        log_level="ERROR")
                        sys.exit()
                self.valid_dict = {}
                for i in tmp_df.index:
                    self.valid_dict[tmp_df.loc[i, "remark"]] = {"seq":tmp_df.loc[i, "seq"], "cmd":tmp_df.loc[i, "cmd"]}
                nd = mobilogger._mobilogrecorder(log_message="Using %s as adaptors." %(adaptor_file),
                log_level="INFO")
        else:
            self.valid_dict = self.default_dict

def str2bool(x):
    return x.lower() in ('true')

class BacPreProcess(object):
    def __init__(self, 
                fastq_path,
                R1_file,
                R2_file,
                output_path, 
                threads, 
                sample_name, 
                pre_process, 
                process_cutadapt, 
                process_fastp, 
                adaptor_list, 
                host_remove, 
                host_ref, 
                download_path, 
                fastp_adapter, 
                qualified_quality_phred, 
                unqualified_percent_limit, 
                n_base_limit, 
                length_required, 
                lib_type, 
                with_CB, 
                star_config, 
                dev_mod, 
                mobilogger):
        self.fastq_path = fastq_path
        self.R1_file = R1_file
        self.R2_file = R2_file
        self.output_path = output_path
        self.threads = threads
        self.sample_name = sample_name
        self.pre_process = pre_process
        self.run_cutadaptor = process_cutadapt
        self.run_fastp = process_fastp
        self.host_remove = host_remove
        self.host_ref = host_ref
        self.download_path = download_path
        self.fastp_adapter = fastp_adapter
        self.qualified_quality_phred = qualified_quality_phred
        self.unqualified_percent_limit = unqualified_percent_limit
        self.n_base_limit = n_base_limit
        self.length_required = length_required
        self.lib_type = lib_type
        self.with_CB = with_CB
        self.star_config = star_config
        self.dev_mod = dev_mod
        self.mobilogger = mobilogger
        self.mobicommandlogger = MobiCommandLogSystem(o_dir=self.mobilogger.working_path, dev_mode=self.mobilogger.dev_mode)
        self.mobiexecutor = CommandExecutor(log_system=self.mobicommandlogger, console_output=False)
        self.adaptor_list = CustomAdaptor(adaptor_list, mobilogger).valid_dict

    def version(self):
        return 1.3
    
    def process_fastqc(self, QC_path:str, threads:int, fastq_path:str, \
                    input_R1:str, input_R2:str, process_R2:bool, \
                    process_type:str):
        ###process raw data fastqc
        raw_QC_path = QC_path + "/%s_QC" %(process_type)
        if not os.path.exists(raw_QC_path):
            os.mkdir(raw_QC_path)
        try_times = 0
        report_file = "null_device"
        while not os.path.exists(raw_QC_path + '/' + report_file) and try_times <= 5:
            if process_R2:
                cmd = "fastqc --noextract -t %s -o %s %s %s" %(str(threads), raw_QC_path, fastq_path + "/" + input_R1, fastq_path + "/" + input_R2)
            else:
                cmd = "fastqc --noextract -t %s -o %s %s " %(str(threads), raw_QC_path, fastq_path + "/" + input_R1)
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "raw-data-qc"}, console_output=False)
            try_times += 1
            temp_file_list = os.listdir(raw_QC_path)
            for z in range(len(temp_file_list)):
                if ".html" in temp_file_list[z]:
                    report_file = temp_file_list[z]
        if not try_times <= 5:
            nd = self.mobilogger._mobilogrecorder(log_message="Run fastqc on raw data failed more than 5 times, plz recheck.",
                log_level="ERROR")
            return 0, 0, False
        for z in range(len(temp_file_list)):
            if process_R2:
                if ".zip" in temp_file_list[z] and "R1" in temp_file_list[z]:
                    temp_R1_file = temp_file_list[z]
            else:
                if ".zip" in temp_file_list[z]:
                    temp_R1_file = temp_file_list[z]
        R1_zip = raw_QC_path + '/' + temp_R1_file
        cmd = "unzip -o " + R1_zip + " -d " + raw_QC_path
        exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "raw-data-qc1-unzip"},console_output=False)
        if exit_code != 0:    
            nd = self.mobilogger._mobilogrecorder(log_message="Run unzip on raw data R1 QC result failed.",
                log_level="ERROR")
        except_R1_qc_file = R1_zip.replace(".zip","") + "/fastqc_data.txt"
        if not os.path.exists(except_R1_qc_file):
            nd = self.mobilogger._mobilogrecorder(log_message="Result files of Fastqc R1 is not found, plz recheck.",
                log_level="ERROR")
        else:
            found_read_length = False
            Raw_reads = 0
            Raw_bases = 0
            with open(except_R1_qc_file, "r") as f:
                while True:
                    line = f.readline()
                    if not line:
                        break
                    line = line.replace("\n","")
                    if "Total Sequences" in line:
                        Raw_reads = int(line.split("\t")[1])
                    if found_read_length:
                        if not ">>END_MODULE" in line:
                            temp_len = line.split("\t")[0]
                            if "-" in temp_len:
                                temp_mean_list = [int(temp_len.split("-")[0]), int(temp_len.split("-")[1])]
                                temp_mean = np.mean(temp_mean_list)
                            else:
                                temp_mean = int(temp_len)
                            Raw_bases += temp_mean * eval(line.split("\t")[1])
                        else:
                            found_read_length = False
                            break
                    if "#Length" in line:
                        found_read_length = True
        if process_R2:
            for z in range(len(temp_file_list)):
                if ".zip" in temp_file_list[z] and "R2" in temp_file_list[z]:
                    temp_R2_file = temp_file_list[z]
            R2_zip = raw_QC_path + '/' + temp_R2_file
            cmd = "unzip -o " + R2_zip + " -d " + raw_QC_path
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "raw-data-qc2-unzip"},console_output=False)
            if exit_code != 0:    
                nd = self.mobilogger._mobilogrecorder(log_message="Run unzip on raw data R2 qc failed.",
                    log_level="ERROR")
            except_R2_qc_file =  R2_zip.replace(".zip","") + "/fastqc_data.txt"
            if not os.path.exists(except_R2_qc_file):
                nd = self.mobilogger._mobilogrecorder(log_message="Result files of Fastqc R1 is not found, plz recheck.",
                    log_level="ERROR")
            else:
                found_read_length = False
                with open(except_R2_qc_file, "r") as f:
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        line = line.replace("\n","")
                        if found_read_length:
                            if not ">>END_MODULE" in line:
                                temp_len = line.split("\t")[0]
                                if "-" in temp_len:
                                    temp_mean_list = [int(temp_len.split("-")[0]), int(temp_len.split("-")[1])]
                                    temp_mean = np.mean(temp_mean_list)
                                else:
                                    temp_mean = int(temp_len)
                                Raw_bases += temp_mean * eval(line.split("\t")[1])
                            else:
                                found_read_length = False
                                break
                        if "#Length" in line:
                            found_read_length = True
        return Raw_reads, Raw_bases, True

    def process_cutadaptor(self, adaptor_list: dict, QC_path:str, threads:int, \
                        output_head:str, fastq_path:str, input_R1:str, \
                            input_R2:str, process_R2:bool):
        cmd_adaptor = ""
        output_path_cutadapt = QC_path + "/cutadapt_out"
        if not os.path.exists(output_path_cutadapt):
            os.mkdir(output_path_cutadapt)
        for a in adaptor_list.keys():
            cmd_adaptor += str(adaptor_list[a]["cmd"] + " " + adaptor_list[a]["seq"] + " ")
        if process_R2:
            cmd = "cutadapt -j %s -m 30 %s-n 5 --json %s -o %s -p %s %s %s" %(str(threads), cmd_adaptor, output_path_cutadapt + "/cutadapt_report.json", 
                                                                    output_path_cutadapt + "/" + output_head + "_combined_clean1_R1.fastq.gz",
                                                                    output_path_cutadapt + "/" + output_head + "_combined_clean1_R2.fastq.gz",
                                                                    fastq_path + "/" + input_R1, fastq_path + "/" + input_R2)
        else:
            cmd = "cutadapt -j %s -m 30 %s-n 5 --json %s -o %s %s " %(str(threads), cmd_adaptor, output_path_cutadapt + "/cutadapt_report.json",
                                                            output_path_cutadapt + "/" + output_head + "_combined_clean1_R1.fastq.gz",
                                                            fastq_path + "/" + input_R1)
        try_times = 0
        while not os.path.exists(output_path_cutadapt + "/" + output_head + "_combined_clean1_R1.fastq.gz") and try_times <= 5:
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "cutadapt"},console_output=False)
            try_times += 1
        if try_times <= 5:
            input_R1 = output_head + "_combined_clean1_R1.fastq.gz"
            input_R2 = output_head + "_combined_clean1_R2.fastq.gz"
            fastq_path = output_path_cutadapt
        else:
            nd = self.mobilogger._mobilogrecorder(log_message="Run cutadapt failed.",
                log_level="ERROR")
            return "", "", "", False
        return input_R1, input_R2, fastq_path, True

    def process_fastp(self, QC_path:str, output_head:str, threads:int, fastp_adaptor_file:str, \
                    unqualified_percent_limit:int, qualified_quality_phred:int, \
                    n_base_limit:int, length_required:int, fastq_path:str, \
                    input_R1:str, input_R2:str, process_R2:bool):
        output_path_fastp = QC_path + "/fastp_out"
        if not os.path.exists(output_path_fastp):
            os.mkdir(output_path_fastp)
        except_out_R1 = output_path_fastp + "/" + output_head + "_combined_clean2_R1.fastq.gz"
        except_out_R2 = output_path_fastp + "/" + output_head + "_combined_clean2_R2.fastq.gz"
        if process_R2:
            cmd = "fastp --in1 %s --in2 %s --out1 %s --out2 %s -h %s -j %s -w %s " %(fastq_path + "/" + input_R1, fastq_path + "/" + input_R2,
                                                                                                            except_out_R1, except_out_R2, 
                                                                                                            output_path_fastp + "/" + output_head + "_fastp-report.html", 
                                                                                                            output_path_fastp + "/" + output_head + "_fastp-report.json", 
                                                                                                            str(threads)) + \
                "-u %s -q %s -n %s " %(unqualified_percent_limit, 
                                                                                            qualified_quality_phred,
                                                                                            n_base_limit)
        else:
            cmd = "fastp --in1 %s --out1 %s -h %s -j %s -w %s  " %(fastq_path + "/" + input_R1,
                                                                                                        except_out_R1, 
                                                                                                        output_path_fastp + "/" + output_head + "_fastp-report.html", 
                                                                                                        output_path_fastp + "/" + output_head + "_fastp-report.json", 
                                                                                                        str(threads)) + \
            "-u %s -q %s -n %s " %(unqualified_percent_limit, 
                                                                                        qualified_quality_phred,
                                                                                        n_base_limit)

        try_times = 0
        while not os.path.exists(except_out_R1) and try_times <= 5:
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "fastp"},console_output=False)
            try_times += 1
        if try_times <= 5:
            input_R1 = output_head + "_combined_clean2_R1.fastq.gz"
            input_R2 = output_head + "_combined_clean2_R2.fastq.gz"
            fastq_path = output_path_fastp
        else:
            nd = self.mobilogger._mobilogrecorder(log_message="Run fastp failed.",
                log_level="ERROR")
            return "", "", "", {}, False
        qc_stat = {}
        if process_R2:
            dd = 2
        else:
            dd = 1
        with open(output_path_fastp + "/" + output_head + "_fastp-report.json", "r") as f:
            json_data = json.load(f)
            qc_stat["Cutadapt_reads"] = json_data['summary']['before_filtering']['total_reads'] / dd
            qc_stat["Cutadapt_bases"] = json_data['summary']['before_filtering']['total_bases']
            qc_stat["Fastp_reads"] = json_data['summary']['after_filtering']['total_reads'] / dd
            qc_stat["Fastp_bases"] = json_data['summary']['after_filtering']['total_bases']
            qc_stat["Q20"] = json_data['summary']['before_filtering']['q20_rate']
            qc_stat["Q30"] = json_data['summary']['before_filtering']['q30_rate']
            qc_stat["GC"] = json_data['summary']['before_filtering']['gc_content']
        return input_R1, input_R2, fastq_path, qc_stat, True

    def rm_adaptor(self, input_R1, input_R2, sample_barcode, sub_process_data, expect_cut, process_R2, ii, fastq_path):
        output_head = self.sample_name
        qc_data = pd.DataFrame(columns=["sample_ID", "sample_barcode", "Raw_reads", "Raw_bases", "Cutadapt_reads", "Cutadapt_bases", \
                                        "Fastp_reads", "Fastp_bases", "Host_unremoved_reads", "Host_unremoved_bases", "Q20", "Q30", "GC"])
        qc_data.loc[0, "sample_ID"] = output_head
        qc_data.loc[0, "sample_barcode"] = sample_barcode
        pre_output_path = self.output_path + "/pre-process/" + sample_barcode
        if not os.path.exists(pre_output_path + "/qc_data.tsv") or not os.path.exists(pre_output_path + "/fastp_out/" + output_head + "_combined_clean2_R2.fastq.gz"):
            if not os.path.exists(pre_output_path):
                os.makedirs(pre_output_path)
            if self.dev_mod:
                nd = self.mobilogger._mobilogrecorder(log_message="Checking raw data quality...",
                    log_level="INFO")
                qc_data.loc[0, "Raw_reads"], qc_data.loc[0, "Raw_bases"], process_stat = \
                    self.process_fastqc(pre_output_path, expect_cut, fastq_path, \
                    input_R1, input_R2, process_R2, process_type="Raw")
                if process_stat:
                    nd = self.mobilogger._mobilogrecorder(log_message="Raw data quality checking is done successfully.",
                        log_level="INFO")
                else:
                    nd = self.mobilogger._mobilogrecorder(log_message="Raw data quality checking failed.",
                        log_level="ERROR")
                    sys.exit()
            else:
                process_stat = True
            if not process_stat:
                return qc_data
            ###process cutadapt
            ###note start time of cutadapt
            if self.run_cutadaptor:
                now = datetime.datetime.now()
                formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
                current_add = sub_process_data.shape[0]
                sub_process_data.loc[current_add, "start"] = formatted_date
                sub_process_data.loc[current_add, "process"] = "cutadapt" 
                nd = self.mobilogger._mobilogrecorder(log_message="Adapter Trimming...",
                    log_level="INFO")
                input_R1, input_R2, fastq_path, process_stat = self.process_cutadaptor(adaptor_list=self.adaptor_list, 
                                                                                QC_path=pre_output_path, \
                                                                                threads=expect_cut, \
                                                                                output_head=output_head, \
                                                                                fastq_path=fastq_path, \
                                                                                input_R1=input_R1, \
                                                                                input_R2=input_R2, \
                                                                                process_R2=process_R2)
                if process_stat:
                    nd = self.mobilogger._mobilogrecorder(log_message="Adapter Trimming is done successfully.",
                        log_level="INFO")
                else:
                    nd = self.mobilogger._mobilogrecorder(log_message="Adapter Trimming failed.",
                        log_level="ERROR")
                    sys.exit()
                if not process_stat:
                    return qc_data
            ###process fastp
            if self.run_fastp:
                nd = self.mobilogger._mobilogrecorder(log_message="Read Filtering ...",
                    log_level="INFO")
                input_R1, input_R2, fastq_path, qc_stat, process_stat = self.process_fastp(QC_path=pre_output_path, \
                                                                        output_head=output_head, \
                                                                        threads=expect_cut, \
                                                                        fastp_adaptor_file=self.fastp_adapter, \
                                                                        unqualified_percent_limit=self.unqualified_percent_limit, \
                                                                        qualified_quality_phred=self.qualified_quality_phred, \
                                                                        n_base_limit=self.n_base_limit, \
                                                                        length_required=self.length_required, \
                                                                        fastq_path=fastq_path, \
                                                                        input_R1=input_R1, \
                                                                        input_R2=input_R2, \
                                                                        process_R2=process_R2)
                if process_stat:
                    nd = self.mobilogger._mobilogrecorder(log_message="Read Filtering  is done successfully.",
                        log_level="INFO")
                else:
                    nd = self.mobilogger._mobilogrecorder(log_message="Read Filtering failed.",
                        log_level="ERROR")
                    sys.exit()
                if not process_stat:
                    return qc_data
                else:
                    qc_data.loc[0, "Cutadapt_reads"] = qc_stat["Cutadapt_reads"]
                    qc_data.loc[0, "Cutadapt_bases"] = qc_stat["Cutadapt_bases"]
                    qc_data.loc[0, "Fastp_reads"] = qc_stat["Fastp_reads"]
                    qc_data.loc[0, "Fastp_bases"] = qc_stat["Fastp_bases"]
                    qc_data.loc[0, "Q20"] = qc_stat["Q20"]
                    qc_data.loc[0, "Q30"] = qc_stat["Q30"]
                    qc_data.loc[0, "GC"] = qc_stat["GC"]
            ###host remove
            if self.host_remove:
                nd = self.mobilogger._mobilogrecorder(log_message="Removing host reads...",
                    log_level="INFO")
                output_path_host = pre_output_path + "/host_remove"
                if not os.path.exists(output_path_host):
                    os.mkdir(output_path_host)
                ref_path = self.host_ref
                if os.path.exists(ref_path + ".1.bt2"):
                    #print("ERROR. No host reference found, plz recheck.")
                    #return
                    cmd = "bowtie2 -x %s -U %s -S %s -p %s --un-gz %s 1> %s 2> %s" %(ref_path, fastq_path + "/" + input_R2, 
                                                                                output_path_host + "/" + output_head + ".sam", 
                                                                                str(expect_cut), 
                                                                                output_path_host + "/" + output_head, 
                                                                                output_path_host + "/" + output_head + "_align.error",
                                                                                output_path_host + "/" + output_head + "_align.out")
                    ##os.system(cmd)
                    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    cmd = "seqkit seq -i %s > %s" %(output_path_host + "/" + output_head, 
                                                    output_path_host + "/" + output_head + "_unmapped.txt")
                    ##os.system(cmd)
                    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    cmd = "sed -i 's/@//' %s" %(output_path_host + "/" + output_head + "_unmapped.txt")
                    ##os.system(cmd)
                    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    cmd = "seqkit grep -j %s -f %s -o %s %s" %(str(expect_cut), output_path_host + "/" + output_head + "_unmapped.txt",
                                                                output_path_host + "/" + output_head + "_combined_clean3_R1.fastq.gz", 
                                                                fastq_path + "/" + input_R1)
                    ##os.system(cmd)
                    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    cmd = "seqkit grep -j %s -f %s -o %s %s" %(str(expect_cut), output_path_host + "/" + output_head + "_unmapped.txt",
                                                                output_path_host + "/" + output_head + "_combined_clean3_R2.fastq.gz", 
                                                                fastq_path + "/" + input_R2)
                    ##os.system(cmd)
                    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                elif os.path.exists(os.path.join(ref_path, "star", "Genome")):
                    unstd = False
                    with gzip.open(os.path.join(fastq_path, input_R2)) as f:
                        line = f.readline().decode('utf8').replace("\n","")
                        if " " in line:
                            unstd = False
                            pass
                        elif line.endswith("/2"):
                            unstd = True
                            postfix1 = "\/1"
                            postfix2 = "\/2"
                    cmd = 'STAR --runThreadN {} --genomeDir {} --outFilterMultimapNmax 10 --quantMode GeneCounts ' \
                    '--outFilterScoreMin 30 --readFilesCommand zcat ' \
                    '--outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 65719476736 --readFilesIn {} ' \
                    '--outFileNamePrefix {} --outTmpDir {} --outReadsUnmapped {} --outBAMsortingThreadN {}'.format(
                    str(expect_cut) , os.path.join(ref_path, "star") ,  
                    os.path.join(fastq_path, input_R2), output_path_host + "/host_remove", os.path.join(output_path_host, "STAR_tmp"), 
                    "Fastx", str(expect_cut))
                    #subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    exit_code = self.mobiexecutor.execute(
                        command=cmd,
                        context={
                            "host-remove": "STAR"
                        },console_output=False
                    )
                    if exit_code != 0:
                        nd = self.mobilogger._mobilogrecorder(log_message="Host remove mapping failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                            log_level="ERROR")
                        sys.exit()
                    else:
                        nd = self.mobilogger._mobilogrecorder(log_message="Host remove mapping done successfully.",
                            log_level="INFO")
                    f_list = os.listdir(output_path_host)
                    for f in f_list:
                        if f.endswith('.out.mate1'):
                            cmd = "mv %s %s " %(os.path.join(output_path_host, f), os.path.join(output_path_host, output_head + "_combined_clean3_R0.fastq"))
                            #subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                            exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "mv"},console_output=False)
                            if exit_code != 0:
                                ns = self.mobilogger._mobilogrecorder(log_message="Host remove mv failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                                    log_level="ERROR")
                                sys.exit()
                            cmd = "pigz -p %s %s " %(str(expect_cut), os.path.join(output_path_host, output_head + "_combined_clean3_R0.fastq"))
                            #subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                            exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "pigz"},console_output=False)
                            if exit_code != 0:
                                ns = self.mobilogger._mobilogrecorder(log_message="Host remove pigz failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                                    log_level="ERROR")
                                sys.exit()
                    cmd = "zcat %s |awk 'NR %% 4 == 1' | cut -d' ' -f1 > %s " %(os.path.join(output_path_host, output_head + "_combined_clean3_R0.fastq.gz"), 
                                                       os.path.join(output_path_host, "unmapped_head.txt"))
                    #subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "unmap-head"},console_output=False)
                    if exit_code != 0:
                        ns = self.mobilogger._mobilogrecorder(log_message="Host remove unmap-head failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                            log_level="ERROR")
                        sys.exit()
                    cmd = "sed -i 's/@//' %s" %(os.path.join(output_path_host, "unmapped_head.txt"))
                    #subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
                    exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "sed"},console_output=False)
                    if exit_code != 0:
                        nd = self.mobilogger._mobilogrecorder(log_message="Host remove sed failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                            log_level="ERROR")
                        sys.exit()
                    if unstd:
                        mask_file1 = os.path.join(output_path_host, "unmapped_head1.txt")
                        cmd = "sed 's/$/%s/' %s > %s" %(postfix1, os.path.join(output_path_host, "unmapped_head.txt"), mask_file1)
                        exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "sed2"},console_output=False)
                        if exit_code != 0:
                            nd = self.mobilogger._mobilogrecorder(log_message="Host remove sed2 failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                                log_level="ERROR")
                            sys.exit()
                        mask_file2 = os.path.join(output_path_host, "unmapped_head2.txt")
                        cmd = "sed 's/$/%s/' %s > %s" %(postfix2, os.path.join(output_path_host, "unmapped_head.txt"), mask_file2)
                        exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "sed2"},console_output=False)
                        if exit_code != 0:
                            nd = self.mobilogger._mobilogrecorder(log_message="Host remove sed2 failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                                log_level="ERROR")
                            sys.exit()
                    else:
                        mask_file1 = os.path.join(output_path_host, "unmapped_head.txt")
                        mask_file2 = os.path.join(output_path_host, "unmapped_head.txt")
                    cmd = "seqkit grep -j %s -f %s -o %s %s" %(str(expect_cut), mask_file1, 
                                            os.path.join(output_path_host, output_head + "_combined_clean3_R1.fastq.gz"), 
                                            os.path.join(fastq_path, input_R1))
                    exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "seqkit"},console_output=False)
                    cmd = "seqkit grep -j %s -f %s -o %s %s" %(str(expect_cut), mask_file2, 
                                            os.path.join(output_path_host, output_head + "_combined_clean3_R2.fastq.gz"), 
                                            os.path.join(fastq_path, input_R2))
                    exit_code = self.mobiexecutor.execute(command=cmd,context={"host-remove": "seqkit"},console_output=False)
                    if exit_code != 0:
                        nd = self.mobilogger._mobilogrecorder(log_message="Host remove seqkit failed. Check Logs/stderr.log or Logs/stdout.log for more information.",
                            log_level="ERROR")
                        sys.exit()
                input_R1 = output_head + "_combined_clean3_R1.fastq.gz"
                input_R2 = output_head + "_combined_clean3_R2.fastq.gz"
                fastq_path = output_path_host
                nd = self.mobilogger._mobilogrecorder(log_message="Host reads removing is done successfully.",
                    log_level="INFO")
            ###process clean data fastqc
            if self.dev_mod:
                nd = self.mobilogger._mobilogrecorder(log_message="Checking clean data quality...",
                    log_level="INFO")
                qc_data.loc[0, "Host_unremoved_reads"], qc_data.loc[0, "Host_unremoved_bases"], process_stat = \
                    self.process_fastqc(QC_path=pre_output_path, \
                                threads=expect_cut, \
                                fastq_path=fastq_path, \
                                input_R1=input_R1, \
                                input_R2=input_R2, \
                                process_R2=process_R2, \
                                process_type="clean")
                if process_stat:
                    nd = self.mobilogger._mobilogrecorder(log_message="Clean data quality checking is done successfully.",
                        log_level="INFO")
                else:
                    nd = self.mobilogger._mobilogrecorder(log_message="Clean data quality checking failed.",
                        log_level="INFO")
                    sys.exit()
            else:
                qc_data.loc[0, "Host_unremoved_reads"] = qc_data.loc[0, "Fastp_reads"]
                qc_data.loc[0, "Host_unremoved_bases"] = qc_data.loc[0, "Fastp_bases"] 
                process_stat = True
            #if not process_stat:
            #    return qc_data
            qc_data.to_csv(pre_output_path + "/qc_data.tsv", sep="\t", index=False)
        else:
            qc_data = pd.read_csv(pre_output_path + "/qc_data.tsv", sep="\t")
        return qc_data, fastq_path, ii

    def process(self):
        valid_mode = ['ploy-T', 'None', 'Nextera']
        expect_cut = self.threads
        fastq_path = self.fastq_path
        input_R1 = self.R1_file
        input_R2 = self.R2_file
        if not self.lib_type in ["Illumina", "Nanopore"]:
            print("Error: The --lib_type input is not valid.")
            return
        else:
            if self.lib_type == "Illumina":
                process_R2 = True
            if self.lib_type == "Nanopore":
                process_R2 = False

        output_head = self.sample_name
        if input_R1 == input_R2 and process_R2:
            nd = self.mobilogger._mobilogrecorder(log_message="R1 and R2 file path are the same. %s %s " %(input_R1, input_R2),
                log_level="ERROR")
            return
        sub_process_data = pd.DataFrame(columns=["start", "end", "process"])
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)
        ###note start time of barcode split
        now = datetime.datetime.now()
        formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
        current_add = sub_process_data.shape[0]
        sub_process_data.loc[current_add, "start"] = formatted_date
        sub_process_data.loc[current_add, "process"] = "barcoding" 
        output_path = self.output_path
        statis_output_path = self.output_path
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        ###process barcode
        if self.star_config["Mobivision-M"]["split_func"] == "GO":
            nd = self.mobilogger._mobilogrecorder(log_message="Barcoding...",
                    log_level="INFO")
            if self.dev_mod:
                dev_cmd = "-dev_mod=true" 
            else:
                dev_cmd = "-dev_mod=false"
            cmd = "%s -fo=%s " %(self.star_config["Mobivision-M"]["go_script"], fastq_path) + \
                "-o=%s -ff=%s -of %s -t=%s " %(output_path, self.star_config["Mobivision-M"]["filter_pattern_file"], self.star_config["Mobivision-M"]["output_pattern_file"], expect_cut) + \
                "-ot=%s -with_CB=%s -ID=%s " %(self.star_config["Mobivision-M"]["output_type"], str(self.with_CB), output_head) + \
                "-e=%s -k=%s %s " %(self.star_config["Mobivision-M"]["encrypt"], self.star_config["Mobivision-M"]["encrypt_key"], dev_cmd) + \
                "-seq_barcode_tag=%s " %(self.star_config["Mobivision-M"]["seq_barcode_tag"]) + \
                "-seq_UMI_tag=%s -split_by=%s" %(self.star_config["Mobivision-M"]["seq_UMI_tag"], self.star_config["Mobivision-M"]["split_by"])
            try_times = 0
            exit_code = 1
            while try_times < 3 and not os.path.exists(output_path + "/" + output_head + "_sample_barcode_stat.tsv") and exit_code != 65:
                #tmp_out = output_path + "/tmp_go_%s.out" %(str(try_times))
                if os.path.exists(output_path + "/split_fastq"):
                    cmd_tmp = "rm -r %s " %(output_path + "/split_fastq")
                    exit_code = self.mobiexecutor.execute(command=cmd_tmp,context={"pre-process": "rm"},console_output=False)
                #with open(tmp_out, "w") as f_o:
                #    subprocess.call(cmd, stdout=f_o, stderr=f_o, shell=True)
                exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "barcoding-%s" %(str(try_times))},console_output=False)
                try_times += 1
            if exit_code != 0:
                nd = self.mobilogger._mobilogrecorder(log_message="Barcoding failed.",
                    log_level="ERROR")
                sys.exit()
            else:
                nd = self.mobilogger._mobilogrecorder(log_message="Barcoding is done successfully. ",
                    log_level="INFO")
            sample_stat = pd.read_csv(output_path + "/" + output_head + "_sample_barcode_stat.tsv", sep="\t")
            #for ii in range(try_times):
            #    cmd = "rm %s " %(output_path + "/tmp_go_%s.out" %(str(ii)))
            #    exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "rm"})
            self.white_list_file = self.star_config["STAR"]["white_list_file"]
            if os.path.exists(self.white_list_file):
                #print("Replace current white list")
                self.mobilogger._mobilogrecorder(log_message="Using custom white list.",
                    log_level="INFO")
                cmd = "cp %s %s " %(self.white_list_file, 
                                    output_path + "/split_fastq/" + output_head + "/white_list.tsv")
                #os.system(cmd)
                exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "replace white-list"},console_output=False)
            ###change _1, _2 to _R1, _R2
            pattern = re.compile(r'([_.])(1|2)([_.])')
            for s in os.path.join(output_path, "split_fastq"):
                if not os.path.isfile(s):
                    continue          # 只处理文件
                new = pattern.sub(lambda m: f'{m.group(1)}R{m.group(2)}{m.group(3)}', s)
                if new != s:        # 真正发生替换才重命名
                    #print(f'{s}  ->  {new}')
                    os.rename(s, new)
        #cutadapt
        if self.pre_process:
            ###note start time of fastqc raw
            now = datetime.datetime.now()
            formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
            current_add = sub_process_data.shape[0]
            sub_process_data.loc[current_add, "start"] = formatted_date
            sub_process_data.loc[current_add, "process"] = "fastqc_raw_data"
            qc_data = pd.DataFrame(columns=["sample_ID", "sample_barcode", "Raw_reads", "Raw_bases", "Cutadapt_reads", "Cutadapt_bases", \
                            "Fastp_reads", "Fastp_bases", "Host_unremoved_reads", "Host_unremoved_bases", "Q20", "Q30", "GC"])

            #all_task = []
            #pool = multiprocessing.Pool(processes = 1)
            for i in sample_stat.index:
                if sample_stat.loc[i, "read_count"] >= 10000:
                    if self.star_config["Mobivision-M"]["split_func"] == "python":
                        input_R1 = output_head + "_" + sample_stat.loc[i, "sample_barcode"] + "_R1.fastq.gz"
                        input_R2 = output_head + "_" + sample_stat.loc[i, "sample_barcode"] + "_R2.fastq.gz"
                    elif self.star_config["Mobivision-M"]["split_func"] == "GO":
                        input_R1 = output_head + "_" + sample_stat.loc[i, "sample_barcode"] + "_S1_L001_R1_001.fastq.gz"
                        input_R2 = output_head + "_" + sample_stat.loc[i, "sample_barcode"] + "_S1_L001_R2_001.fastq.gz"   
                    temp_sample_barcode = sample_stat.loc[i, "sample_barcode"]
                    #temp_args = [input_R1, input_R2, temp_sample_barcode, sub_process_data, expect_cut, process_R2, i, sample_stat.loc[i, "fastq_path"]]
                    #all_task.append(pool.apply_async(self.rm_adaptor, temp_args))
                    stat = self.rm_adaptor(input_R1, input_R2, temp_sample_barcode, sub_process_data, expect_cut, process_R2, i, sample_stat.loc[i, "fastq_path"])
                    sample_stat.loc[i, "fastq_path_clean"] = stat[1]
                    qc_data = pd.concat([qc_data, stat[0]], axis=0)
                else:
                    sample_stat.loc[i, "fastq_path_clean"] = "NA"
                    temp_add = qc_data.shape[0]
                    qc_data.loc[temp_add, "sample_ID"] = output_head
                    qc_data.loc[temp_add, "sample_barcode"] = sample_stat.loc[i, "sample_barcode"]
                    qc_data.loc[temp_add, "Raw_reads"] = sample_stat.loc[i, "read_count"]
                    for j in qc_data.columns.tolist():
                        if pd.isna(qc_data.loc[temp_add, j]):
                            qc_data.loc[temp_add, j] = 0
        else:
            qc_data = pd.DataFrame(columns=["sample_ID", "Raw_reads", "Raw_bases", "Cutadapt_reads", "Cutadapt_bases", \
                                            "Fastp_reads", "Fastp_bases", "Q20", "Q30", "GC"], index = [0], data = 0)
        qc_data.to_csv(output_path + "/" + output_head + "_qc_stat.tsv", sep="\t", index=False)
        #sample_stat.to_csv(output_path + "/" + output_head + "_cleaned_sample_stat.tsv", sep="\t", index=False)
        ###note end time of fastqc clean data
        now = datetime.datetime.now()
        formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
        sub_process_data.loc[current_add, "end"] = formatted_date
        ###final output
        final_fastq_path = output_path + "/clean_data"
        if os.path.exists(final_fastq_path):
            cmd = "rm -r %s " %(final_fastq_path)
        os.makedirs(output_path + "/clean_data")
        for i in sample_stat.index:
            temp_sample_barcode = sample_stat.loc[i, "sample_barcode"]
            temp_clean_data_path = final_fastq_path + "/" + temp_sample_barcode
            os.makedirs(temp_clean_data_path)
            cmd = "mv %s %s " %(sample_stat.loc[i, "fastq_path_clean"] + "/" + output_head + "*_R1*.fastq.gz", 
                                temp_clean_data_path + "/" + output_head + "_" + temp_sample_barcode + "_S0_L001_R1_001.fastq.gz")
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "mv R1"},console_output=False)
            cmd = "mv %s %s " %(sample_stat.loc[i, "fastq_path_clean"] + "/" + output_head + "*_R2*.fastq.gz", 
                                temp_clean_data_path + "/" + output_head + "_" + temp_sample_barcode + "_S0_L001_R2_001.fastq.gz")
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "mv R2"},console_output=False)
            cmd = "cp %s %s " %(sample_stat.loc[i, "fastq_path"] + "/white_list.tsv", temp_clean_data_path + "/white_list.tsv")
            exit_code = self.mobiexecutor.execute(command=cmd,context={"pre-process": "mv white-list"},console_output=False)
        ###final
        with open(self.output_path + "/Job_done.flag", "w") as f:
            f.write("The pre_process process done succeccfully.")

        now = datetime.datetime.now()
        formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
        sub_process_data.loc[current_add, "end"] = formatted_date
        temp_name = self.output_path + "/sub_process_annnote.tsv"
        sub_process_data.to_csv(temp_name, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--fastq_path', type=str, help="The path of fastq files.")
    parser.add_argument('-o', '--output_path', type=str, help="The path to store splited fastq.gz.files.")
    parser.add_argument('-t', '--threads', type=int, help="The number of threads. Default value is 8.", default=8)
    parser.add_argument('-id', '--sample_name', type=str, help="Sample name use in output file name.")
    parser.add_argument('-pre', '--pre_process', type=str2bool, help="If pre-process should be run. Default valye is True.", default='True')
    parser.add_argument('-al', '--adaptor_list', type=str, help="The path of adaptor list. Only used in pre-process.")
    parser.add_argument('--host_remove', type=str2bool, help="If host should be removed. Default value is False.", default="False")
    parser.add_argument('--host_ref', type=str, help="The path of host reference. Only avaliable in host remove process. \
                        Default value is NA", default="NA")
    parser.add_argument('--download_path', type=str, help="The path to save import results. Defaul value is '', which means skip this process.", default="")
    parser.add_argument('--fastp_adapter', type=str, help="The path of fastp_adapter. \
                        Defaul value is /share/home/sc/Projects/microbe-seq_v0.1/fastp_adapter.fasta.", default="/share/home/sc/Projects/microbe-seq_v0.1/fastp_adapter.fasta")
    parser.add_argument('--qualified_quality_phred', type=int, help="Fastp argument. The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified", default=15)
    parser.add_argument('--unqualified_percent_limit', type=int, help="Fastp argument. how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40 percentage (int [=40])", default=40)
    parser.add_argument('--n_base_limit', type=int, help="Fastp argument. If one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])", default=5)
    parser.add_argument('--length_required', type=int, help="Fastp argument. Reads shorter than length_required will be discarded, default is 15. (int [=15])", default=15)
    parser.add_argument('--lib_type', type=str, help="The type of input library. Valid values are 'Illumina' and 'Nanopore'. Default value is 'Illumina", default="Illumina")
    parser.add_argument('--filter_file', type=str, help="The filter pattern file")
    parser.add_argument('--with_CB', type=str2bool, help="If CB include in lib.", default="True")
    parser.add_argument('--star_config', type=str, help="The path of extra argument in json. The default setting is 'NA'.", default='NA')
    parser.add_argument('--dev_mod', type=str2bool, help="If turn on the development mode.", default="False")
    args = parser.parse_args()
    print("All input argments are:")
    for k,v in vars(args).items():
        print(k,'=',v)
    with open(args.star_config, 'r', encoding='utf-8') as file:
        config_data = json.load(file)
    mp = BacPreProcess(fastq_path = args.fastq_path, 
                output_path=args.output_path, 
                threads=args.threads, 
                sample_name=args.sample_name, 
                pre_process=args.pre_process, 
                adaptor_list=args.adaptor_list, 
                host_remove=args.host_remove, 
                host_ref=args.host_ref, 
                download_path=args.download_path, 
                fastp_adapter=args.fastp_adapter, 
                qualified_quality_phred=args.qualified_quality_phred, 
                unqualified_percent_limit=args.unqualified_percent_limit, 
                n_base_limit=args.n_base_limit, 
                length_required=args.length_required, 
                lib_type=args.lib_type, 
                filter_file=args.filter_file, 
                with_CB = args.with_CB, 
                star_config = config_data,
                dev_mod=args.dev_mod)
    mp.process()