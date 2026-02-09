import os
import time
import sys
import json
import shutil
import gzip
import traceback
sys.path.append(os.path.dirname(__file__))
from mako_report import ExportReport
from dataexp import Data_InjectTool
from re_assign_multi import call_cell_by_threshold, call_cell_by_STAR
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor

if os.path.exists(os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m")):
    resource_dir = os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m", "resources")
else:
    resource_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources")

class ReCallCell:
    def __init__(self, analysis_dir:str, o_dir:str, call_mtx:str, filter_args:str, sample_ID:str, ver_kit:str, threads:int, nosecondary:bool, \
                 keep_bam:bool, keep_unmapped:bool, config:dict, multiplet_method:str, mobilogger: MobiLoggingSystem, run_cmd="demno", 
                 Temperature = 2):
        self.o_dir = o_dir
        if not mobilogger is None:
            self.mobilogger = mobilogger
        elif os.path.exists(self.o_dir):
            self.mobilogger = MobiLoggingSystem(o_dir=o_dir, dev_mode=False)
            self.mobilogger._mobilogrecorder(log_message="The output path is already existed. Won't overwrite.",
                log_level="ERROR")
            sys.exit()
        else:
            self.mobilogger = MobiLoggingSystem(o_dir=o_dir, dev_mode=False)
        self.mobicommandlogger = MobiCommandLogSystem(o_dir=self.mobilogger.working_path, dev_mode=self.mobilogger.dev_mode)
        self.mobiexecutor = CommandExecutor(log_system=self.mobicommandlogger, console_output=False)
        self.mobilogger._mobilogrecorder(log_message="Re-call microbes from an analysised resutl started.",
            log_level="INFO")
        if not os.path.exists(analysis_dir):
            self.mobilogger._mobilogrecorder(log_message="The analysis dir not found: %s. Plz recheck." %(analysis_dir),
                            log_level="ERROR")
            sys.exit()
        else:
            self.analysis_dir = analysis_dir
        if call_mtx == None:
            if os.path.exists(os.path.join(self.analysis_dir, "raw_re-assigned_cell_gene_matrix")):
                self.fetched_mtx_dir = os.path.join(self.analysis_dir, "raw_re-assigned_cell_gene_matrix")
            elif os.path.exists(os.path.join(self.analysis_dir, "raw_cell_gene_matrix")):
                self.fetched_mtx_dir = os.path.join(self.analysis_dir, "raw_cell_gene_matrix")
            else:
                self.mobilogger._mobilogrecorder(log_message= "No mtx fold found in %s. Plz assign with -c manually." %(self.analysis_dir),
                    log_level="ERROR")
                sys.exit()
        else:
            self.fetched_mtx_dir = os.path.join(self.analysis_dir, call_mtx)
            if not os.path.exists(self.fetched_mtx_dir):
                self.mobilogger._mobilogrecorder(log_message= "Mtx fold %s not found in %s. Plz recheck." %(call_mtx, self.analysis_dir),
                    log_level="ERROR")
                sys.exit()
        self.mobilogger._mobilogrecorder(log_message= "Using %s as raw mtx." %(self.fetched_mtx_dir),
            log_level="INFO")
        self.filter_mtx_dir = os.path.join(self.o_dir, "filtered_cell_gene_matrix/")
        self.tmp_path = os.path.join(o_dir, "_STARtmp2")
        expect_file_gz = os.path.join(self.analysis_dir, "raw_cell_gene_matrix", "matrix.mtx.gz")
        expect_file = os.path.join(self.analysis_dir, "raw_cell_gene_matrix", "matrix.mtx")
        if not os.path.exists(expect_file) and not os.path.exists(expect_file_gz):
            self.mobilogger._mobilogrecorder(log_message= "Mtx file not found in %s. Plz recheck." %(self.analysis_dir),
                log_level="ERROR")
            sys.exit()
        self.report_json = os.path.join(analysis_dir, "report.json")
        if not os.path.exists(self.report_json):
            self.mobilogger._mobilogrecorder(log_message="The report json not found in %s." %(self.report_json),
                log_level="ERROR")
            sys.exit()
        with open(self.report_json,'r') as load_f:
            self.report_json = json.load(load_f)
        self.filter_args = filter_args
        self.final_map_stats_flie = os.path.join(analysis_dir, 'map_stat.tsv')
        shutil.copy(self.final_map_stats_flie, os.path.join(self.o_dir, 'map_stat.tsv'))
        if sample_ID is None:
            try:
                self.sample_id = self.report_json["sample"]["id_ori"].strip()
            except IndexError:
                #prog_runlog(time.strftime("%Y-%m-%d %H:%M:%S\t", time.localtime()) + "The sample id doesn't found in report.json. The json file maybe correpted.")
                self.mobilogger._mobilogrecorder(log_message="The sample id doesn't found in report.json. The json file maybe correpted.",
                    log_level="ERROR")
                sys.exit()
        else:
            self.sample_id = sample_ID
        self.cell_stat_file = os.path.join(self.analysis_dir, "barcode_info.tsv.gz")
        if ver_kit != "":
            self.ver_kit = ver_kit
        else:
            if "ver_kit" in self.report_json.keys():
                self.ver_kit = self.report_json["ver_kit"]
            else:
                self.ver_kit = "Unknown"
        self.run_cmd = run_cmd
        self.threads = threads
        self.nosecondary = nosecondary
        self.keep_bam = keep_bam
        self.keep_unmapped = keep_unmapped
        self.config = config
        self.multiplet_method = multiplet_method
        os.makedirs(self.filter_mtx_dir)
        ###cp raw mtx
        tmp_base = os.path.basename(self.fetched_mtx_dir)
        shutil.copytree(self.fetched_mtx_dir, os.path.join(self.o_dir, tmp_base))
        if os.path.exists(os.path.join(self.analysis_dir, "raw_re-assigned_cell_gene_matrix")) and \
        not os.path.exists(os.path.join(self.o_dir, "raw_re-assigned_cell_gene_matrix")):
            shutil.copytree(os.path.join(self.analysis_dir, "raw_re-assigned_cell_gene_matrix"), os.path.join(self.o_dir, "raw_re-assigned_cell_gene_matrix"))
        if os.path.exists(os.path.join(self.analysis_dir, "raw_cell_gene_matrix")) and \
        not os.path.exists(os.path.join(self.o_dir, "raw_cell_gene_matrix")):
            shutil.copytree(os.path.join(self.analysis_dir, "raw_cell_gene_matrix"), os.path.join(self.o_dir, "raw_cell_gene_matrix"))
        #if os.path.exists(os.path.join(self.analysis_dir, "gene_type.json")):
        #    shutil.copy(os.path.join(self.analysis_dir, "gene_type.json"), self.o_dir)
        self.fetched_mtx_dir = os.path.join(self.o_dir, os.path.basename(self.fetched_mtx_dir))
        self.mobilogger = MobiLoggingSystem(o_dir=o_dir, dev_mode=False)
        self.Temperature = Temperature
        return
    
    def export_json(self, final_map_stats_flie, raw_mtx_dir, final_mtx, last_json, gene_type_file):
        pre_file = 'NA'
        Exp_test = Data_InjectTool(sample_ID=self.sample_id, 
                            output_path=self.o_dir,
                            pre_stat_file=pre_file,
                            input_stat_file='NA',
                            filter_stats_file='NA',
                            saturation_file=final_map_stats_flie, 
                            raw_mtx_path=raw_mtx_dir,
                            filter_mtx_path=final_mtx,
                            cell_stat_file=self.cell_stat_file, 
                            reference_json_path='NA',
                            kit=self.ver_kit, 
                            run_cmd=self.run_cmd, 
                            threads=int(self.threads),
                            run_cluster=(not self.nosecondary),
                            last_json=last_json,
                            multiplet_method=self.multiplet_method, 
                            gene_type_file=gene_type_file, 
                            dev_mod=False, 
                            Temperature = self.Temperature)
        report_json = Exp_test.export_json_microbe()
        return report_json

    def process(self):
        self.mobilogger._mobilogrecorder(log_message="Cell Calling by %s " %(self.filter_args),
            log_level="INFO")
        if not "min_UMI" in self.filter_args and not "min_reads" in self.filter_args:
            filter_mtx_dir, info_flag, info = call_cell_by_STAR(mtx_dir=self.fetched_mtx_dir, tmp_dir=self.tmp_path, filter_args=self.filter_args, filter_mtx_dir=self.filter_mtx_dir, mobiexecutor=self.mobiexecutor)
            self.mobilogger._mobilogrecorder(log_message=info, log_level=info_flag)
            if info_flag == "ERROR":
                sys.exit()
        else:
            ###filter by UMI or genes
            filter_mtx_dir = call_cell_by_threshold(mtx_dir=self.fetched_mtx_dir, cell_stat_file=self.cell_stat_file, filter_args=self.filter_args, filter_mtx_dir=self.filter_mtx_dir, threads=self.threads, mobiexecutor=self.mobiexecutor)
            self.mobilogger._mobilogrecorder(log_message="Cell Calling by hard filter is done successfully.", log_level="INFO")
        ###generate json
        self.mobilogger._mobilogrecorder(log_message="Making report... " ,
            log_level="INFO")
        try_times = 0
        done = False
        if not os.path.exists(os.path.join(self.analysis_dir, "gene_type.json")):
            gene_dict = {}
            with gzip.open(os.path.join(self.filter_mtx_dir, "features.tsv.gz"), "r") as f:
                for line in f:
                    tmp_line = line.decode('utf-8').replace("\n","")
                    tmp_gene = tmp_line.split("\t")[1]
                    gene_dict[tmp_gene] = {"type":"unknown","species":"unknown"}
            with open(os.path.join(self.o_dir, "gene_type.json"), "w") as json_file:
                json.dump(gene_dict, json_file, indent=4)
        else:
            cmd = "cp %s %s" %(os.path.join(self.analysis_dir, "gene_type.json"), os.path.join(self.o_dir, "gene_type.json"))
            os.system(cmd)
        while try_times <= 2:
            #p = multiprocessing.Process(target=self.export_json, args=(self.final_map_stats_flie, 
            #                                                            os.path.join(self.o_dir, "raw_cell_gene_matrix"), 
            #                                                            self.filter_mtx_dir, 
            #                                                            os.path.join(self.analysis_dir, "report.json"), 
            #                                                            os.path.join(self.o_dir, "gene_type.json")))
            #try_times += 1
            #p.start()
            #p.join()
            #if p.exitcode != 0:
            #    self.mobilogger._mobilogrecorder(log_message="Making report failed. Retry in 1 minute." ,
            #        log_level="WARNING")
            #    time.sleep(60)
            #else:
            #    self.mobilogger._mobilogrecorder(log_message="Making report is done successfully." ,
            #        log_level="INFO")
            #    done = True
            #    break
            try_times += 1
            try:
                tmp_json = self.export_json(self.final_map_stats_flie, os.path.join(self.o_dir, "raw_cell_gene_matrix"), 
                                            self.filter_mtx_dir, os.path.join(self.analysis_dir, "report.json"), 
                                            os.path.join(self.o_dir, "gene_type.json"))
            except Exception as e:
                self.mobilogger._mobilogrecorder(log_message="Making report failed. Retry in 1 minute." ,
                    log_level="WARNING")
                traceback.print_exc()
                time.sleep(60)
            else:
                self.mobilogger._mobilogrecorder(log_message="Making report is done successfully." ,
                    log_level="INFO")
                done = True
                break
        if not done:
            self.mobilogger._mobilogrecorder(log_message="Making report failed." ,
                log_level="ERROR")
        ###mako report
        p = ExportReport(template_file=os.path.join(resource_dir, "report", "report_template.html"), 
                    json_file=os.path.join(self.o_dir, "report.json"), 
                    output_file=os.path.join(self.o_dir, self.sample_id + "_report.html"),  
                    jquery=os.path.join(resource_dir, "report", "js", "jquery-latest.min.js"), 
                    plotly=os.path.join(resource_dir, "report", "js", "plotly-latest.min.js"),
                    favicon_file=os.path.join(resource_dir, "report", "image", "favicon.ico"), 
                    web_logo=self.config["Mobivision-M"]["logo_path"], 
                    web_back=os.path.join(resource_dir, "report", "image","back.png"))
        p.process()
        ###re-arrange other files
        i_list = os.listdir(self.analysis_dir)
        del_files = ["filtered_gene_type_stats.csv", "processed_gene_type.tsv", "raw_gene_type_stats.csv"]
        for i in i_list:
            if (i.endswith(".bam") or i.endswith(".bai")) and self.keep_bam:
                shutil.copyfile(os.path.join(self.analysis_dir, i), os.path.join(self.o_dir, i))
            elif i.endswith("fastq.gz") and self.keep_unmapped:
                shutil.copyfile(os.path.join(self.analysis_dir, i), os.path.join(self.o_dir, i))
        for i in os.listdir(self.o_dir):
            if i in del_files:
                os.remove(os.path.join(self.o_dir, i))
            elif i.endswith("_matrix") and os.path.isdir(os.path.join(self.o_dir, i)):
                tmp_dir = os.path.join(self.o_dir, i)
                j_list = os.listdir(tmp_dir)
                for j in j_list:
                    if j.endswith(".tsv") or j.endswith(".mtx"):
                        cmd = "gzip %s" %(os.path.join(tmp_dir, j))
                        exit_code = self.mobiexecutor.execute(command=cmd, context={"re_call_microbe": "gunzip"},console_output=False)
                    else:
                        os.remove(os.path.join(tmp_dir, j))
        self.mobilogger._mobilogrecorder(log_message="Re-call microbes from an analysised resutl finished." ,
                log_level="INFO")
        return

if __name__ == "__main__":
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', '--analysis_dir', type=str, help='The analysis results path.')
    parse.add_argument('-o', '--output_dir', type=str, help="The output path.")
    parse.add_argument('-f', '--filter_args', type=str, help="The filter argument.")
    parse.add_argument('-t', '--threads', type=int, help="The number of threads.", default=4)
    parse.add_argument('-s', '--sample_ID', type=str, help="The id of the result file user defined.")
    parse.add_argument('--nosecondary', help="Don't run the scanpy analysis.", default=False, action='store_true')
    parse.add_argument('--kit', type=str, help="Set the version of MobiMicrobe RNA-seq kit. The default setting is 'v1.0'.", default="v1.0")
    parse.add_argument('--keep_bam', default=False, action='store_true', help="Whether to keep the original bam files.")
    parse.add_argument('--keep_unmapped', default=False, action="store_true", help="Whether to keep the unmapped files.")
    args = parse.parse_args()
    p = ReCallCell(analysis_dir=args.analysis_dir, 
                   o_dir=args.output_dir, 
                   call_mtx=None, 
                   filter_args=args.filter_args, 
                   sample_ID=args.sample_ID, 
                   ver_kit=args.kit, 
                   threads=args.threads, 
                   nosecondary=args.nosecondary, 
                   keep_bam=args.keep_bam, 
                   keep_unmapped=args.keep_unmapped, 
                   mobilogger=None)
    p.process()