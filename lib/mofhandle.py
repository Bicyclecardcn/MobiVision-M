import shutil
import os
import stat
import subprocess
from glob import glob
import json
import csv

class file_rearrange:
    def __init__(self, raw_fdir, filter_mtx, raw_mtx, fetched_mtx, summary_dir, results_id, geneFeature_Output, output_path, keep_bam, 
                 dev_mod, keep_qc=True, matrix_files = ["barcodes.tsv", "features.tsv", "matrix.mtx", "UniqueAndMult-Uniform.mtx"]):
    
        self.real_rawd = os.path.abspath(raw_fdir)
        self.filter_dir = filter_mtx
        self.raw_dir = raw_mtx
        self.fetch_dir = fetched_mtx
        self.star_sf = summary_dir
        #self.result_dir = self.real_rawd + "/" + results_id + "_outs"
        self.result_dir = output_path
        self.sample_id  = results_id
        self.mtxfiles = matrix_files
        os.makedirs(self.result_dir, 0o755, exist_ok=True)
        
        self.gene_featureDir = geneFeature_Output
        self.real_dir = os.path.abspath(self.result_dir)
        self.keep_bam = keep_bam
        self.dev_mod = dev_mod
        self.keep_qc = keep_qc
        
    def fileMove(self):
        self.interfile_deletion()
        o_res, e_res = self.resultfile_move()
    
    def dev_move(self, dir1, dir2, file):
        if os.path.exists(os.path.join(dir1, file)):
            if self.dev_mod:
                shutil.move(os.path.join(dir1, file), dir2)
            else:
                os.remove(os.path.join(dir1, file))

    def interfile_deletion(self):
        os.chdir(self.real_rawd)
        mapdir_files = self.list_file_fullpath("map_result")
        ## star align Summary.csv file will be stored in star_sf
        up_dir = os.path.dirname(self.real_rawd)
        pre_process_dir = os.path.join(up_dir, 'pre_process')
        star_tmp_dir = os.path.join(up_dir, '_STARtmp')
        if os.path.exists(star_tmp_dir):
            cmd = "rm -r %s " %(star_tmp_dir)
            os.system(cmd)
        if self.dev_mod:
            dev_path = self.real_rawd + "/dev_files"
            if not os.path.exists(dev_path):
                os.makedirs(dev_path)
            if os.path.exists(pre_process_dir):
                cmd = "mv %s %s " %(pre_process_dir, dev_path)
                os.system(cmd)
        elif self.keep_qc:
            shutil.move(pre_process_dir, self.result_dir)
        else:
            if os.path.exists(pre_process_dir):
                cmd = "rm -r %s " %(pre_process_dir)
                os.system(cmd)
        for f in mapdir_files:
            if f.endswith("unique.bam") or f.endswith(".bai"):
                if self.keep_bam:
                    shutil.copy(f, self.result_dir)
                else:
                    os.remove(f)
            elif f.endswith(".out") and not os.path.isdir(f):
                if self.dev_mod:
                    shutil.copy(f, dev_path)
            elif f.endswith(".out.mate1"):
                shutil.copy(f, self.result_dir + "/Unmapped_R2.fastq")
                cmd = "gzip %s " %(self.result_dir + "/Unmapped_R2.fastq")
                os.system(cmd)
            elif f.endswith(".out.mate2"):
                shutil.copy(f, self.result_dir + "/Unmapped_R1.fastq")
                cmd = "gzip %s " %(self.result_dir + "/Unmapped_R1.fastq")
                os.system(cmd)
            elif f.endswith("map_stat.tsv"):
                shutil.copy(f, self.result_dir)
            ###BacDrop keep rRNA
            ##elif f.endswith("gene_type_stats.csv"):
            ##    shutil.copy(f, self.result_dir)
            ##elif f.endswith("stat.tsv"):
            ##    shutil.copy(f, self.result_dir)
            elif self.dev_mod:
                if f.endswith("bam"):
                    shutil.copy(f, dev_path)
                elif f.endswith(".csv") or f.endswith(".tsv"):
                    shutil.copy(f, dev_path)
                elif os.path.isdir(f):
                    shutil.copytree(f, dev_path + "/" + f.split("/")[-1])
                    #cmd = "mv %s %s " %(f, dev_path)
                    #os.system(cmd)
                else:
                    shutil.copy(f, dev_path)
        shutil.copytree(self.filter_dir, os.path.join(self.result_dir, "filtered_cell_gene_matrix"))
        shutil.copytree(self.raw_dir, os.path.join(self.result_dir, "raw_cell_gene_matrix")) 
        os.chmod(os.path.join(self.result_dir, "filtered_cell_gene_matrix") , stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH + stat.S_IXOTH)
        os.chmod(os.path.join(self.result_dir, "raw_cell_gene_matrix") , stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH + stat.S_IXOTH)
        if self.fetch_dir != None:
            shutil.copytree(self.fetch_dir, os.path.join(self.result_dir, "raw_re-assigned_cell_gene_matrix")) 
            os.chmod(os.path.join(self.result_dir, "raw_re-assigned_cell_gene_matrix") , stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH + stat.S_IXOTH)
        sumdir_files = self.list_file_fullpath("summary_data")
        for f in sumdir_files:
            if f.endswith(".html") or f.endswith("h5ad"):
                shutil.copy(f, self.result_dir)
            elif f.endswith(".json"):
                ## generate sample_id + summary.csv file
                self.convert_csv2(f,  self.sample_id + "_summary.csv" , self.star_sf)
                shutil.move(self.sample_id + "_summary.csv", self.result_dir)
                shutil.copy(f, self.result_dir)
            elif self.dev_mod:
                if f.endswith("bam"):
                    shutil.copy(f, dev_path)
                elif f.endswith(".csv") or f.endswith(".tsv"):
                    shutil.copy(f, dev_path)
                elif os.path.isdir(f):
                    shutil.copy(f, dev_path)
        if os.path.exists("run_analysis_cmds.txt"):
            shutil.move("run_analysis_cmds.txt", self.result_dir)
        tmp_files = os.listdir(up_dir)
        for i in tmp_files:
            if i.endswith('fastq.del.gz'):
                os.remove(os.path.join(up_dir , i))
        check_list = ["tmp.log", "processed_gene_type.tsv", "clean_score.tsv", "gene_correlation.tsv", "gene_score.tsv"]
        for i in check_list:
            if self.dev_mod:
                self.dev_move(self.result_dir, dev_path, i)
            else:
                tmp_dir = os.path.join(self.result_dir, i)
                if os.path.exists(tmp_dir):
                    os.remove(tmp_dir)
            
    def resultfile_move(self):
        out ,err = "complete", ""        

        os.chdir(self.result_dir)
        
        test_gzres = self.check_and_gzip("filtered_cell_gene_matrix")
        test_gzres = self.check_and_gzip("raw_cell_gene_matrix")
        if os.path.exists(os.path.join(self.result_dir, "raw_re-assigned_cell_gene_matrix")):
            test_gzres = self.check_and_gzip("raw_re-assigned_cell_gene_matrix")

        f_bam = glob("*.out.bam")
        f_bai = glob("*out.bam.bai")
        
        os.chdir("../")
        shutil.rmtree("map_result")
        shutil.rmtree("summary_data")
        return out, err
        
    def check_and_gzip(self, dir_name):
        
        out_s = {}
        files = self.list_file_fullpath(dir_name)
        for f in files:
            basef = os.path.basename(f)
            if basef in self.mtxfiles:
                ret_co = subprocess.check_call("gzip %s"%(f), shell=True)
                out_s[f] = ret_co
            else:
                os.remove(f)
                
        # print(out_s)
        return out_s
     
    def list_file_fullpath(self, filedir):
        allf = os.listdir(filedir)

        resL = []
        for f in allf:
            resL.append(os.path.join(filedir, f))

        return resL 

    def convert_csv2(self, in_json, outfile , star_summaryfile, special_value=""):
        '''
         convert report json data to summary csv result
        '''
        
        frac_unique = 0
        with open(star_summaryfile, "r") as sum_fh:
            for line in sum_fh:
                if line.startswith("Fraction of Unique Reads in Cells"):
                    frac_unique = round(float(line.rstrip().split(",")[1]) * 100.0, 2)
    
        with open(in_json, 'r', encoding="utf-8") as fh:
            f_info = json.load(fh)
        
        rawInfo = []
        out_dict  = {}
        ## get_dict_allkeys
    
        for keyname, v in f_info.items():
            if keyname in ["mapping_table_data", "cell_table_data", "sample", "summary_tab"]:
                if isinstance(v,dict):
                    indict_a = v
                    for key, value in indict_a.items():
                        if isinstance(value,str) and value != special_value and key != "description":
                            rawInfo.append(key)
                            out_dict[key] = value
                        elif isinstance(value, dict):
                            for i in value:
                                if value[i] != special_value:
                                    rawInfo.append(i)
                                    out_dict[i]= value[i]
                elif isinstance(v, list):
                    inlist_a = v
                    for i in inlist_a:
                        if isinstance(i, str):
                            continue
                        #print(i)
                        elif isinstance(i, dict):
                            if i["name"] == "Median Genes per Cell" :
                                rawInfo.append(i["name"])
                                out_dict[i["name"]] =  i["metric"]
                    
                else:
                    continue
    
        fieldInfo = []
        for i in rawInfo:
            if i not in fieldInfo:
                fieldInfo.append(i)
                
        ## add information of "Fraction of Unique Reads in Cells"
        fieldInfo.append("Fraction of Unique Reads in Cells")        
        out_dict["Fraction of Unique Reads in Cells"]  =  str(frac_unique)
        
        for i in out_dict:
            tmp_v = out_dict[i]
            out_dict[i] = tmp_v.replace("\\n", "")
    
        with open( outfile, 'w',  encoding="utf-8", newline='') as csvfile:
        #写入表头名称
            writer = csv.DictWriter(csvfile, delimiter=',', fieldnames=fieldInfo )
            writer.writeheader()
        # 单行写入
            writer.writerow(out_dict)


if __name__ == "__main__":
    import argparse

    parse = argparse.ArgumentParser(prog = 'Demo')

    parse.add_argument('-d', '--origindir', help='orginal directory store convert_fastq ...')
    parse.add_argument("--filter_mtx", type=str, help="The path of filterd mtx.")
    parse.add_argument("--raw_mtx", type=str, help="The path of raw mtx.")
    parse.add_argument("--summary_file", type=str, help="The path of summary_dir.")
    parse.add_argument('-i', '--in_id', help='input id of sample')
    parse.add_argument('--intron_option', type=str, help="The intron option")
    parse.add_argument('-o', '--output_path', type=str, help="The output path")
    parse.add_argument('--dev_mod', type=str, help="If activate development mode")

    args = parse.parse_args()
    test_check = file_rearrange(args.origindir, args.filter_mtx, args.raw_mtx, args.summary_dir, args.in_id, args.intron_option, args.output_path, args.dev_mod)

    test_check.fileMove()
    

