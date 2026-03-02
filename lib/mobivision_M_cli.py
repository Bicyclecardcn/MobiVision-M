#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# @Date:  2024-0201 09:37
# @Author: sc
# @Filename: mobivis_CLI_M
# @Develop env: vscode

import os
import sys
import argparse
import textwrap
sys.path.append(os.path.dirname(__file__))
from mobivisionlogging import MobiLoggingSystem

def find_duplicates(lst):
    duplicates = []
    count_dict = {}
    for i in lst:
        if i in count_dict:
            count_dict[i] += 1
            duplicates.append(i)
        else:
            count_dict[i] = 1
    return duplicates

def find_exist(lst):
    return_list = []
    for i in lst:
        if not os.path.exists(i):
            return_list.append(i)
    return return_list

class ArgSetting:
    def __init__(self, versionN, log_file):
        self.version_num = versionN
        self.log_name = log_file
        self.args = None

    def get_all_args(self):
        parser = argparse.ArgumentParser(prog='', add_help=False, usage=argparse.SUPPRESS, description="",
                                         formatter_class=self.CustomHelpFormatter, epilog=textwrap.dedent('\n'))
        parser._optionals.title = "FLAGS"
        # define the title for main subcommands of MobiVision-M
        subparsers = parser.add_subparsers(dest="subcmds", metavar="", title="SUBCOMMANDS")
        # quantify
        scmd_quan = subparsers.add_parser('quantify', help='Count gene expression reads from a single sample',add_help=False,
                                          usage=argparse.SUPPRESS)
        scmd_quan._optionals.title = 'OPTIONS'
        scmd_quan.add_argument(
             '-f', '--fastqDir', 
             dest="input_dir", 
             action='store',
             help="The directory of fastq files.")
        scmd_quan.add_argument(
             '-t', '--threads', metavar='<NUM>', 
             dest='core_num', 
             default="12", 
             action='store',
             help='The number of threads.')
        scmd_quan.add_argument(
             '-i', '--index_dir', 
             action='store',
             help='Path of folder containing MobiVision-M-compatible transcriptome reference.')
        scmd_quan.add_argument(
             '-o', '--ouput_dir', 
             action='store', 
             dest="output_dir",
             help='The path of the output files. Use the current working directory as output path if not assigned.', 
             default=None)
        scmd_quan.add_argument(
             '-s', '--sample_ID', 
             default=None, 
             action='store', 
             help='The id of the result file user defined.')
        scmd_quan.add_argument(
             '--cellnumber',
             help='Force cell number for cell filter.',
             dest="topc",
             type=int, 
             default=None
        )    ## default=3000,
        scmd_quan.add_argument(
             '--cr2.2',
             help='Set CellRanger2.2 algorithm for cell filter. If not designated, the EmptyDrops algorithm will be used.',
             default=False,
             dest="empty_cr_params",
             action='store_true'
        )    ## defalut use EmptyDrops_CR , otherwise CellRanger2.2
        scmd_quan.add_argument(
            '--hard_filter', 
            help="Use a minium gene or UMI threshold to call microbes.", 
            default=None
        )
        scmd_quan.add_argument("--temperature", 
                        help="The Temperature of softmax.", 
                        default=2, 
                        type=int)
        scmd_quan.add_argument(
             '--kit', metavar='<Kit Version>',
             help="Set the version of scRNA-seq kit. The default setting is 'Unknown'.",
             default="v1.0",
             dest="kitV",
             type=str
        )
        scmd_quan.add_argument(
            '--config',
            help="",
            default="NA",
            dest="config",
            type=str
        )
        scmd_quan.add_argument(
            '--ksmiaa', 
            '--dev', 
            dest='ksmiaa', 
            help="Enable development mode (same as --ksmiaa)", 
            default=False, 
            action='store_true'
        )
        scmd_quan.add_argument(
            '--UMI_adjust', 
            help="The method of UMI adjust.", 
            type=str, 
            default="no_adjust", 
            choices=["no_adjust", "step_1", "step_1_and_2"]
        )
        scmd_quan.add_argument(
            "--test_run", 
            help="Run small demo data of quantify.", 
            default=False, 
            action='store_true'
        )
        scmd_quan.add_argument(
            "--nosecondary", 
            help="Don't run the scanpy analysis.", 
            default=False, 
            action='store_true'
        )
        scmd_quan.add_argument(
            "--keep_bam", 
            help="Keep the bam file.", 
            default=False, 
            action='store_true'
        )
        scmd_quan.add_argument(
            '--keep_unmap_reads', 
            help="Keep the unmapped reads.",
            default=False,
            action='store_true'
        )
        scmd_quan.add_argument(
            '--multiplet_method', 
            help="The method to detect multiplet.", 
            choices=["scaled_softmax", "majority", "auto"], 
            type=str, 
            default="auto"
        )
        scmd_quan.add_argument(
            "--qc_only", 
            help="Only run QC process.", 
            default=False, 
            action='store_true'
        )
        scmd_quan.add_argument(
            "--host_remove", 
            help="If remove host reads.",
            default=False, 
            action='store_true'
        )
        scmd_quan.add_argument(
            "--host_reference", 
            help="The path of host reference.", 
            default="NA"
            )
        scmd_quan.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                                  help="Show this information.")
        # mkindex
        scmd_mkindex = subparsers.add_parser('mkindex', add_help=False,
                                             help='Prepare a reference for MobiVision-M software',
                                             usage=argparse.SUPPRESS)
        scmd_mkindex._optionals.title = 'OPTIONS'
        scmd_mkindex.add_argument('-n', '--nameOfSpecies', metavar='<GENOME>', dest="genome", action='append', 
                                  help="Unique genome name [a-zA-Z0-9_]+. If the index is created with multiple genomes, Plz specify the -n argument multiple times. For example, if one specifies -n <genome1> and another -n <genome2>, the output index folder will be named <genome1>_and_<genome2>.")
        scmd_mkindex.add_argument('--input_file',  metavar='<INPUT_FILE>', dest="input_file", action='store', default=None, 
                                  help="The Tab-separated txt file contains the list of fasta and gtf file.")
        scmd_mkindex.add_argument('-f', '--fastaPath', metavar='<FASTA>', dest="fasta", action='append', 
                                  help="Genome fasta file path. If the index is created with multiple FASTA files, Plz specify the -f argument multiple times. The supported formats of genome fasta files include '.fasta', '.fa.gz' and '.fna.gz'.")
        scmd_mkindex.add_argument('-g', '--gtfPath', metavar='<GTF>', dest="gtf", action='append', 
                                  help="GTF file path for genome. If the index is created with multiple GTFs, Plz specify the -g argument multiple times. If -n and -f arguments were specified multiple times, -g argument needs to be specified an equal amount of times and the specified order must be same.")

        scmd_mkindex.add_argument('-r', '--referenceVerString', metavar='<REF_VERSION>', dest="refv", action='store',
                                  help="The version of output genome index, e.g. \"ref_v1.0\", \"ref-2022-06\".")
        scmd_mkindex.add_argument('-m', '--memoryUsed',  metavar='<MEM_GB>', dest="memgb", action='store', default=64, type=int,
                                  help="Maximum memory (GB) used when building index files with STAR.")
        scmd_mkindex.add_argument('-o', '--output_dir', type=str, help="The output name.", default=None)
        scmd_mkindex.add_argument(
            '--ksmiaa', 
            '--dev',
            dest='ksmiaa',
            default=False, 
            action='store_true', 
            help='Enable development mode (same as --ksmiaa)'
        )
        scmd_mkindex.add_argument(
            "--test_run", 
            help="Run small demo data of quantify.", 
            default=False, 
            action='store_true'
        )
        scmd_mkindex.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                                  help="Show this information.")
        #re-call cell
        scmd_rcell = subparsers.add_parser('rcmicrobe', add_help=False,
                                             help='re-count microbe form a MobiVision-M result.',
                                             usage=argparse.SUPPRESS)
        scmd_rcell.add_argument('-i', '--analysis_dir', type=str, help="The path of analysised MobiVision-M results.")
        scmd_rcell.add_argument('-o', '--output_dir', type=str, help="The path of output.", default=None)
        scmd_rcell.add_argument('-c', '--call_mtx', type=str, help="Re-call microbes by which mtx.", default=None)
        scmd_rcell.add_argument(
            '-t', '--threads', metavar='<NUM>', 
            dest='core_num', 
            default="12", 
            action='store',
            help='The number of threads.')
        scmd_rcell.add_argument(
             '--kit', 
             help="Set the version of scRNA-seq kit. The default setting is 'Unknown'.",
             default="v1.0",
             type=str
        )
        scmd_rcell.add_argument(
            '-s', '--sample_ID', 
            default=None, 
            action='store', 
            help='The id of the result file user defined.')
        scmd_rcell.add_argument(
            '--cellnumber', metavar='<TOP_CELL>',
            help='Force cell number for cell filter.',
            dest="topc",
            type=int, 
            default=None
        )
        scmd_rcell.add_argument(
            '--cr2.2',
            help='Set CellRanger2.2 algorithm for cell filter. If not designated, the EmptyDrops algorithm will be used.',
            default=False,
            dest="empty_cr_params",
            action='store_true'
        )
        scmd_rcell.add_argument(
            '--hard_filter', 
            help="Use a minium gene or UMI threshold to call microbes.", 
            default=None
        )
        scmd_rcell.add_argument(
            '--config',
            help="",
            default="NA",
            dest="config",
            type=str
        )
        scmd_rcell.add_argument(
            "--nosecondary", 
            help="Don't run the scanpy analysis.", 
            default=False,
            action='store_true'
        )
        scmd_rcell.add_argument(
            "--keep_bam", 
            help="Keep the bam file.", 
            default=False, 
            action='store_true'
        )
        scmd_rcell.add_argument(
            '--keep_unmap_reads', 
            help="Keep the unmapped reads.",
            default=False,
            action='store_true'
        )
        scmd_rcell.add_argument(
            '--multiplet_method', 
            help="The method to detect multiplet.", 
            choices=["scaled_softmax", "majority", "auto"], 
            type=str, 
            default="auto"
        )
        scmd_rcell.add_argument("--temperature", 
                        help="The Temperature of softmax.", 
                        default=2, 
                        type=int)
        scmd_rcell.add_argument(
            '--ksmiaa', 
            '--dev', 
            dest='ksmiaa',
            help="Enable development mode (same as --ksmiaa)", 
            default=False, 
            action='store_true'
        )
        scmd_rcell.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                                  help="Show this information.")
        scmd_mkindex._optionals.title = 'OPTIONS'
        args = parser.parse_args()
        return args

    def out_commandsLog(self, commandStr):
        outf = open(self.log_name, "w")
        outf.write(commandStr + "\n")
        outf.close()
        return commandStr

    def CustomHelpFormatter(self, prog):
        return argparse.HelpFormatter(prog, max_help_position=100, width=200)

def run_MobiVisionM():
    sys.path.append(os.path.dirname(__file__))
    m_args = ArgSetting("v1.3.2", "run_analysis_cmds.txt")
    cmd_args = m_args.get_all_args()
    ###check output dir
    if cmd_args.output_dir is None:
        if cmd_args.subcmds != "mkindex":
            print("An output path is required.")
            sys.exit()
        else:
            if cmd_args.genome == None:
                print("Name of speceis is required.")
                sys.exit()
            out_abs_path = os.path.join(os.getcwd(), '_and_'.join(cmd_args.genome))
    else:
        if os.path.isabs(cmd_args.output_dir):
            out_abs_path = cmd_args.output_dir
        else:
            out_abs_path = os.path.join(os.getcwd(), cmd_args.output_dir)         
    # check the if the path of output dir exists  -- cmd_args.output_dir  --> out_abs_path
    if os.path.exists(out_abs_path):
        mobilogger = MobiLoggingSystem(o_dir=out_abs_path, dev_mode=False)
        mobilogger._mobilogrecorder(log_message="Current working path is %s" %(os.getcwd()),
            log_level="INFO")
        mobilogger._mobilogrecorder(log_message="Using %s as output path." %(out_abs_path),
            log_level="INFO")  
        mobilogger._mobilogrecorder(log_message="The output path is already existed. Won't overwrite. %s" %(out_abs_path),
            log_level="ERROR")
        sys.exit()
    else:
        os.makedirs(out_abs_path, mode=0o755)
        mobilogger = MobiLoggingSystem(o_dir=out_abs_path, dev_mode=False)
        mobilogger._mobilogrecorder(log_message="Current working path is %s" %(os.getcwd()),
            log_level="INFO")
        mobilogger._mobilogrecorder(log_message="Using %s as output path." %(out_abs_path),
            log_level="INFO")  
    if cmd_args.subcmds == "quantify":
        from quantify import MobivisionProcessM
        if cmd_args.test_run:
            if os.path.exists(os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m")):
                cmd_args.input_dir = os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m", "demo")
                index_dir = os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m", "demo", "demo_ref")
                cmd_args.sample_ID = "demo"
            else:
                script_path = os.path.dirname(os.path.abspath(__file__))
                parent_path = os.path.dirname(script_path)
                cmd_args.input_dir = os.path.join(parent_path, "demo")
                index_dir = os.path.join(parent_path, "demo", "demo_ref")
                cmd_args.sample_ID = "demo"
        else:
            ## check fastq files and return real input data directory;
            if not os.path.exists(cmd_args.input_dir):
                mobilogger._mobilogrecorder(log_message="The path of input_dir is not exists. Plz recheck. %s" %(cmd_args.input_dir),
                    log_level="ERROR")
                sys.exit()
            ## index dir is the directory of reference index.
            index_dir = os.path.abspath(cmd_args.index_dir)
            if not os.path.exists(index_dir):
                mobilogger._mobilogrecorder(log_message="The path of index is not exists. Plz recheck. %s" %(index_dir), 
                    log_level="ERROR")
                sys.exit()
            elif not os.path.exists(os.path.join(index_dir, "star")) or \
                not os.path.exists(os.path.join(index_dir, "fasta")) or \
                not os.path.exists(os.path.join(index_dir, "genes")) or \
                not os.path.exists(os.path.join(index_dir, "reference.json")):
                mobilogger._mobilogrecorder(log_message="The path of index is not valid or not complete. Plz recheck or re-build. %s" %(index_dir), 
                    log_level="ERROR")
                sys.exit()
        os.chdir(out_abs_path)
        # check host
        if cmd_args.host_remove and not os.path.exists(cmd_args.host_reference):
            mobilogger._mobilogrecorder(log_message="The host reference is not exist. %s" %(cmd_args.host_reference), 
                log_level="ERROR")
            sys.exit()
        ## output the command line into the "run_analysis_cmds.txt"
        all_cmd = m_args.out_commandsLog( " ".join(["MobiVision-M", ' '.join(sys.argv[1:])]))
        try:
            core_num = int(cmd_args.core_num)
        except TypeError:
            mobilogger._mobilogrecorder(log_message="The threads must be a positive integer larger than 0 (larger than 8 recommanded). cmd_args.core_num" %(cmd_args.host_reference), 
                log_level="ERROR")
            sys.exit()
        if core_num <= 0:
            mobilogger._mobilogrecorder(log_message="The threads must be a positive integer larger than 0 (larger than 8 recommanded). cmd_args.core_num" %(cmd_args.host_reference), 
                log_level="ERROR")
            sys.exit()
        mp = MobivisionProcessM(fastq_path=cmd_args.input_dir, 
                            core_num=core_num, 
                            ref_path=index_dir, 
                            out_path=out_abs_path, 
                            vers_n=str(m_args.version_num), 
                            intron_incOpt="included", 
                            topcell_n=cmd_args.topc , 
                            emptydrop_opt=cmd_args.empty_cr_params, 
                            harder_filter=cmd_args.hard_filter, 
                            vers_kit=cmd_args.kitV, 
                            run_cmd=all_cmd, 
                            new_sampleid=cmd_args.sample_ID, 
                            with_CB=False, 
                            star_config=cmd_args.config, 
                            dev_mod=cmd_args.ksmiaa,
                            UMI_adjust=cmd_args.UMI_adjust, 
                            nosecondary=cmd_args.nosecondary, 
                            keep_bam=cmd_args.keep_bam, 
                            keep_unmapped=cmd_args.keep_unmap_reads, 
                            qc_only=cmd_args.qc_only, 
                            multiplet_method=cmd_args.multiplet_method, 
                            host_remove=cmd_args.host_remove, 
                            host_reference=cmd_args.host_reference, 
                            Temperature=cmd_args.temperature, 
                            mobilogger=mobilogger)
        mp.process()
        with open(os.path.join(out_abs_path, "Job_done.flag"), "w") as f:
            f.write("The quantify process done succeccfully.")
    ## mkdindex command initialize
    elif cmd_args.subcmds == "mkindex":
        from mkindex import IndexTool
        if cmd_args.test_run:
            in_genomes = ["demo"]
            if os.path.exists(os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m")):
                input_fasta_files = [os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m", "demo", "demo.fasta")]
                input_genes_files = [os.path.join(os.environ["CONDA_PREFIX"], "share", "mobivision-m", "demo", "demo.gtf")]
            else:
                script_path = os.path.dirname(os.path.abspath(__file__))
                parent_path = os.path.dirname(script_path)
                input_fasta_files = [os.path.join(parent_path, "demo", "demo.fasta")]
                input_genes_files = [os.path.join(parent_path, "demo", "demo.gtf")]
        else:
            if cmd_args.input_file != None:
                if os.path.exists(cmd_args.input_file):
                    import pandas as pd
                    mobilogger._mobilogrecorder(log_message="Using input file...", 
                        log_level="INFO")
                    i_df = pd.read_csv(cmd_args.input_file, sep="\t")
                    input_fasta_files = []
                    input_genes_files = []
                    in_genomes = []
                    for i in i_df.index:
                        input_fasta_files.append(i_df.loc[i, "fasta"])
                        input_genes_files.append(i_df.loc[i, "gtf"])
                        in_genomes.append(i_df.loc[i, "name"])
                else:
                    mobilogger._mobilogrecorder(log_message="Input file not found. %s" %(cmd_args.input_file), 
                        log_level="ERROR")
                    sys.exit()
            else:
                if cmd_args.fasta is None:
                    mobilogger._mobilogrecorder(log_message="No fasta file provided. Plz recheck.", 
                        log_level="ERROR")
                    sys.exit()
                if cmd_args.gtf is None:
                    mobilogger._mobilogrecorder(log_message="No gtf file provided. Plz recheck.", 
                        log_level="ERROR")
                    sys.exit()
                input_fasta_files, input_genes_files, in_genomes = cmd_args.fasta, cmd_args.gtf, cmd_args.genome
            tmp_dup = find_duplicates(in_genomes)
            if len(tmp_dup) > 0:
                mobilogger._mobilogrecorder(log_message="Duplicated name found in input genomes. Plz recheck. %s" %(tmp_dup), 
                    log_level="ERROR")
                sys.exit()
            tmp_dup = find_duplicates(input_fasta_files)
            none_exists = find_exist(input_fasta_files)
            if len(tmp_dup) > 0:
                mobilogger._mobilogrecorder(log_message="Duplicated name found in input fasta files. Plz recheck. %s" %(tmp_dup), 
                    log_level="ERROR")
                sys.exit()
            if len(none_exists) > 0:
                mobilogger._mobilogrecorder(log_message="Fasta files not found. Plz recheck. %s" %(none_exists), 
                    log_level="ERROR")
                sys.exit()
            tmp_dup = find_duplicates(input_genes_files)
            none_exists = find_exist(input_genes_files)
            if len(tmp_dup) > 0:
                mobilogger._mobilogrecorder(log_message="Duplicated name found in input gene fils. Plz recheck. %s" %(tmp_dup), 
                    log_level="ERROR")
                sys.exit()
            if len(none_exists) > 0:
                mobilogger._mobilogrecorder(log_message="Genes files not found. Plz recheck. %s" %(none_exists), 
                    log_level="ERROR")
                sys.exit()
        referenceUtil=IndexTool("STAR", in_genomes, input_fasta_files, input_genes_files, cmd_args.memgb, out_abs_path, mobilogger)
        referenceUtil.index_Builder(cmd_args.refv, "MobiVision-M v" + str(m_args.version_num))
        with open(os.path.join(out_abs_path, "Job_done.flag"), "w") as f:
            f.write("The quantify process done succeccfully.")
    elif cmd_args.subcmds == "rcmicrobe":
        from re_call_cell import ReCallCell
        from quantify import update_config, process_filter_arguments
        ###config
        default_config = update_config(star_config=cmd_args.config,
                                with_CB=False, 
                                species_number=1,  
                                UMI_adjust="no_adjust", 
                                mobilogger=mobilogger)
        ###apply filter val
        cell_filter_val, check_flag, check_info = process_filter_arguments(top_cell=cmd_args.topc,
                                                                    cr_flag=cmd_args.empty_cr_params, 
                                                                    hard_filter=cmd_args.hard_filter, 
                                                                    config=default_config)
        mobilogger._mobilogrecorder(log_message=check_info, 
            log_level=check_flag)
        if check_flag == "ERROR":
            sys.exit()
        ## check the output dir of quantify result and get absolute path of output dir
        if cmd_args.output_dir is None:
            out_abs_path = os.path.abspath("./")
        else:
            out_abs_path = os.path.abspath(cmd_args.output_dir)            
        os.chdir(out_abs_path)
        try:
            core_num = int(cmd_args.core_num)
        except TypeError:
            mobilogger._mobilogrecorder(log_message="The threads must be a positive integer larger than 0 (larger than 8 recommanded). cmd_args.core_num" %(cmd_args.host_reference), 
                log_level="ERROR")
            sys.exit()
        if core_num <= 0:
            mobilogger._mobilogrecorder(log_message="The threads must be a positive integer larger than 0 (larger than 8 recommanded). cmd_args.core_num" %(cmd_args.host_reference), 
                log_level="ERROR")
            sys.exit()
        p = ReCallCell(analysis_dir=cmd_args.analysis_dir, 
                        o_dir=out_abs_path, 
                        call_mtx=cmd_args.call_mtx, 
                        filter_args=cell_filter_val, 
                        sample_ID=cmd_args.sample_ID, 
                        ver_kit=cmd_args.kit, 
                        threads=core_num, 
                        nosecondary=cmd_args.nosecondary, 
                        keep_bam=cmd_args.keep_bam, 
                        keep_unmapped=cmd_args.keep_unmap_reads, 
                        multiplet_method=cmd_args.multiplet_method, 
                        config=default_config, 
                        mobilogger=mobilogger, 
                        Temperature=cmd_args.temperature)
        p.process()
        with open(os.path.join(out_abs_path, "Job_done.flag"), "w") as f:
            f.write("The quantify process done succeccfully.")
    return        

def CustomHelp(arglist):
    usage_info = {
                "quantify": '''Process Gene Expression data for High-throughput Single Microbe RNA-Seq kit
USAGE:   
MobiVision quantify [OPTIONS] -f <INPUT_DIR> -i <INDEX_PATH>


OPTIONS:
  -f, --fastqDir
                        The directory of fastq files. It should contain both fastq files of R1 and R2 with specific tag in file names such as "R1" for R1 fastq file.
                        Gzip fastq files are accepted. Check the user mannul for more information.
  -t, --threads
                        The number of threads.
  -i, --index_dir
                        Path of folder contains MobiVision-M compatible reference. Normaly made from sun-command mkindex.
  -o, --output_dir
                        The path of the output fold.
  -s, --sample_ID
                        The id of the result file user defined.
  --cellnumber          Force microbe number for microbe filter.
  --cr2.2               Set CellRanger2.2 algorithm for microbe filter. If not designated, the EmptyDrops algorithm will be used.
  --hard_filter         Use a hard threshold to call microbes. Such as "min_UMI:2000" or "min_reads:5000".
  --UMI_adjust          The method of UMI adjustion. Choose from "no_adjust", "step_1" and "step_1_and_2". Default is "no_adjust".
  --multiplet_method    The method to detect multiplet. Choose form "scaled_softmax", "majority" and "auto". Default is "auto".
  --nosecondary         Don't run the cluster analysis.
  --keep_bam            Keep the bam file of alignment.
  --keep_unmap_reads    Keep the un-aligned reads.
  --host_remove         Remove host reads from data.
  --host_reference      The reference of host. Made form sub-command mkindex.
  --qc_only             Only run the QC process.
  --test_run            Run the small demo data of quantify.
  --temperature         Temperature of softmax
  --kit                 Set the version of RNA kit used to construct the library. The default setting is 'v1.0'.
  --config              The config file with more detailed arguments. Check the user mannul for more information. Default is "NA"
  -h, --help            Show this information.''' ,
                "mkindex": '''Index preparation tool for MobiVision-M.

Build a MobiVision-M-compatible index folder with genome FASTA and gene GTF files of single or multiple species.
The command "mkindex" should be preceded by "MobiVision-M"

USAGE:
MobiVision-M mkindex [OPTIONS] -n <GENOME> -f <FASTA_PATH> -g <GTF_PATH>

OPTIONS:
  -n, --nameOfSpecies
                        Unique genome name [a-zA-Z0-9_]+. If the index is created with multiple genomes, Plz specify the -n
                        argument multiple times. For example, if one specifies -n <genome1> and another -n <genome2>, the output
                        index folder will be named <genome1>_and_<genome2> if -o is not assigned.

  -f, --fastaPath
                        Genome fasta file path. If the index is created with multiple FASTA files, Plz specify the -f argument
                        multiple times. The supported formats of genome fasta files include '.fasta', '.fa.gz' and '.fna.gz'.
  -g, --gtfPath
                        GTF file path for genome. If the index is created with multiple GTFs, Plz specify the -g argument
                        multiple times. If -n and -f arguments were specified multiple times, -g argument needs to be specified 
                        an equal amount of times and the specified order must be same.              
  --inut_file
                        he Tab-separated txt file contains the list of fasta and gtf file. Default value is NA, which means no 
                        file provided. Refer to the user manual for detailed information.
  -r, --referenceVerString
                        The version of output genome index, e.g. "ref-v1.0", "ref-2022-06".
  -m, --memoryUsed
                        Maximum memory (GB) used when building index files with STAR.
  -o, --output_dir
                        The name of output fold. if not assigned, it will be will be named as <genome1>_and_<genome2>.
  --test_run            Run a small demo data of mkindex.
  -h, --help            Show this information.''',
                "rcmicrobe": '''Re-call microbes form a MobiVision-M result.

Performance a different microbe-calling method on a finished MobiVision-M result without pre-process and alignment.
The command "rcmicrobe" should be preceded by "MobiVision-M"

USAGE:
MobiVision-M rcmicrobe [OPTIONS] -i <Analysised Results> -o <Output Dir>

OPTIONS:
  -i, --analysis_dir 
                        The fold path of analysised MobiVision-M results.
  -o, --output_dir
                        The fold path of output results. Won't overwrite an existed path.
  -c, --call_mtx        Re-call microbes form which matrix. Default is raw_re-assigned_cell_gene_matrix.
  -t, --threads 
                        The number of threads. 
  -s, --sample_ID
                        The id of the result file user defined. 
  --cellnumber          Force microbe number for microbe filter.
  --cr2.2               Set CellRanger2.2 algorithm for microbe filter. If not designated, the EmptyDrops algorithm will be used.
  --hard_filter         Use a hard threshold to call microbes. Such as "min_UMI:2000" or "min_reads:5000".
  --multiplet_method    The method to detect multiplet. Choose form "scaled_softmax", "majority" and "auto". Default is "auto".
  --nosecondary         Don't run the cluster analysis.
  --keep_bam            Keep the bam file of alignment.
  --keep_unmap_reads       Keep the un-aligned reads.
  --kit                 Set the version of RNA kit used to construct the library. The default setting is 'v1.0'.
  --temperature         Temperature of softmax
  --config              The config file with more detailed arguments. Check mannul for more information. Default is "NA"
  -h, --help            Show this information.''', 
                 "MobiVision_help": '''MobiVision-M-v1.3.2
Process Single Microbe Sequencing data


USAGE:
    MobiVision-M <SUBCOMMAND>


FLAGS:
    -h, --help     Print help information
    -V, --version  Print version information

SUBCOMMANDS:
      
    quantify      Count gene expression reads from a single sample.
    mkindex       Prepare a reference for MobiVision-M software'
    rcmicrobe     Re-call microbes form an analysised MobiVision-M result.'''          

}
    if len(arglist) == 2:
        if arglist[1] == "-h" or arglist[1] == "--help":
            #print help information for all subcommands of MobiVision-M
            print(usage_info["MobiVision_help"])
            sys.exit(0)
        elif arglist[1] in ["quantify", "mkindex", "rcmicrobe"]:
            #Prints help information or the help of the given subcommand(s).
            print(usage_info[arglist[1]])
            sys.exit(0)
        elif arglist[1] == "-V" or arglist[1] == "-v" or arglist[1] == "--version":
            print("MobiVision-M version v1.3.2")
            sys.exit(0)
        else:
            print("invalid choice: '{}' (choose from 'quantify', 'mkindex')".format(arglist[1]))
            sys.exit(0)            
    else:
        if len(arglist) == 1 and os.path.basename(arglist[0]) == "MobiVision-M":
            #print help information for all subcommands of MobiVision-M
            print(usage_info["MobiVision_help"])
            sys.exit(0)
        elif len(arglist) == 3:
            if arglist[1] in ["quantify", "mkindex", "rcmicrobe"] and arglist[2] in ["-h", "--help"]:
                print(usage_info[arglist[1]])
                sys.exit(0)
            else:
                print("invalid choice: '{}' (choose from 'quantify', 'mkindex')".format(arglist[1]))
                sys.exit(0)

    


