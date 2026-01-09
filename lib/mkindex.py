import hashlib
import csv
import os, sys
import json
import re
import math
import collections
import stat
sys.path.append(os.path.dirname(__file__))
from cmdutil import run_command_safely
from mobivisionlogging import MobiLoggingSystem, MobiCommandLogSystem
from mobivisionexecutor import CommandExecutor

class SeqFile:
    """ a class for fasta file process,  which will be used in indexbuilder, ... .
    """
    def __init__(self, filenames = []):
        self.fn = filenames
        
    def parse_fasta(self):
        if len(self.fn) < 1:
            return "No seq files parsed, please check the input parameters."
        
        num_chrs = 0
        
        first_fafile = self.fn[0]
        fh = open( first_fafile, 'r')
        #sequences = {}
        for line in fh:
            if line.startswith(">"):
                name = line.rstrip("\n")
                #sequences[name] = ""
                num_chrs += 1
            else:
                continue
                #sequences[name] = sequences[name] + line.rstrip("\n")
        #return ((i[0]+"\n"+i[1]) for i in sequences.items())
        return num_chrs
    
    def add_chr_prefix(self, chr_prefixes, out_faname):
        """
        add prefix for chromosomes in fasta file of each species 
        """
        num_chrs = 0
        with open(out_faname, 'w') as f:
            for genome_prefix, in_fasta_fn in zip(chr_prefixes, self.fn):
                with open(in_fasta_fn, 'r') as g:
                    for line in g:
                        line = line.strip()
                        if line.startswith('>'):
                            line = '>' + genome_prefix + '_' + line[1:]
                            num_chrs += 1
                        f.write(line + '\n')
        return num_chrs
        


class IndexTool:
    """
    mkindex subcommand utility class for mobivision
    """

    def __init__(self, star_bin_path, genomes, in_fasta_fns, in_gtf_fns, memgb=None, output_name="NA", mobilogger=None):

        self.align_tool = star_bin_path

        self.mem_gb = memgb
        self.genomes = genomes
        self.in_fasta_fns = in_fasta_fns
        self.in_gtf_fns = in_gtf_fns
        if output_name == "NA":
            self.out_dir = self.get_outdir_from_names()
        else:
            self.out_dir = output_name
        if mobilogger == None:
            if os.path.exists(output_name):
                self.mobilogger = MobiLoggingSystem(o_dir=output_name, dev_mode=False)
                self.mobilogger._mobilogrecorder(log_message="The output path is already existed. Won't overwrite.",
                    log_level="ERROR")
                sys.exit()
            else:
                os.makedirs(output_name, mode=0o755)
                self.mobilogger = MobiLoggingSystem(o_dir=output_name, dev_mode=False)
        else:
            self.mobilogger = mobilogger
        self.mobicommandlogger = MobiCommandLogSystem(o_dir=self.mobilogger.working_path, dev_mode=False)
        self.mobiexecutor = CommandExecutor(log_system=self.mobicommandlogger, console_output=False)
        self.format_genome_prefnames()

    def index_Builder(self, ref_version, mobi_version, mobilogger=None, num_threads=1):

        """
        initialize star index builder, and create merged gtf
        """

        #print("Creating new reference folder at %s" % self.out_dir)
        self.mobilogger._mobilogrecorder(log_message="Creating new reference folder at %s" % self.out_dir,
            log_level="INFO")
        #print("Writing genome fasta and genes GTF file into reference folder...")
        self.mobilogger._mobilogrecorder(log_message="Writing genome fasta and genes GTF file into reference folder...",
            log_level="INFO")
        fafile = os.path.join(self.out_dir, "fasta/genome.fa")
        os.mkdir(os.path.dirname(fafile))

        gtffile = os.path.join(self.out_dir, "genes/genes.gtf")
        os.mkdir(os.path.dirname(gtffile))
        genefile = os.path.join(self.out_dir, "genes/gene_info.json")
        
        ## check if the fasta files is gzipped, if yes, then uncompress it
        for file_n, f in enumerate(self.in_fasta_fns):
            if f.endswith(".gz"):
                ##run_command_safely('gzip', ['-dc',  out_genomename])
                new_uncomp_file = os.path.join(self.out_dir, os.path.basename(f)[:-3] )
                cmd = 'gzip -dc %s > %s'%(f, new_uncomp_file)
                exit_code = self.mobiexecutor.execute(
                    command=cmd,
                    context={
                        "mkindex": "gzip"
                    }, console_output=False
                )
                if exit_code != 0:
                    self.mobilogger._mobilogrecorder(log_message="Gzip failed on file: %s" %(f),
                        log_level="ERROR")
                self.in_fasta_fns[file_n] = new_uncomp_file
                
        ## check if the gtf files is gzipped, if yes, then uncompress it
        for file_n, f in enumerate(self.in_gtf_fns):
            if f.endswith(".gz"):
                new_uncomp_file = os.path.join(self.out_dir, os.path.basename(f)[:-3] )
                cmd = 'gzip -dc %s > %s'%(f, new_uncomp_file)
                exit_code = self.mobiexecutor.execute(
                    command=cmd,
                    context={
                        "mkindex": "gzip"
                    }, console_output=False
                )
                if exit_code != 0:
                    self.mobilogger._mobilogrecorder(log_message="Gzip failed on file: %s" %(f),
                        log_level="ERROR")
                self.in_gtf_fns[file_n] = new_uncomp_file
        
        ## create result genome fasta file and genes gtf file.
        genome_num_chrs = self.merge_genome_infofiles(fafile, gtffile, genefile, self.mobilogger)
        #print("...done\n")
        self.mobilogger._mobilogrecorder(log_message="Writing genome fasta and genes GTF file into reference folder is done successfully.",
            log_level="INFO")
        
        ## remove the temporary uncompressed file in output directory.
        for f in os.listdir(self.out_dir):
            if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fna') or f.endswith('.gtf'):
                uncompress_file_path = os.path.join(self.out_dir, f)
                os.remove(uncompress_file_path)

        #print("Computing hash of genome FASTA file and GTF file...")
        self.mobilogger._mobilogrecorder(log_message="Computing hash of genome FASTA file and GTF file...",
            log_level="INFO")
        fasta_hash = self.compute_hash(fafile)
        gtf_hash = self.compute_hash(gtffile)
        self.mobilogger._mobilogrecorder(log_message="Computing hash of genome FASTA file and GTF file is done successfully.",
            log_level="INFO")
        #print("...done\n")

        jsonInfo_data = {
            "fasta_hash": fasta_hash,
            "genomes": self.genomes,
            "gtf_hash": gtf_hash,
            "input_fasta_files": [os.path.basename(x) for x in self.in_fasta_fns],
            "input_gtf_files": [os.path.basename(x) for x in self.in_gtf_fns],
            "mem_gb": self.mem_gb,
            "mobi_version": mobi_version,  #### "mkindex_v1.0"  - buildtools_version,
            "threads": num_threads,
            "ref_version": ref_version
        }

        #print("Writing genome metadata JSON file into reference folder...")
        self.mobilogger._mobilogrecorder(log_message="Writing genome metadata JSON file into reference folder...",
            log_level="INFO")
        metadata_json = os.path.join(self.out_dir, "reference.json")
        with open(metadata_json, 'w') as f:
            json.dump(jsonInfo_data, f, sort_keys=True, indent=4)
        #print("...done\n")
        self.mobilogger._mobilogrecorder(log_message="Writing genome metadata JSON file into reference folder is done successfully.",
            log_level="INFO")

        #print("Indexing genome FASTA file...")
        self.mobilogger._mobilogrecorder(log_message="Indexing genome FASTA file...", 
            log_level="INFO")
        #subprocess.check_call(["samtools", "faidx", fafile])
        cmd = "samtools faidx %s" %(fafile)
        exit_code = self.mobiexecutor.execute(
            command=cmd,
            context={
                "mkindex": "samtools faidx"
            }, console_output=False
        )
        if exit_code != 0:
            self.mobilogger._mobilogrecorder(log_message="Indexing failed on file: %s" %(fafile),
                log_level="ERROR")
        #print("...done\n")
        self.mobilogger._mobilogrecorder(log_message="Indexing genome FASTA file is done successfully.",
            log_level="INFO")

        genome_size_b = float(os.path.getsize(fafile))
        genome_size_gb = float(genome_size_b) / float(10 ** 9)

        sa_index_n_bases = min(14, int(math.log(genome_size_b, 2) / 2 - 1))
        chr_bin_n_bits = min(18, int(math.log(genome_size_b / genome_num_chrs, 2)))

        if self.mem_gb is None:
            sa_sparse_d = None
            limit_ram = None
        else:
            # Total memory = SA memory + SA index memory + Genome memory
            # SA memory = (8 * genome size) / genomeSAsparseD
            # SA index memory = 8*(4**genomeSAindexNbases)
            # Genome memory = genome size
            sa_index_mem_gb = float(8 * (4 ** sa_index_n_bases)) / float(10 ** 9)
            genome_mem_gb = genome_size_gb

            min_mem_gb = int(genome_mem_gb + sa_index_mem_gb + 3)
            if self.mem_gb < min_mem_gb:
                #sys.exit(
                #    "WARNING: STAR requires at least %d GB of memory when aligning reads to your reference.\nPlease start again with --memgb=%d." % (
                #        min_mem_gb, min_mem_gb))
                self.mobilogger._mobilogrecorder(log_message="STAR requires at least %d GB of memory when aligning reads to your reference. Please start again with --memgb=%d." %(min_mem_gb, min_mem_gb),
                    log_level="ERROR")
                sys.exit()

            limit_ram = self.mem_gb * 1024 ** 3

            # 2 GB of buffer space
            self.mem_gb = max(1, self.mem_gb - 2)

            sa_sparse_d = float(8 * genome_size_gb) / (float(self.mem_gb) - genome_mem_gb - sa_index_mem_gb)
            sa_sparse_d = max(1, int(math.ceil(sa_sparse_d)))

        args = [self.align_tool, '--runMode', 'genomeGenerate', '--genomeDir', os.path.join(self.out_dir, "star"),
                '--runThreadN', str(num_threads), '--genomeFastaFiles', fafile,
                '--sjdbGTFfile', gtffile]
        if limit_ram is not None:
            args += ['--limitGenomeGenerateRAM', str(limit_ram)]
        if sa_sparse_d is not None:
            args += ['--genomeSAsparseD', str(sa_sparse_d)]
        if sa_index_n_bases is not None:
            args += ['--genomeSAindexNbases', str(sa_index_n_bases)]
        if chr_bin_n_bits is not None:
            args += ['--genomeChrBinNbits', str(chr_bin_n_bits)]
        
        self.mobilogger._mobilogrecorder(log_message="Building index...", 
            log_level="INFO")
        #subprocess.check_call(args, start_new_session=True)
        cmd = " ".join(args)
        exit_code = self.mobiexecutor.execute(
            command=cmd,
            context={
                "mkindex": "STAR index"
            }, console_output=False
        )
        if exit_code != 0:
            self.mobilogger._mobilogrecorder(log_message="Building index failed.",
                log_level="ERROR")
            sys.exit()

        os.chmod(self.out_dir , stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH + stat.S_IXOTH)
        for root, dirs, files in os.walk(self.out_dir):
            for dir in dirs:
                dp = os.path.join(root, dir)
                os.chmod(dp , stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH + stat.S_IXOTH)
        self.mobilogger._mobilogrecorder(log_message="Building index is done successfully.", 
            log_level="INFO")


    def format_genome_prefnames(self):

        """ format output name for input genomes names """

        if len(self.genomes) > 1:
            max_length = max([len(g) for g in self.genomes])
            self.genome_prefixes = []
            for genome in self.genomes:
                genome_prefix = genome
                if len(genome_prefix) < max_length:
                    genome_prefix += '_' * (max_length - len(genome_prefix))
                assert genome_prefix not in self.genome_prefixes
                self.genome_prefixes.append(genome_prefix)
        else:
            self.genome_prefixes = self.genomes

    def merge_genome_infofiles(self, out_genomename, out_gtfname, out_gene, mobilogger):

        """ merge input files ,and output name is set with out_genomename and out_gtfname
        """

        # write fasta file
        if len(self.genomes) > 1:
            files_handle = SeqFile(self.in_fasta_fns)
            n_chrs = files_handle.add_chr_prefix(self.genome_prefixes, out_genomename)
        else:
            run_command_safely('cp', [self.in_fasta_fns[0], out_genomename])
            files_handle = SeqFile([out_genomename])
            n_chrs = files_handle.parse_fasta()
        # write gtf file
        gene_dict = {}
        with open(out_gtfname, 'w') as f:
            writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            for genome_prefix, in_gtf_fn in zip(self.genome_prefixes, self.in_gtf_fns):
                if len(self.genomes) > 1:
                    prefix_func = lambda s: '%s_%s' % (genome_prefix, s)
                else:
                    prefix_func = lambda s: s
                transcript_to_chrom = {}
                cross_chrom_transcripts = set()
                with open(in_gtf_fn, 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    for i, row in enumerate(reader):
                        add_key = "NAZZZ"
                        add_info = {"species":genome_prefix, "type":"unknown"}
                        if len(row) == 0:
                            continue
                        if row[0].startswith('#'):
                            writer.writerow(row)
                            continue
                        chrom = prefix_func(row[0])
                        row[0] = chrom
                        properties = collections.OrderedDict()
                        pattern = re.compile('(\S+?)\s*"(.*?)"')
                        for m in re.finditer(pattern, row[8]):
                            key = m.group(1)
                            value = m.group(2)
                            properties[key] = value
                        ## row[2] is annotation
                        if row[2] == 'exon' and (
                                'transcript_id' not in properties or ';' in properties['transcript_id'] or re.search(
                            r'\s', properties['transcript_id']) is not None or \
                                'gene_id' not in properties or ';' in properties['gene_id']):
                            sys.exit(
                                "Invalid character in attributes columns in GTF line %d: %s\n\nYou could use \"sed -n '%dp' %s\"\n%s" % (
                                    i + 1, in_gtf_fn, i + 1, in_gtf_fn, "Please fix your GTF and start again."))
                        if 'transcript_id' in properties and properties['transcript_id'] != "":
                            properties['transcript_id'] = prefix_func(properties['transcript_id'])
                            curr_tx = properties['transcript_id']
                            if curr_tx in transcript_to_chrom and transcript_to_chrom[
                                curr_tx] != chrom and curr_tx != "":
                                # ignore recurrences of a transcript on different chromosomes - it will break the STAR index
                                cross_chrom_transcripts.add(curr_tx)
                                continue
                            transcript_to_chrom[curr_tx] = chrom
                        if 'gene_id' in properties:
                            properties['gene_id'] = prefix_func(properties['gene_id'])
                        if 'gene_name' in properties:
                            properties['gene_name'] = prefix_func(properties['gene_name'])
                            add_key = properties['gene_name']
                        if 'gene_biotype' in properties:
                            add_info["type"] = properties['gene_biotype']
                        if 'transcript_biotype' in properties:
                            add_info["type"] = properties['transcript_biotype']
                        if add_key != "NAZZZ":
                            gene_dict[add_key] = add_info
                        properties_str = []
                        for key, value in properties.items():
                            properties_str.append('%s "%s"' % (key, value))
                        row[8] = '; '.join(properties_str)
                        writer.writerow(row)
                if len(cross_chrom_transcripts) > 0:
                    #print("WARNING: The following transcripts appear on multiple chromosomes in the GTF:")
                    #print('\n'.join(list(cross_chrom_transcripts)) + '\n')
                    #print(
                    #    "This can indicate a problem with the reference or annotations. Only the first chromosome will be counted.")
                    mobilogger._mobilogrecorder(log_message="The following transcripts appear on multiple chromosomes in the GTF:", 
                        log_level="WARNING")
                    mobilogger._mobilogrecorder(log_message='\n'.join(list(cross_chrom_transcripts)), 
                        log_level="WARNING")
                    mobilogger._mobilogrecorder(log_message="This can indicate a problem with the reference or annotations. Only the first chromosome will be counted.", 
                        log_level="WARNING")
        if len(self.genomes) >= 1:
            with open(out_gene, "w") as json_file:
                json.dump(gene_dict, json_file, indent=4)
        return n_chrs

    def get_outdir_from_names(self):
        return '_and_'.join(self.genomes)

    #    @staticmethod
    def compute_hash(self, filename, block_size_bytes=2 ** 20):

        """ compute hash value of file """

        digest = hashlib.sha1()
        with open(filename, 'rb') as f:
            for chunk in iter(lambda: f.read(block_size_bytes), b''):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def write_empty_json(filename):
        with open(filename, 'w') as f:
            json.dump({}, f)

