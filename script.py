import argparse
import subprocess 
import yaml
import os  

parser = argparse.ArgumentParser(description='Process sequencing data with Snakemake.')  

parser.add_argument('--type', choices=['se', 'pe'], help='Input type: single-end (se) or paired-end (pe)')  
parser.add_argument('--threads', type=int, default=2, help='Number of threads to use for processing, default:2')  
parser.add_argument('--remove-intermediate_output', action='store_true', help='Enable temporary file handling(default: False)')  
  
parser.add_argument('--hg38-ref', type=str, default="Database/hg38/hg38.fa", help='Path to the hg38 reference genome FASTA file, default:"Database/hg38/hg38.fa"')  
parser.add_argument('--bwa-ref', type=str, default="Database/hg38/bwa_index/hg38.fa", help='Path to the BWA-indexed hg38 reference genome, default:"Database/hg38/bowtie2_build/hg38"')  
parser.add_argument('--bowtie2-ref', type=str, default="Database/hg38/bowtie2_build/hg38", help='Path to the Bowtie2-indexed hg38 reference genome, default:"Database/hg38/bowtie2_build/hg38"')  
parser.add_argument('--kraken2-ref', type=str, default="Database/k2_pluspf_20240605", help='Path to the Kraken2 database, default:"Database/k2_pluspf_20240605"') 

parser.add_argument('--confidence', type=float,default="0.0",help=' Confidence score threshold (default: 0.0) when using Kraken2; must be in [0, 1].')
VERSION = "1.0.0"
parser.add_argument('--version',action='store_true',help="show program's version number and exit")  
group = parser.add_mutually_exclusive_group(required=False)  
group.add_argument('-i0','--input0', type=str,help='Single input file for single-end or first input file for paired-end (if --type=pe and single file mode)')  
group.add_argument('-i1','--input1', type=str,help='Read1 input file for paired-end (if --type=pe)')  
parser.add_argument('-i2', '--input2',type=str, help='Read2 input file for paired-end (if --type=pe)')  
group.add_argument('-l','--samples-list', type=argparse.FileType('r'),help='File containing sample names of input files (one per line)')  
parser.add_argument('--input-directory', type=str, help='Directory containing input files (if using --list)')  
parser.add_argument('--output-directory', type=str,default="snakemake-out",help='Output directory,default="snakemake-out"')
parser.add_argument('--ext1', type=str,default=".R1.fastq.gz",help='Suffix of the R1 input filename(if using --list),default=".R1.fastq.gz"',)
parser.add_argument('--ext2', type=str,default=".R2.fastq.gz",help='Suffix of the R2 input filename(if using --list),default=".R2.fastq.gz"')

parser.add_argument('--cores', type=int, default=1, help="Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the number of available CPU cores. In case of cluster/cloud execution, this argument sets the maximum number of cores requested from the cluster or cloud scheduler. (See https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources-remote-execution for more info)This number is available to rules via workflow.cores.")
parser.add_argument('--configfile', type=str,default="script.yaml" ,help="Specify or overwrite the config file of the workflow (see the docs). Values specified in JSON or YAML format are available in the global config dictionary inside the workflow. default:'script.yaml'")
parser.add_argument('--jobs', type=int, help="Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores (it is though recommended to use --cores in that case). Note: Set to 'unlimited' to allow any number of parallel jobs.")
parser.add_argument('--keep-going', action='store_true', help='Continue with other jobs even if some fail(default: False)')
parser.add_argument('-n','--dry-run', action='store_true', help='Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs. (default: False)')
parser.add_argument('--verbose', action='store_true', help='Increase verbosity of output(default: False)')
parser.add_argument('--workflow-profile', type=str, help="ath (relative to current directory) to workflow specific profile folder to use for configuring Snakemake with parameters specific for this workflow (like resources). If this flag is not used, Snakemake will by default use 'profiles/default' if present (searched both relative to current directory and relative to Snakefile, in this order). For skipping any workflow specific profile provide the special value 'none'. Settings made in the workflow profile will override settings made in the general profile (see --profile). The profile folder has to contain a file 'config.yaml'. This file can be used to set default values for command line options in YAML format. For example, '--executor slurm' becomes 'executor: slurm' in the YAML file.It is advisable to use the workflow profile to set or overwrite e.g. workflow specific resources like the amount of threads of a particular rule or the amount of memory needed. Note that in such cases, the arguments may be given as nested YAML mappings in the profile, e.g. 'set-threads: myrule: 4' instead of 'set-threads: myrule=4'.")  
parser.add_argument('--list-rules',  action='store_true', help='Show available rules in given Snakefile. (default: False)')  
parser.add_argument('--dag', action='store_true', help='Print the directed acyclic graph of jobs in the dot language.')  
parser.add_argument('--summary', '-S', action='store_true', help='Print a summary of all files created by the workflow.')              
parser.add_argument('--nolock', action='store_true', help='Do not lock the working directory (default: False)')  
               
args = parser.parse_args()
  
config = {
  'database': {  
    'hg38_ref': args.hg38_ref,
    'bwa_ref': args.bwa_ref,
    'bowtie2_ref': args.bowtie2_ref,
    'kraken2_ref': args.kraken2_ref
   },  
  'output_directory': {
    'bwa_outdir': "/s1_bwa_samtools",
    'fastp_bwa_outdir': '/s2_fastp_bwa',
    'bbduk_outdir': '/s3_bbduk',
    'fastp_bbduk_outdir': '/s4_fastp_bbduk',
    'kneaddata_outdir': '/s5_kneaddata',
    'fastp_kneaddata_outdir': '/s6_fastp_kneaddata',
    'taxa_annotation_outdir': '/s7_taxa_annotation',
    'fp_report_outdir': '/fp_report'
    },
  'parameters': {
    'temp': args.remove_intermediate_output,
    'threads': args.threads,
    'confidence': args.confidence
  },  
  'samples': {
    'ext1': args.ext1,
    'ext2': args.ext2,
    'input0': args.input0 if args.type == 'se' and args.input0 else None,
    'input1': args.input1 if args.type == 'pe' and args.input1 else None,
    'input2': args.input2 if args.type == 'pe' and args.input2 else None,
    'samples_list': args.samples_list.name if args.samples_list else None,
    'samples_indir': args.input_directory if args.input_directory else None,
    'outdir': args.output_directory
  }
  }
with open(args.configfile, "w") as f:
  yaml.dump(config, f,default_flow_style=False)

snakefile = None
if args.type == 'pe' and args.samples_list:
  snakefile = 'pe-list.py'
if args.type == 'pe' and args.input1:
  snakefile = 'pe-i.py'
if args.type == 'se' and args.samples_list:
  snakefile = 'se-list.py'
if args.type == 'se' and args.input0:
  snakefile = 'se-i.py'
  
if args.version:  
  print(f"{VERSION}")  
  exit(0)  

if args.workflow_profile== 'none':  
  args.workflow_profile = None  

snakemake_cmd = ['snakemake', '-s', snakefile]
if args.jobs:
  snakemake_cmd += ['--jobs', str(args.jobs)]
if args.keep_going:
  snakemake_cmd.append('--keep-going')
if args.dry_run:
  snakemake_cmd.append('--dry-run')
if args.verbose:
  snakemake_cmd.append('--verbose')
if args.configfile:
  snakemake_cmd += ['--configfile', args.configfile]  
if args.workflow_profile:
  snakemake_cmd += ['--workflow-profile', args.workflow_profile]    
if args.list_rules:
  snakemake_cmd.append('--list-rules')
if args.dag:
  snakemake_cmd.append('--dag')
if args.summary:
  snakemake_cmd.append('--summary')
if args.nolock:
  snakemake_cmd.append('--nolock')
  
subprocess.run(snakemake_cmd)