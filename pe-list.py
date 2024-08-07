import pandas as pd
import csv
import os
import json

configfile: "script.yaml"
workdir: os.getcwd()

with open(config['samples']["samples_list"], 'r') as file:
  reader = csv.reader(file)  
  Sample = []  
  for row in reader:
    Sample.append([row[i] for i in range(len(row))])
Samples = [''.join(code) for code in Sample]
print(Samples)

Threads= config['parameters']["threads"]

Outdir1= config['samples']['outdir']+config['output_directory']["bwa_outdir"]
Outdir2= Outdir1+config['output_directory']["fastp_bwa_outdir"]
Outdir3= config['samples']['outdir']+config['output_directory']["bbduk_outdir"]
Outdir4= Outdir3+config['output_directory']["fastp_bbduk_outdir"]
Outdir5= config['samples']['outdir']+config['output_directory']["kneaddata_outdir"]
Outdir6= Outdir5+config['output_directory']["fastp_kneaddata_outdir"]
Outdir7= config['samples']['outdir']+config['output_directory']["taxa_annotation_outdir"]
Outdir8= config['samples']['outdir']+config['output_directory']["fp_report_outdir"]

wildcard_constraints:  
    sample=r"\w+" 

rule all:
  input:
    Outdir8+"/coverage.txt",
    Outdir8+"/fp_bwa.txt",
    Outdir8+"/fp_bbduk.txt",
    Outdir8+"/fp_kneaddata.txt",
    expand(Outdir7+"/{sample}_taxa.count",sample=Samples),
    expand(Outdir7+"/{sample}_taxa.percent",sample=Samples),
    Outdir7+"/all_count_result",
    os.path.dirname(config['database']['bwa_ref']),
    os.path.dirname(config['database']['bowtie2_ref'])
    
    
####step0 : arrangement
rule bwa_index:
  input:
    config['database']["hg38_ref"]
  output:
    directory(os.path.dirname(config['database']['bwa_ref']))
  params:
    os.path.basename(config['database']['bwa_ref'])
  shell:
    "bwa index {input} {output}/{params}"

rule bowtie2_build:  
  input:  
    config['database']["hg38_ref"]  
  output:  
    directory(os.path.dirname(config['database']['bowtie2_ref']) ) 
  params:
    os.path.basename(config['database']['bowtie2_ref'])
  shell:  
    "bowtie2-build -f {input} {output}/{params}" 
    
#### step1 : alignment + select unmap    
rule select_unmap_ppmi:
  input:
    r1=config['samples']["samples_indir"]+"/{sample}/{sample}"+config['samples']["ext1"],
    r2=config['samples']["samples_indir"]+"/{sample}/{sample}"+config['samples']["ext2"]
  output:
    Outdir1+"/alignment/{sample}.bam"
  params:  
    bwa_ref=config['database']['bwa_ref']
  shell:
    "bwa mem -k 19 -K 10000000 -t {Threads} {params.bwa_ref} {input.r1} {input.r2}  > {output}.unsorted.sam && "
    "samtools sort --threads {Threads} {output}.unsorted.sam -o {output} && "
    "rm {output}.unsorted.sam"

rule samtools:
  input:
    Outdir1+"/alignment/{sample}.bam"
  output:
    Outdir1+"/{sample}.R1.fastq.gz",
    Outdir1+"/{sample}.R2.fastq.gz",
    Outdir8+"/{sample}.coverage.txt",
    temp(Outdir8+"/{sample}.depth.txt")
  run:
    shell("""
      samtools fastq -f 12 {input}  -1 {output[0]} -2 {output[1]} &&
      samtools depth -a {input} > {output[3]}
      """)
    def calculate_average_depth(depth_file):  
      total_depth = 0  
      num_positions = 0  
      with open(str(depth_file), 'r') as f:  
        for line in f:  
          chrom, pos, depth = line.strip().split('\t')  
          depth = int(depth)  
          if depth > 0: 
            total_depth += depth  
            num_positions += 1  
      if num_positions > 0:  
        return total_depth / num_positions  
      else:  
        return 0  
    average_depth = calculate_average_depth(output[3])
    with open(output[2], 'w') as out:
      out.write(f"Sample: {wildcards.sample}\tAverage depth: {average_depth}")

rule coverage:
  input:
    expand(Outdir8+"/{sample}.coverage.txt",sample=Samples)
  output:
    Outdir8+"/coverage.txt"
  run:
    with open(output[0], 'w') as outfile:  
      for fname in input:  
        with open(fname, 'r') as infile:  
          outfile.write(infile.read() + '\n')
    shell("rm {input}")  
 
#### step2 : fastp_1bwa   
rule fastp_1bwa:
  input:
    Outdir1+"/{sample}.R1.fastq.gz",
    Outdir1+"/{sample}.R2.fastq.gz"
  output:
    Outdir2+"/{sample}/out.R1.fastq.gz",
    Outdir2+"/{sample}/out.R2.fastq.gz",
    Outdir2+"/{sample}/fastp_out.html",
    Outdir2+"/{sample}/fastp_out.json"
  shell: 
    "fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -h {output[2]} -j {output[3]} "
   
rule fastp_bwa_stats:
  input:
    expand(Outdir2+"/{sample}/fastp_out.json",sample=Samples)
  output:
    Outdir8+"/fp_bwa.txt"
  run:
    with open(output[0], 'w') as out:
      out.write("samples\ttotal_bases\treads\tQ20_rate\tQ30_rate\tGC_content\n")
      for json_file in input:
        with open(str(json_file)) as f:
          data = json.load(f)
          samples = os.path.basename(os.path.dirname(os.path.abspath(json_file))) 
          size = data["summary"]["after_filtering"]["total_bases"]
          reads = data["summary"]["after_filtering"]["total_reads"]
          Q20 = data["summary"]["after_filtering"]["q20_rate"]
          Q30 = data["summary"]["after_filtering"]["q30_rate"]
          GC = data["summary"]["after_filtering"]["gc_content"]
          out.write(f"{samples}\t{size}\t{reads}\t{Q20}\t{Q30}\t{GC}\n")
 
#### step3 : quality control
rule quality_filter:
  input:
    Outdir2+"/{sample}/out.R1.fastq.gz",
    Outdir2+"/{sample}/out.R2.fastq.gz"
  output:
    Outdir3+"/{sample}_qc1.R1.fastq.gz",
    Outdir3+"/{sample}_qc1.R2.fastq.gz"
  shell:
    "bbduk.sh in1={input[0]} in2={input[1]} "
    "out1={output[0]} out2={output[1]} "
    "threads={Threads} "
    "maq=10 "
    "overwrite=t"
    
rule quality_trim:
  input:
    Outdir3+"/{sample}_qc1.R1.fastq.gz",
    Outdir3+"/{sample}_qc1.R2.fastq.gz"
  output:
    Outdir3+"/{sample}_qc2.R1.fastq.gz",
    Outdir3+"/{sample}_qc2.R2.fastq.gz"
  shell: "bbduk.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} threads={Threads} qtrim=r1 trimq=10 overwrite=t"

rule remove_low_complexity:
  input:
    Outdir3+"/{sample}_qc2.R1.fastq.gz",
    Outdir3+"/{sample}_qc2.R2.fastq.gz"
  output:
    Outdir3+"/{sample}_qc3.R1.fastq.gz",
    Outdir3+"/{sample}_qc3.R2.fastq.gz"
  shell:
    "bbduk.sh in1={input[0]} in2={input[1]} "
    "out1={output[0]} out2={output[1]} "
    "threads={Threads} "
    "entropy=0.6 "
    "entropywindow=50 "
    "entropyk=5 "
    "overwrite=t"

#### step4 : fastp_2bbduk
use rule fastp_1bwa as fastp_2bbduk with:
  input:
    Outdir3+"/{sample}_qc3.R1.fastq.gz",
    Outdir3+"/{sample}_qc3.R1.fastq.gz"
  output:
    Outdir4+"/{sample}/out.R1.fastq.gz",
    Outdir4+"/{sample}/out.R2.fastq.gz",
    Outdir4+"/{sample}/fastp_out.html",
    Outdir4+"/{sample}/fastp_out.json" 
 
use rule fastp_bwa_stats as fastp_bbduk_stats with:
  input:
    expand(Outdir4+"/{sample}/fastp_out.json",sample=Samples)
  output:
    Outdir8+"/fp_bbduk.txt" 

#### step5 : kneaddata
rule remove_host_reads_quality_control:
  input:
    Outdir4+"/{sample}/out.R1.fastq.gz",
    Outdir4+"/{sample}/out.R2.fastq.gz"
  output:
    directory(Outdir5+"/{sample}"),
    temp(Outdir5+"/{sample}/{sample}_paired_1.fastq"),
    temp(Outdir5+"/{sample}/{sample}_paired_2.fastq"),
    Outdir5+"/{sample}/{sample}_kneaddata_paired_1.fastq",
    Outdir5+"/{sample}/{sample}_kneaddata_paired_2.fastq"
  params:  
    bowtie2_ref=config['database']['bowtie2_ref']
  shell:
    "kneaddata -i1 {input[0]} -i2 {input[1]} "
    "-db {params.bowtie2_ref} -t {Threads} -p 2 --serial --bypass-trf --remove-intermediate-output "
    "-o {output[0]} --output-prefix {wildcards.sample} && "
    "cat {output[1]} > {output[3]} && "
    "cat {output[2]} > {output[4]}"

#### step6 : fastp_3kneaddata   
use rule fastp_1bwa as fastp_3kneaddata with:
  input:
    Outdir5+"/{sample}/{sample}_kneaddata_paired_1.fastq",
    Outdir5+"/{sample}/{sample}_kneaddata_paired_2.fastq"
  output:
    Outdir6+"/{sample}/out.R1.fastq.gz",
    Outdir6+"/{sample}/out.R2.fastq.gz",
    Outdir6+"/{sample}/fastp_out.html",
    Outdir6+"/{sample}/fastp_out.json" 
 
use rule fastp_bwa_stats as fastp_kneaddata_stats with:
  input:
    expand(Outdir6+"/{sample}/fastp_out.json",sample=Samples)
  output:
    Outdir8+"/fp_kneaddata.txt" 

#### step7 : taxa_annotation
rule kraken2_kreport:
  input:
    Outdir6+"/{sample}/out.R1.fastq.gz",
    Outdir6+"/{sample}/out.R2.fastq.gz"
  output: 
    Outdir7+"/{sample}.kreport"
  params:
    Kraken2_db=config['database']["kraken2_ref"],
    confidence=config['parameters']['confidence']
  shell: 
    "kraken2 --db {params.Kraken2_db} --threads {Threads} --confidence {params.confidence} --output - --report {output} "
    "--paired {input[0]} {input[1]}"
    
rule kreport2mpa_taxa:
  input:
    Outdir7+"/{sample}.kreport"
  output: 
    Outdir7+"/{sample}_taxa.count",
    Outdir7+"/{sample}_taxa.percent"
  shell:
    """ 
    kreport2mpa.py -r {input} -o {output[0]} \
    --display-header --no-intermediate-ranks --read_count
    kreport2mpa.py -r {input} -o {output[1]} \
    --display-header --no-intermediate-ranks --percentages 
    """
    
#### step8 : classification
rule count_result:
  input:
    expand(Outdir7+"/{sample}_taxa.count",sample=Samples)
  params:
    config['parameters']['temp']
  output:
    directory(Outdir7+"/all_count_result")
  run:
    shell("mkdir {output} && " 
    "cp {input} {output}/")
    file_dir = str(output)
    files = os.listdir(file_dir)
    print(files)
    data_path = os.path.join(file_dir, files[0])
    DATA = pd.read_csv(data_path, sep='\t')

    for file in files[1:]:
      data_path = os.path.join(file_dir, file)
      data = pd.read_csv(data_path, sep='\t')
      DATA = pd.merge(DATA, data, on="#Classification", how="outer")
      
    data_0 = DATA.fillna(0)
    data_0['count'] = data_0.iloc[:, 1:].apply( sum, axis=1)
    data_0 = data_0.rename(columns={"#Classification": "Classification"})
    data_1 = data_0
    
    result_dir = str(output)
    os.makedirs(result_dir, exist_ok=True)

    def write_filtered_data(df, pattern, exclude_pattern, result_file):
      filtered_df = df[df['Classification'].str.contains(pattern)]
      if exclude_pattern:
        filtered_df = filtered_df[~filtered_df['Classification'].str.contains(exclude_pattern)]
      filtered_df.to_csv(result_file, index=False)
    
    shell("mkdir {output}/Bacteria")
    bacteria = data_1[data_1['Classification'].str.contains("k__Bacteria")]
    write_filtered_data(bacteria, "k__", "p__|c__|o__|f__|g__|s__", os.path.join(result_dir, "Bacteria/All_k_bacteria.csv"))
    write_filtered_data(bacteria, "p__", "c__|o__|f__|g__|s__", os.path.join(result_dir, "Bacteria/All_p_bacteria.csv"))
    write_filtered_data(bacteria, "c__", "o__|f__|g__|s__", os.path.join(result_dir, "Bacteria/All_c_bacteria.csv"))
    write_filtered_data(bacteria, "o__", "f__|g__|s__", os.path.join(result_dir, "Bacteria/All_o_bacteria.csv"))
    write_filtered_data(bacteria, "f__", "g__|s__", os.path.join(result_dir, "Bacteria/All_f_bacteria.csv"))
    write_filtered_data(bacteria, "g__", "s__", os.path.join(result_dir, "Bacteria/All_g_bacteria.csv"))
    write_filtered_data(bacteria, "s__", None, os.path.join(result_dir, "Bacteria/All_s_bacteria.csv"))

    shell("mkdir {output}/Viruses")
    viruses = data_1[data_1['Classification'].str.contains("k__Viruses")]
    write_filtered_data(viruses, "k__", "p__|c__|o__|f__|g__|s__", os.path.join(result_dir, "Viruses/All_k_viruses.csv"))
    write_filtered_data(viruses, "p__", "c__|o__|f__|g__|s__", os.path.join(result_dir, "Viruses/All_p_Viruses.csv"))
    write_filtered_data(viruses, "c__", "o__|f__|g__|s__", os.path.join(result_dir, "Viruses/All_c_Viruses.csv"))
    write_filtered_data(viruses, "o__", "f__|g__|s__", os.path.join(result_dir, "Viruses/All_o_Viruses.csv"))
    write_filtered_data(viruses, "f__", "g__|s__", os.path.join(result_dir, "Viruses/All_f_Viruses.csv"))
    write_filtered_data(viruses, "g__", "s__", os.path.join(result_dir, "Viruses/All_g_Viruses.csv"))
    write_filtered_data(viruses, "s__", None, os.path.join(result_dir, "Viruses/All_s_Viruses.csv"))

    shell("mkdir {output}/Archaea")
    archaea = data_1[data_1['Classification'].str.contains("k__Archaea")]
    write_filtered_data(archaea, "k__", "p__|c__|o__|f__|g__|s__", os.path.join(result_dir, "Archaea/All_k_archaea.csv"))
    write_filtered_data(archaea, "p__", "c__|o__|f__|g__|s__", os.path.join(result_dir, "Archaea/All_p_Archaea.csv"))
    write_filtered_data(archaea, "c__", "o__|f__|g__|s__", os.path.join(result_dir, "Archaea/All_c_Archaea.csv"))
    write_filtered_data(archaea, "o__", "f__|g__|s__", os.path.join(result_dir, "Archaea/All_o_Archaea.csv"))
    write_filtered_data(archaea, "f__", "g__|s__", os.path.join(result_dir, "Archaea/All_f_Archaea.csv"))
    write_filtered_data(archaea, "g__", "s__", os.path.join(result_dir, "Archaea/All_g_Archaea.csv"))
    write_filtered_data(archaea, "s__", None, os.path.join(result_dir, "Archaea/All_s_Archaea.csv"))

    shell("mkdir {output}/Eukaryota_Fungi")
    fungi = data_1[data_1['Classification'].str.contains("k__Fungi")]
    write_filtered_data(fungi, "k__", "p__|c__|o__|f__|g__|s__", os.path.join(result_dir, "Eukaryota_Fungi/All_k_fungi.csv"))
    write_filtered_data(fungi, "p__", "c__|o__|f__|g__|s__", os.path.join(result_dir, "Eukaryota_Fungi/All_p_Fungi.csv"))
    write_filtered_data(fungi, "c__", "o__|f__|g__|s__", os.path.join(result_dir, "Eukaryota_Fungi/All_c_Fungi.csv"))
    write_filtered_data(fungi, "o__", "f__|g__|s__", os.path.join(result_dir, "Eukaryota_Fungi/All_o_Fungi.csv"))
    write_filtered_data(fungi, "f__", "g__|s__", os.path.join(result_dir, "Eukaryota_Fungi/All_f_Fungi.csv"))
    write_filtered_data(fungi, "g__", "s__", os.path.join(result_dir, "Eukaryota_Fungi/All_g_Fungi.csv"))
    write_filtered_data(fungi, "s__", None, os.path.join(result_dir, "Eukaryota_Fungi/All_s_Fungi.csv"))

    shell("mkdir {output}/Eukaryota_Metazoa")
    Metazoa = data_1[data_1['Classification'].str.contains("k__Metazoa")]
    write_filtered_data(Metazoa, "k__", "p__|c__|o__|f__|g__|s__", os.path.join(result_dir, "Eukaryota_Metazoa/All_k_Metazoa.csv"))
    write_filtered_data(Metazoa, "p__", "c__|o__|f__|g__|s__", os.path.join(result_dir, "Eukaryota_Metazoa/All_p_Metazoa.csv"))
    write_filtered_data(Metazoa, "c__", "o__|f__|g__|s__", os.path.join(result_dir, "Eukaryota_Metazoa/All_c_Metazoa.csv"))
    write_filtered_data(Metazoa, "o__", "f__|g__|s__", os.path.join(result_dir, "Eukaryota_Metazoa/All_o_Metazoa.csv"))
    write_filtered_data(Metazoa, "f__", "g__|s__", os.path.join(result_dir, "Eukaryota_Metazoa/All_f_Metazoa.csv"))
    write_filtered_data(Metazoa, "g__", "s__", os.path.join(result_dir, "Eukaryota_Metazoa/All_g_Metazoa.csv"))
    write_filtered_data(Metazoa, "s__", None, os.path.join(result_dir, "Eukaryota_Metazoa/All_s_Metazoa.csv"))

    shell("mkdir {output}/other_Eukaryota")
    eukaryota = data_1[data_1['Classification'].str.contains("k__Eukaryota") & ~data_1['Classification'].str.contains("k__Fungi|k__Metazoa")]
    write_filtered_data(eukaryota, "k__", "p__|c__|o__|f__|g__|s__", os.path.join(result_dir, "other_Eukaryota/All_k_other_Eukaryota.csv"))
    write_filtered_data(eukaryota, "p__", "c__|o__|f__|g__|s__", os.path.join(result_dir, "other_Eukaryota/All_p_other_Eukaryota.csv"))
    write_filtered_data(eukaryota, "c__", "o__|f__|g__|s__", os.path.join(result_dir, "other_Eukaryota/All_c_other_Eukaryota.csv"))
    write_filtered_data(eukaryota, "o__", "f__|g__|s__", os.path.join(result_dir, "other_Eukaryota/All_o_other_Eukaryota.csv"))
    write_filtered_data(eukaryota, "f__", "g__|s__", os.path.join(result_dir, "other_Eukaryota/All_f_other_Eukaryota.csv"))
    write_filtered_data(eukaryota, "g__", "s__", os.path.join(result_dir, "other_Eukaryota/All_g_other_Eukaryota.csv"))
    write_filtered_data(eukaryota, "s__", None, os.path.join(result_dir, "other_Eukaryota/All_s_other_Eukaryota.csv"))
    
    shell("""
          rm {output}/*.count
          """)
    shell("""      
           if {params};then
             rm -r {Outdir1} {Outdir2} {Outdir3} {Outdir4} {Outdir5} {Outdir6}
           fi
          """)