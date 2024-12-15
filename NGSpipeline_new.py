import os
import subprocess
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def check_and_install_tool(tool_name, install_command):
    """Check if a tool is installed, and install it if not."""
    if tool_name == "featureCounts":
        try:
            subprocess.run([tool_name, "-v"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            print(f"{tool_name} is already installed.")
        except subprocess.CalledProcessError:
            print(f"{tool_name} is not installed. Installing...")
            try:
                subprocess.run(install_command, shell=True, check=True)
                print(f"{tool_name} installed successfully.")
            except Exception as e:
                print(f"Failed to install {tool_name}. Error: {e}")
                raise
    else:
            try:
                subprocess.run([tool_name, "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                print(f"{tool_name} is already installed.")
            except subprocess.CalledProcessError:
                print(f"{tool_name} is not installed. Installing...")
                try:
                    subprocess.run(install_command, shell=True, check=True)
                    print(f"{tool_name} installed successfully.")
                except Exception as e:
                    print(f"Failed to install {tool_name}. Error: {e}")
                    raise

def fetch_data(directory, filename_list):
    for file in filename_list:
        file_path = os.path.join(directory, file)
        if os.path.exists(file_path):
            print(f"File '{file}' already exists")
        else:
            try:
                subprocess.run(f"prefetch --output-directory {directory} {file}", shell=True, check=True)
            except Exception as e:
                print(f"File download failed.\nError: {e}")

def fastq_dump(directory, filename_list):
    print("Starting fastq-dump...")
    for file in filename_list:
        file_path = os.path.join(directory, file)
        if os.path.exists(f'{file_path}_1.fastq'):
            print(f"{file}.fastq files already exist")
        else:
            command = f"fasterq-dump {directory}{file}/{file}.sra --outdir {directory}"
            subprocess.run(command, shell=True, check=True)

def quality_check_1(directory, filename_list):
    print("Generating Quality Report...")
    for file in filename_list:
        file_path = os.path.join(directory, file)
        for end in [1, 2]:  # Assumes paired-end data
            output_html = f"{file_path}_{end}_fastqc.html"
            if os.path.exists(output_html):
                print(f"{output_html} already exists")
            else:
                command = f"fastqc {directory}{file}_{end}.fastq -o {directory}"
                subprocess.run(command, shell=True, check=True)

def trimming(filename_list, directory, threads, trimmomatic_path):
    for file in filename_list:
        file_path = os.path.join(directory, file)
        if os.path.exists(f'{file_path}_trim_forward_paired.fastq'):
            print('trimmed files already present')
        else:
            cmd = f"java -jar {trimmomatic_path} PE -threads {threads} {file_path}_1.fastq {file_path}_2.fastq {file_path}_trim_forward_paired.fastq {file_path}_trim_forward_unpaired.fastq {file_path}_trim_reverse_paired.fastq {file_path}_trim_reverse_unpaired.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"
            subprocess.run(cmd,shell=True,check=True)

def quality_check_2(directory, filename_list):
    print("Generating Quality Report...")
    for file in filename_list:
        file_path = os.path.join(directory, file)
        for end in ['_trim_forward_paired','_trim_reverse_paired']:  
            output_html = f"{file_path}{end}_fastqc.html"
            if os.path.exists(output_html):
                print(f"{output_html} already exists")
            else:
                command = f"fastqc {file_path}{end}.fastq -o {directory}"
                subprocess.run(command, shell=True, check=True)

def indexing_genome(directory, reference_genome_dna,threads):
    print("Starting genome indexing...")
    index_path = os.path.join(directory, "index_gene.1.ht2")
    if os.path.exists(index_path):
        print("Reference genome index already exists")
    else:
        command = f"hisat2-build -p {threads} {reference_genome_dna} {directory}index_gene"
        subprocess.run(command, shell=True, check=True)

def run_hisat2(directory, filename_list,threads):
    print("Running alignment mapping...")
    for file in filename_list:
        file_path = os.path.join(directory, file)
        sam_file = f"{file_path}.sam"
        if os.path.exists(sam_file):
            print(f"Aligned SAM file {sam_file} already exists")
        else:
            cmd1 = f"hisat2 -x {directory}index_gene -U {file_path}_trim_forward_paired.fastq,{file_path}_trim_reverse_paired.fastq -S {sam_file} -p {threads}"
            subprocess.run(cmd1, shell=True, check=True)

        bam_file = f"{file_path}.bam"
        if os.path.exists(bam_file):
            print(f"Converted BAM file {bam_file} already exists")
        else:
            cmd2 = f"samtools view -S -b {sam_file} > {bam_file}"
            subprocess.run(cmd2, shell=True, check=True)

        sorted_bam = f"{file_path}_sort.bam"
        if os.path.exists(sorted_bam):
            print(f"Sorted BAM file {sorted_bam} already exists")
        else:
            cmd3 = f"samtools sort {bam_file} -o {sorted_bam}"
            subprocess.run(cmd3, shell=True, check=True)

def feature_count(reference_genome_annotation, directory, filename_list, threads):
    print("Starting feature count...")
    for file in filename_list:
        file_path = os.path.join(directory, file)
        output_file = f"{file_path}_fc"
        if os.path.exists(output_file):
            print(f"Feature count file {output_file} already exists")
        else:
            cmd = f"featureCounts -a {reference_genome_annotation} -o {output_file} -T {threads} {file_path}_sort.bam -g db_xref"
            subprocess.run(cmd, shell=True, check=True)

def featurecount_matrix(input):
    gene_id = []
    gene_exp = []
    with open(input) as file:
        for e, line in enumerate(file):
            if e < 2: 
                continue
            line_elements = line.split()
            gene_id.append(line_elements[0])
            gene_exp.append(line_elements[-1])
    return np.array(gene_exp), np.array(gene_id)

def matrix_creation(directory, filename_list,expression_level):
    if os.path.exists(f'{directory}filtered_matrix.csv'):
        print("final matrix file is already present")
    else:
        gene_exp_matrix = []
        gene_id = None
        for file in filename_list:
            file_path = os.path.join(directory, f"{file}_fc")
            gene_exp, gene_id = featurecount_matrix(file_path)
            gene_exp_matrix.append(gene_exp)

        gene_exp_matrix = np.array(gene_exp_matrix).T
        matrix = pd.DataFrame(gene_exp_matrix, index=gene_id, columns=filename_list)
        matrix.to_csv(os.path.join(directory, "unfiltered_matrix.csv"))
        unfiltered = pd.read_csv(os.path.join(directory, "unfiltered_matrix.csv"),index_col=0)
        filtered_matrix  = unfiltered[unfiltered.sum(axis=1) > expression_level]
        filtered_matrix = filtered_matrix
        filtered_matrix.to_csv(os.path.join(directory, "filtered_matrix.csv"))

def DEGS_analysis(directory,metadata_table):
    if os.path.exists(f'{directory}Degs_result.csv'):
        print("DESeq file is already present.")
    else:
        cmd = "pip install pydeseq2"
        subprocess.run(cmd,shell=True,check=True)
        counts_df = pd.read_csv(f"{directory}1filtered_matrix.csv", index_col=0).T
        counts_df = counts_df.round().astype(int)
        clinical_df = pd.read_csv(metadata_table, index_col=0)

        dds = DeseqDataSet(counts=counts_df, metadata=clinical_df, design_factors=["condition"])

        dds.deseq2()

        contrast = ["condition", "treated", "untreated"]  
        ds = DeseqStats(dds=dds, contrast=contrast)

        results_df = ds.summary()
        results_df = pd.DataFrame(results_df)
        print(results_df)
        results_df.to_csv(f"{directory}Degs_result.csv")

def main():
    tools = {
        "prefetch": "sudo apt-get install -y sra-toolkit",
        "fastqc": "sudo apt-get install -y fastqc",
        "hisat2": "sudo apt-get install -y hisat2",
        "samtools": "sudo apt-get install -y samtools",
        "featureCounts": "sudo apt-get install -y subread"
    }
    for tool, install_cmd in tools.items():
        check_and_install_tool(tool, install_cmd)
    fetch_data(directory, filename_list)
    fastq_dump(directory, filename_list)
    quality_check_1(directory, filename_list)
    trimmomatic_path = input("Enter path of your trimmomatic file you have downloaded")
    trimming(filename_list, directory, threads,trimmomatic_path)
    quality_check_2(directory, filename_list)
    while True:
        proceed = input("----Please Analyse Quality Report and overwrite existing QC file if needed.----\nContinue? (yes/no)\n")
        if proceed.lower() == 'yes':
            indexing_genome(directory, reference_genome_dna,threads)
            run_hisat2(directory, filename_list,threads)
            feature_count(reference_genome_annotation, directory, filename_list,threads)
            expression_level = int(input("Enter total expression level on the basis you want to filter feature count matrix.\n"))
            matrix_creation(directory, filename_list,expression_level)
            metadata_table = input("Enter path of your metadata table. (row should contain sample id, column should contain Conditions)")
            DEGS_analysis(directory,metadata_table)
            break
        elif proceed.lower() == 'no':
            print("Pipeline terminated as requested.")
            return
        else:
            print("Invalid input. Please TRY AGAIN.")

def get_valid_path(prompt):
    while True:
        path = input(prompt)
        if os.path.exists(path):
            return path
        else:
            print(f"File or directory '{path}' not found. Please try again.")

if __name__ == "__main__":
    print('''This is an NGS Pipeline. It will use 10 CPU threads.
    Requirements:
    1. java version 11
    2. Make a text file containing all your SRR IDs and provide the path.
    3. Provide the path of the directory where you want outputs.
    4. Provide the path of the reference genome (fna) file.
    5. Provide the path of the reference genome annotation file.
    6. you should download trimmomatic jar file''')
    
    # INPUTS
    textfile = input("Enter path of text file contain sample ID (SRR ID).\n")
    input_directory = input("Enter path of output directory.\n")
    directory = os.path.join(input_directory, '')
    reference_genome_dna = input("Enter path of reference genome.\n")
    reference_genome_annotation = input("Enter path of reference genome's annotation file.\n")
    threads = int(input("Enter no. of threads you want this pipeline use...\n"))

    filename_list = []
    with open(textfile, 'r') as f:
        for line in f:
            filename = line.strip()
            filename_list.append(filename)

    main()
