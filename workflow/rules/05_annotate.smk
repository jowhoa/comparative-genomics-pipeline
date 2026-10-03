# 1. Fetch paired-end RNA-seq from SRA
rule download_rnaseq:
    output:
        r1="data/rna/raw/{sra_id}_1.fastq",
        r2="data/rna/raw/{sra_id}_2.fastq"
    threads: 6
    resources:
        mem_mb=8000,
        runtime=1440
    conda:
        "../envs/sra_tools.yaml"
    shell:
        "fasterq-dump {wildcards.sra_id} --split-files -O data/rna/raw/ -e {threads}"

# 2. Trim paired-end RNA-seq reads
rule trim_rnaseq:
    input:
        r1="data/rna/raw/{sra_id}_1.fastq",
        r2="data/rna/raw/{sra_id}_2.fastq"
    output:
        trimmed_r1="data/rna/trimmed/{sra_id}_1_trimmed.fastq",
        trimmed_r2="data/rna/trimmed/{sra_id}_2_trimmed.fastq",
        html="results/qc/fastp/{sra_id}_fastp.html"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=1440
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.trimmed_r1} -O {output.trimmed_r2} -h {output.html} --thread {threads}"

# 3. Align, Sort, and Index BAM
rule map_and_sort_rnaseq:
    input:
        genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta",
        r1="data/rna/trimmed/{sra_id}_1_trimmed.fastq",
        r2="data/rna/trimmed/{sra_id}_2_trimmed.fastq"
    output:
        bam="results/annotation/{sample}_{assembly_id}/{sra_id}_mapped.bam",
        bai="results/annotation/{sample}_{assembly_id}/{sra_id}_mapped.bam.bai"
    threads: 24
    resources:
        mem_mb=64000,
        runtime=1440
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        minimap2 -ax splice -t {threads} {input.genome} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """
# 4. Download and Combine Protein FASTA
rule prepare_protein_evidence:
    output:
        combined_prot="data/proteins/Phytophthora_combined_proteins.fasta"
    conda:
        "../envs/ncbi_datasets.yaml"
    shell:
        """
        mkdir -p data/proteins/raw
        datasets download genome accession GCA_020800215.2 GCA_031305395.1 GCF_018691715.1 \
            --include protein \
            --filename data/proteins/raw/phytophthora_proteins.zip
        unzip -o data/proteins/raw/phytophthora_proteins.zip -d data/proteins/raw/
        cat data/proteins/raw/ncbi_dataset/data/*/*.faa > {output.combined_prot}
        """

# 5. Execute BRAKER3
rule run_braker:
    input:
        genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta",
        bams=expand("results/annotation/{{sample}}_{{assembly_id}}/{sra_id}_mapped.bam", sra_id=config["sra_accessions"]),
        proteins="data/proteins/Phytophthora_combined_proteins.fasta"
    output:
        gff3="results/annotation/{sample}_{assembly_id}/braker/braker.gff3",
        aa="results/annotation/{sample}_{assembly_id}/braker/braker.aa"
    params:
        species=config["species_name"]
    threads: 24
    resources:
        mem_mb=256000,
        runtime=4320
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        BAM_LIST=$(echo {input.bams} | tr ' ' ',')
        
        braker.pl \
            --genome={input.genome} \
            --prot_seq={input.proteins} \
            --bam=$BAM_LIST \
            --species={params.species}_{wildcards.sample}_{wildcards.assembly_id} \
            --workingdir=results/annotation/{wildcards.sample}_{wildcards.assembly_id}/braker \
            --threads={threads}
        """
