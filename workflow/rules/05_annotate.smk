rule run_hite:
    input:
        genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    output:
        te_lib="results/hite/{sample}_{assembly_id}/confident_TE.cons.fa"
    conda:
        "../envs/hite.yaml"
    threads: 48
    resources:
        mem_mb=64000,
        runtime=1440
    log:
        "results/logs/slurm/rule_run_hite_{sample}_{assembly_id}.log"
    shell:
        """
        python /quobyte/fzmichelgrp/undergrads/Jonas/comparative_project/HiTE/main.py \
            --genome {input.genome} \
            --thread {threads} \
            --out_dir results/hite/{wildcards.sample}_{wildcards.assembly_id}/ > {log} 2>&1
        """
rule repeat_masker:
    input:
        genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta",
        te_lib="results/hite/{sample}_{assembly_id}/confident_TE.cons.fa"
    output:
        masked_genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta.masked"
    conda:
        "../envs/repeatmasker.yaml"
    threads: 48
    resources:
        mem_mb=32000,
        runtime=1440
    log:
        "results/logs/slurm/rule_repeat_masker_{sample}_{assembly_id}.log"
    shell:
        """
        RepeatMasker -pa {threads} -xsmall -lib {input.te_lib} \
        -dir results/scaffolds/{wildcards.sample}_{wildcards.assembly_id}/ \
        {input.genome} > {log} 2>&1
        """
# 1. Fetch paired-end RNA-seq from SRA
rule download_rnaseq:
    output:
        r1="data/rna/raw/{sra_id}_1.fastq",
        r2="data/rna/raw/{sra_id}_2.fastq"
    log:
        "results/logs/slurm/rule_download_rnaseq_{sra_id}.log"
    threads: 6
    resources:
        mem_mb=8000,
        runtime=1440
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        fasterq-dump {wildcards.sra_id} --split-files -O data/rna/raw/ -e {threads} > {log} 2>&1
        """

# 2. Trim paired-end RNA-seq reads
rule trim_rnaseq:
    input:
        r1="data/rna/raw/{sra_id}_1.fastq",
        r2="data/rna/raw/{sra_id}_2.fastq"
    output:
        trimmed_r1="data/rna/trimmed/{sra_id}_1_trimmed.fastq",
        trimmed_r2="data/rna/trimmed/{sra_id}_2_trimmed.fastq",
        html="results/qc/fastp/{sra_id}_fastp.html"
    log:
        "results/logs/slurm/rule_trim_rnaseq_{sra_id}.log"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=1440
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.trimmed_r1} -O {output.trimmed_r2} \
            -h {output.html} --thread {threads} > {log} 2>&1
        """

# 3. Align, Sort, and Index BAM
rule map_and_sort_rnaseq:
    input:
        genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta",
        r1="data/rna/trimmed/{sra_id}_1_trimmed.fastq",
        r2="data/rna/trimmed/{sra_id}_2_trimmed.fastq"
    output:
        bam="results/annotation/{sample}_{assembly_id}/{sra_id}_mapped.bam",
        bai="results/annotation/{sample}_{assembly_id}/{sra_id}_mapped.bam.bai"
    log:
        "results/logs/slurm/rule_map_and_sort_rnaseq_{sample}_{assembly_id}_{sra_id}.log"
    threads: 24
    resources:
        mem_mb=64000,
        runtime=1440
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        # Pipe stderr of minimap2 and samtools to the log file so we catch any mapping errors
        minimap2 -ax splice -t {threads} {input.genome} {input.r1} {input.r2} 2> {log} | \
        samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        
        samtools index {output.bam} 2>> {log}
        """

# 4. Download and Combine Protein FASTA
rule prepare_protein_evidence:
    output:
        combined_prot="data/proteins/Phytophthora_combined_proteins.fasta"
    log:
        "results/logs/slurm/rule_prepare_protein_evidence.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/ncbi_datasets.yaml"
    shell:
        """
        mkdir -p data/proteins/raw
        datasets download genome accession GCA_020800215.2 GCA_031305395.1 GCF_018691715.1 \
            --include protein \
            --filename data/proteins/raw/phytophthora_proteins.zip > {log} 2>&1
            
        unzip -o data/proteins/raw/phytophthora_proteins.zip -d data/proteins/raw/ >> {log} 2>&1
        cat data/proteins/raw/ncbi_dataset/data/*/*.faa > {output.combined_prot} 2>> {log}
        """

# 5. Execute BRAKER3
rule run_braker:
    input:
        genome="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta.masked",
        bams=expand("results/annotation/{{sample}}_{{assembly_id}}/{sra_id}_mapped.bam", sra_id=config["sra_accessions"]),
        proteins="data/proteins/Phytophthora_combined_proteins.fasta"
    output:
        gtf = "results/annotation/{sample}_{assembly_id}/braker/braker.gtf",
        codingseq = "results/annotation/{sample}_{assembly_id}/braker/braker.codingseq"
    log:
        "results/logs/slurm/rule_run_braker_{sample}_{assembly_id}.log"
    params:
        species=config["species_name"]
    threads: 48
    resources:
        mem_mb=256000,
        runtime=4320
    shell:
        """
        # 1. Safely split all inputs by space, filter strictly for .bam, and join with commas
        BAM_LIST=$(echo {input.bams} | tr ' ' '\\n' | grep '\\.bam$' | tr '\\n' ',')
        BAM_LIST=${{BAM_LIST%,}}

        # 2. Define directories
        LOCAL_CONFIG="$PWD/results/annotation/{wildcards.sample}_{wildcards.assembly_id}/augustus_config"
        BRAKER_DIR="results/annotation/{wildcards.sample}_{wildcards.assembly_id}/braker"

        # 3. FORCE CLEANUP: Delete dirty state from previous failed runs
        rm -rf $LOCAL_CONFIG $BRAKER_DIR
        mkdir -p $LOCAL_CONFIG

        # 4. ISOLATION FIX: Prevent local Python packages in ~/.local from contaminating the container
        export APPTAINERENV_PYTHONNOUSERSITE=1

        # 5. Evaluate the wildcard INSIDE the container and explicitly bind the Quobyte drive
        apptainer exec -B $PWD:$PWD /quobyte/fzmichelgrp/undergrads/Jonas/braker3.sif bash -c "cp -r /opt/Augustus/config/* $LOCAL_CONFIG/"

        # 6. Force Apptainer to use your local writeable directory instead of its internal one
        export APPTAINERENV_AUGUSTUS_CONFIG_PATH="$LOCAL_CONFIG"

        # 7. Execute BRAKER3
        apptainer exec -B $PWD:$PWD -B ~/.gm_key:/root/.gm_key /quobyte/fzmichelgrp/undergrads/Jonas/braker3.sif braker.pl \
            --genome={input.genome} \
            --prot_seq={input.proteins} \
            --bam=$BAM_LIST \
            --species=Phytophthora_infestans_{wildcards.sample}_{wildcards.assembly_id} \
            --workingdir=$BRAKER_DIR \
            --threads={threads} > {log} 2>&1
        """
