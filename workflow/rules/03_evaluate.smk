rule run_quast:
    """
    Evaluates assembly contiguity (N50, total length).
    """
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta"
    output:
        html="results/03_evaluate/quast/{sample}_{assembly_id}/{sample}_{assembly_id}_report.html",
        tsv="results/03_evaluate/quast/{sample}_{assembly_id}/{sample}_{assembly_id}_report.tsv"
    log:
        "results/logs/slurm/rule_quast_{sample}_{assembly_id}.log"
    conda:
        "../envs/quast.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=360
    params:
        outdir="results/03_evaluate/quast/{sample}_{assembly_id}"
    shell:
        """
        # 1. Run QUAST (it will create the generic report.html and report.tsv)
        quast.py {input} -o {params.outdir} -t {threads} --large > {log} 2>&1

        # 2. Rename the generic files to your specific output targets
        mv {params.outdir}/report.html {output.html}
        mv {params.outdir}/report.tsv {output.tsv}
        """

rule run_busco:
    """
    Assesses genome completeness using conserved single-copy orthologs.
    """
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta"
    output:
        summary="results/03_evaluate/busco/{sample}_{assembly_id}/short_summary.txt"
    log:
        "results/logs/slurm/rule_busco_{sample}_{assembly_id}.log"
    conda:
        "../envs/busco.yaml"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=1440
    params:
        outdir="results/03_evaluate/busco",
        lineage="stramenopiles_odb10",
        prefix="{sample}_{assembly_id}"
    shell:
        """
        # BUSCO writes to the current directory by default if --out_path isn't strictly handled,
        # so we force it into the correct results folder.
        busco -i {input.fasta} -o {params.prefix} --out_path {params.outdir} \
            -l {params.lineage} -m geno -c {threads} -f > {log} 2>&1
            
        # Move the summary to match Snakemake's exact expected output path
        cp {params.outdir}/{params.prefix}/short_summary.specific.{params.lineage}.{params.prefix}.txt {output.summary}
        """

rule meryl_count:
    """
    Builds the exact k-mer database required by Merqury from the filtered HiFi reads.
    """
    input:
        reads="results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        meryldb=directory("results/00_qc/meryl/{sample}.meryl")
    conda:
        "../envs/merqury.yaml"
    log:
        "results/logs/slurm/rule_meryl_{sample}.log"
    threads: 16
    resources:
        mem_mb=32000,
        runtime=120
    shell:
        """
        meryl k=21 count threads={threads} memory=32 {input.reads} output {output.meryldb} > {log} 2>&1
        """

rule run_merqury:
    """
    Calculates Assembly Quality Value (QV) and k-mer completeness.
    """
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta",
        meryldb="results/00_qc/meryl/{sample}.meryl"
    output:
        qv="results/03_evaluate/merqury/{sample}_{assembly_id}/{sample}_{assembly_id}.qv",
        completeness="results/03_evaluate/merqury/{sample}_{assembly_id}/{sample}_{assembly_id}.completeness.stats"
    log:
        "results/logs/slurm/rule_merqury_{sample}_{assembly_id}.log"
    conda:
        "../envs/merqury.yaml"
    threads: 8
    resources:
        mem_mb=48000,
        runtime=120
    params:
        outdir="results/03_evaluate/merqury/{sample}_{assembly_id}",
        prefix="{sample}_{assembly_id}"
    shell:
        """
        # Capture the absolute path of your main project folder before moving
        ROOT_DIR=$(pwd)

        mkdir -p {params.outdir}
        cd {params.outdir}

        # Use ROOT_DIR to ensure the script always finds the exact files
        merqury.sh $ROOT_DIR/{input.meryldb} $ROOT_DIR/{input.fasta} {params.prefix} > $ROOT_DIR/{log} 2>&1
        """

rule run_blast:
    input:
        assembly="results/01_assembly/{sample}.fasta",
        reference=config["reference_genome"]
    output:
        blast_out="results/01_assembly/QC/blast_{sample}/{sample}.blast.out"
    log:
        "results/logs/slurm/rule_blast_{sample}.log"
    threads: 24
    resources:
        mem_mb=32000,
        runtime=360
    conda:
        "../envs/blast.yaml"
    params:
        db_prefix="results/01_assembly/QC/blast_{sample}/ref_db"
    shell:
        """
        makeblastdb -in {input.reference} -dbtype nucl -out {params.db_prefix} > {log} 2>&1
        blastn -query {input.assembly} -db {params.db_prefix} -outfmt 6 -out {output.blast_out} -num_threads {threads} >> {log} 2>&1
        """

rule generate_dotplot:
    input:
        blast="results/01_assembly/QC/blast_{sample}/{sample}.blast.out"
    output:
        plot="results/01_assembly/QC/blast_{sample}/{sample}_dotplot.png"
    log:
        "results/logs/slurm/rule_generate_dotplot_{sample}.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        # 1. Automatically shrink the massive BLAST file
        awk '$4 > 10000' {input.blast} > {input.blast}.tmp_filtered

        # 2. Pass the files directly to your Python script with explicit flags
        python workflow/scripts/plot_alignments.py -i {input.blast}.tmp_filtered -o {output.plot} > {log} 2>&1

        # 3. Clean up the temporary file to save storage space
        rm {input.blast}.tmp_filtered
        """
