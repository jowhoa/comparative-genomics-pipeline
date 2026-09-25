rule run_quast:
    input:
        "results/01_assembly/{sample}.fasta"
    output:
        "results/01_assembly/QC/quast_{sample}/report.txt"
    conda:
        "envs/qc.yaml"
    shell:
        "quast.py {input} -o results/01_assembly/QC/quast_{wildcards.sample} --large"

rule run_busco:
    input:
        "results/01_assembly/{sample}.fasta"
    output:
        "results/01_assembly/QC/busco_{sample}/short_summary.txt"
    params:
        lineage=config["busco_lineage"]
    conda:
        "envs/qc.yaml"
    shell:
        """
        # 1. Run the normal BUSCO command
        busco -i {input} -o busco_{wildcards.sample} --out_path results/01_assembly/QC -l eukaryota_odb10 -m genome --f>

        # 2. Rename the dynamically generated file to match what Snakemake expects
        mv results/01_assembly/QC/busco_{wildcards.sample}/short_summary.specific.*.txt {output}
        """

rule run_blast:
    input:
        assembly="results/01_assembly/{sample}.fasta",
        reference=config["reference_genome"]
    output:
        blast_out="results/01_assembly/QC/blast_{sample}/{sample}.blast.out"
    params:
        # We tell it to build the temporary database inside the results folder
        db_prefix="results/01_assembly/QC/blast_{sample}/ref_db"
    threads:
        8
    conda:
        "envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input.reference} -dbtype nucl -out {params.db_prefix}

        blastn -query {input.assembly} -db {params.db_prefix} -outfmt 6 -out {output.blast_out} -num_threads {threads}
        """

rule generate_dotplot:
    input:
        blast="results/01_assembly/QC/blast_{sample}/{sample}.blast.out"
    output:
        plot="results/01_assembly/QC/blast_{sample}/{sample}_dotplot.png"
    shell:
        """
        # 1. Automatically shrink the massive BLAST file
        awk '$4 > 10000' {input.blast} > {input.blast}.tmp_filtered

        # 2. Pass the files directly to your Python script
        python make_dotplot.py {input.blast}.tmp_filtered {output.plot}

        # 3. Clean up the temporary file to save storage space
        rm {input.blast}.tmp_filtered
        """
