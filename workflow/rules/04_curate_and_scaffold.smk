rule run_kraken2:
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta"
    output:
        kraken_out="results/qc/kraken/{sample}_{assembly_id}.kraken",
        kraken_report="results/qc/kraken/{sample}_{assembly_id}.report"
    log:
        "results/logs/slurm/rule_kraken2_{sample}_{assembly_id}.log"
    params:
        db=config["kraken_db"]
    threads: 24
    resources:
        mem_mb=64000,
        runtime=360
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 --db {params.db} \
            --threads {threads} \
            --report {output.kraken_report} \
            --output {output.kraken_out} \
            {input.fasta} > {log} 2>&1
        """

rule extract_clean_contigs:
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta",
        kraken_out="results/qc/kraken/{sample}_{assembly_id}.kraken",
        kraken_report="results/qc/kraken/{sample}_{assembly_id}.report"
    output:
        clean_fasta="results/assembly_clean/{sample}_{assembly_id}_clean.fasta"
    log:
        "results/logs/slurm/rule_extract_clean_contigs_{sample}_{assembly_id}.log"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        extract_kraken_reads.py \
            -k {input.kraken_out} \
            -r {input.kraken_report} \
            -s {input.fasta} \
            -o {output.clean_fasta} \
            -t 4762 0 \
            --include-children > {log} 2>&1
        """

rule run_ragtag:
    input:
        clean_fasta="results/assembly_clean/{sample}_{assembly_id}_clean.fasta",
        reference=config["reference_genome"]
    output:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    log:
        "results/logs/slurm/rule_ragtag_{sample}_{assembly_id}.log"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120
    conda:
        "../envs/scaffold.yaml"
    shell:
        """
        # RagTag automatically creates its own output directory. 
        # We point it to the results/scaffolds folder.
        ragtag.py scaffold {input.reference} {input.clean_fasta} \
            -o results/scaffolds/{wildcards.sample}_{wildcards.assembly_id} \
            -t {threads} > {log} 2>&1
        """

rule find_telomeres:
    input:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    output:
        telomere_plot="results/telomeres/{sample}_{assembly_id}/{sample}_{assembly_id}_telomeric_repeat_windows.tsv"
    log:
        "results/logs/slurm/rule_find_telomeres_{sample}_{assembly_id}.log"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/scaffold.yaml"
    shell:
        """
        # tidk searches for the TTTAGGG motif in windows across your new chromosomes
        tidk search \
            --string TTTAGGG \
            --dir results/telomeres/{wildcards.sample}_{wildcards.assembly_id} \
            --output {wildcards.sample}_{wildcards.assembly_id} \
            --extension tsv \
            {input.scaffolds} > {log} 2>&1
        """

rule make_blast_db_reference:
    input:
        reference=config["reference_genome"]
    output:
        ref_db="results/rbh/dbs/reference/reference.ndb"
    log:
        "results/logs/slurm/rule_make_blast_db_reference.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        makeblastdb -in {input.reference} -dbtype nucl -out results/rbh/dbs/reference/reference > {log} 2>&1
        """

rule make_blast_db_scaffolds:
    input:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    output:
        scaffold_db="results/rbh/dbs/{sample}_{assembly_id}/{sample}_{assembly_id}.ndb"
    log:
        "results/logs/slurm/rule_make_blast_db_scaffolds_{sample}_{assembly_id}.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        makeblastdb -in {input.scaffolds} -dbtype nucl -out results/rbh/dbs/{wildcards.sample}_{wildcards.assembly_id}/{wildcards.sample}_{wildcards.assembly_id} > {log} 2>&1
        """

rule run_reciprocal_blast:
    input:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta",
        reference=config["reference_genome"],
        scaffold_db="results/rbh/dbs/{sample}_{assembly_id}/{sample}_{assembly_id}.ndb",
        ref_db="results/rbh/dbs/reference/reference.ndb"
    output:
        q_to_r="results/rbh/blast/{sample}_{assembly_id}_to_ref.tsv",
        r_to_q="results/rbh/blast/ref_to_{sample}_{assembly_id}.tsv"
    log:
        "results/logs/slurm/rule_reciprocal_blast_{sample}_{assembly_id}.log"
    threads: 24
    resources:
        mem_mb=64000,
        runtime=1440
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        blastn -query {input.scaffolds} -db results/rbh/dbs/reference/reference \
            -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads {threads} -out {output.q_to_r} > {log} 2>&1
            
        blastn -query {input.reference} -db results/rbh/dbs/{wildcards.sample}_{wildcards.assembly_id}/{wildcards.sample}_{wildcards.assembly_id} \
            -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads {threads} -out {output.r_to_q} >> {log} 2>&1
        """

# 1. Generate the strict 1:1 ortholog table (No plotting)
rule parse_rbh:
    input:
        q_to_r="results/rbh/blast/{sample}_{assembly_id}_to_ref.tsv",
        r_to_q="results/rbh/blast/ref_to_{sample}_{assembly_id}.tsv"
    output:
        rbh_table="results/rbh/{sample}_{assembly_id}_reciprocal_best_hits.tsv"
    log:
        "results/logs/slurm/rule_parse_rbh_{sample}_{assembly_id}.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        python workflow/scripts/find_rbh.py -q {input.q_to_r} -r {input.r_to_q} -o {output.rbh_table} > {log} 2>&1
        """

# 2. Plot the whole-genome alignments using the raw data
rule plot_raw_alignments:
    input:
        # Feeding the raw, unfiltered BLAST alignments to the script instead of the RBH table
        raw_alignments="results/rbh/blast/{sample}_{assembly_id}_to_ref.tsv"
    output:
        dotplot="results/rbh/{sample}_{assembly_id}_dotplot.png"
    log:
        "results/logs/slurm/rule_plot_raw_alignments_{sample}_{assembly_id}.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        python workflow/scripts/plot_alignments.py -i {input.raw_alignments} -o {output.dotplot} > {log} 2>&1
        """
