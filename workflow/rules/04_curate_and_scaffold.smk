rule run_kraken2:
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta"
    output:
        kraken_out="results/qc/kraken/{sample}_{assembly_id}.kraken",
        kraken_report="results/qc/kraken/{sample}_{assembly_id}.report"
    params:
        db=config["kraken_db"]
    threads: 24
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 --db {params.db} \
            --threads {threads} \
            --report {output.kraken_report} \
            --output {output.kraken_out} \
            {input.fasta}
        """

rule extract_clean_contigs:
    input:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta",
        kraken_out="results/qc/kraken/{sample}_{assembly_id}.kraken",
        kraken_report="results/qc/kraken/{sample}_{assembly_id}.report"
    output:
        clean_fasta="results/assembly_clean/{sample}_{assembly_id}_clean.fasta"
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
            --include-children
        """
rule run_ragtag:
    input:
        clean_fasta="results/assembly_clean/{sample}_{assembly_id}_clean.fasta",
        reference=config["reference_genome"]
    output:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    threads: 8
    conda:
        "../envs/scaffold.yaml"
    shell:
        """
        # RagTag automatically creates its own output directory. 
        # We point it to the results/scaffolds folder.
        ragtag.py scaffold {input.reference} {input.clean_fasta} \
            -o results/scaffolds/{wildcards.sample}_{wildcards.assembly_id} \
            -t {threads}
        """

rule find_telomeres:
    input:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    output:
        telomere_plot="results/telomeres/{sample}_{assembly_id}/{sample}_{assembly_id}_telomeric_repeat_windows.tsv"
    threads: 4
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
            {input.scaffolds}
        """

rule make_blast_db_reference:
    input:
        reference=config["reference_genome"]
    output:
        ref_db="results/rbh/dbs/reference/reference.ndb"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        makeblastdb -in {input.reference} -dbtype nucl -out results/rbh/dbs/reference/reference
        """

rule make_blast_db_scaffolds:
    input:
        scaffolds="results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    output:
        scaffold_db="results/rbh/dbs/{sample}_{assembly_id}/{sample}_{assembly_id}.ndb"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        makeblastdb -in {input.scaffolds} -dbtype nucl -out results/rbh/dbs/{wildcards.sample}_{wildcards.assembly_id}/{wildcards.sample}_{wildcards.assembly_id}
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
    threads: 16
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        blastn -query {input.scaffolds} -db results/rbh/dbs/reference/reference \
            -outfmt 6 -num_threads {threads} -out {output.q_to_r}
            
        blastn -query {input.reference} -db results/rbh/dbs/{wildcards.sample}_{wildcards.assembly_id}/{wildcards.sample}_{wildcards.assembly_id} \
            -outfmt 6 -num_threads {threads} -out {output.r_to_q}
        """

# 1. Generate the strict 1:1 ortholog table (No plotting)
rule parse_rbh:
    input:
        q_to_r="results/rbh/blast/{sample}_{assembly_id}_to_ref.tsv",
        r_to_q="results/rbh/blast/ref_to_{sample}_{assembly_id}.tsv"
    output:
        rbh_table="results/rbh/{sample}_{assembly_id}_reciprocal_best_hits.tsv"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        "python workflow/scripts/find_rbh.py -q {input.q_to_r} -r {input.r_to_q} -o {output.rbh_table}"

# 2. Plot the whole-genome alignments using the raw data
rule plot_raw_alignments:
    input:
        # Feeding the raw, unfiltered BLAST alignments to the script instead of the RBH table
        raw_alignments="results/rbh/blast/{sample}_{assembly_id}_to_ref.tsv"
    output:
        dotplot="results/rbh/{sample}_{assembly_id}_dotplot.png"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        "python workflow/scripts/plot_alignments.py -i {input.raw_alignments} -o {output.dotplot}"
