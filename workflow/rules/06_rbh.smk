rule download_reference:
    output:
        "data/p_infestans_ref/p_infestans_ref.aa"
    shell:
        """
        mkdir -p data/p_infestans_ref
        wget -O data/p_infestans_ref/ref.faa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/142/945/GCF_000142945.1_ASM14294v1/GCF_000142945.1_ASM14294v1_protein.faa.gz"
        gunzip -c data/p_infestans_ref/ref.faa.gz > {output}
        rm -f data/p_infestans_ref/ref.faa.gz
        """

rule extract_braker_proteins:
    input:
        gtf = "results/annotation/{sample}_{assembly_id}/braker/braker.gtf",
        genome = "results/scaffolds/{sample}_{assembly_id}/ragtag.scaffold.fasta"
    output:
        aa = "results/annotation/{sample}_{assembly_id}/braker/braker.aa"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        gffread {input.gtf} -g {input.genome} -y {output.aa}
        """

rule filter_isoforms:
    input:
        braker_aa = "results/annotation/{sample}_{assembly_id}/braker/braker.aa",
        ref_aa = "data/p_infestans_ref/p_infestans_ref.aa"
    output:
        braker_longest = "results/rbh/{sample}_{assembly_id}/braker_longest.aa",
        ref_longest = "results/rbh/{sample}_{assembly_id}/ref_longest.aa"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        """
        python workflow/scripts/filter_isoforms.py {input.braker_aa} {output.braker_longest} braker
        python workflow/scripts/filter_isoforms.py {input.ref_aa} {output.ref_longest} ncbi
        """

rule run_braker_ref_blast:
    input:
        braker_longest = "results/rbh/{sample}_{assembly_id}/braker_longest.aa",
        ref_longest = "results/rbh/{sample}_{assembly_id}/ref_longest.aa"
    output:
        b2r = "results/rbh/{sample}_{assembly_id}/braker_to_ref.tsv",
        r2b = "results/rbh/{sample}_{assembly_id}/ref_to_braker.tsv"
    conda:
        "../envs/rbh_tools.yaml"
    threads: 48
    shell:
        """
        DB_DIR="results/rbh/{wildcards.sample}_{wildcards.assembly_id}"
        
        makeblastdb -in {input.ref_longest} -dbtype prot -out $DB_DIR/ref_db
        makeblastdb -in {input.braker_longest} -dbtype prot -out $DB_DIR/braker_db

        blastp -query {input.braker_longest} -db $DB_DIR/ref_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 -num_threads {threads} -out {output.b2r}

        blastp -query {input.ref_longest} -db $DB_DIR/braker_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 -num_threads {threads} -out {output.r2b}
        """

rule analyze_stringency:
    input:
        b2r = "results/rbh/{sample}_{assembly_id}/braker_to_ref.tsv",
        r2b = "results/rbh/{sample}_{assembly_id}/ref_to_braker.tsv",
        braker_longest = "results/rbh/{sample}_{assembly_id}/braker_longest.aa"
    output:
        report = "results/rbh/{sample}_{assembly_id}/rbh_stringency_report.txt"
    conda:
        "../envs/rbh_tools.yaml"
    shell:
        "python workflow/scripts/analyze_rbh_stringency.py {input.b2r} {input.r2b} {input.braker_longest} {output.report}"
