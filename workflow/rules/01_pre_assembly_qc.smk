rule filter_hifi_adapters:
    input:
        fastq=lambda wildcards: samples_df.loc[wildcards.sample, "file_path"]
    output:
        filt_fastq="results/00_qc/reads/{sample}.filt.fastq.gz",
        blocklist="results/00_qc/reads/{sample}.contaminant.blocklist",
        stats="results/00_qc/reads/{sample}.adapter_stats.txt"
    log:
        "results/logs/slurm/rule_filter_hifi_adapters_{sample}.log"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=1440
    conda:
        "../envs/read_qc.yaml"
    params:
        db="workflow/scripts/HiFiAdapterFilt/DB/pacbio_vectors_db"
    shell:
        """
        set -euo pipefail
        mkdir -p results/00_qc/reads

        echo "Starting chunked native adapter filtration on local SSD..." > {log}

        LOCAL_FASTA="$TMPDIR/{wildcards.sample}_temp.fasta"
        LOCAL_BLAST="$TMPDIR/{wildcards.sample}_blast.tmp"
        LOCAL_CHUNKS="$TMPDIR/{wildcards.sample}_chunks"

        # 1. Convert FASTQ to FASTA to SSD
        seqkit fq2fa {input.fastq} > $LOCAL_FASTA 2>> {log}

        # 2. Slice the 11GB file into 20 safe, bite-sized pieces
        mkdir -p $LOCAL_CHUNKS
        seqkit split2 -p 20 $LOCAL_FASTA -O $LOCAL_CHUNKS 2>> {log}

        # 3. Loop BLAST over each chunk independently, appending (>>) to the master tmp file
        touch $LOCAL_BLAST
        for chunk in $LOCAL_CHUNKS/*.fasta*; do
            echo "Processing $chunk..." >> {log}
            blastn -query "$chunk" -db {params.db} \
                -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 \
                -dust no -soft_masking true -evalue 0.1 -num_threads {threads} \
                -outfmt 6 >> $LOCAL_BLAST 2>> {log}
        done

        # 4. Extract IDs to Quobyte
        cut -f1 $LOCAL_BLAST | sort | uniq > {output.blocklist}

        # 5. Filter original reads
        if [ -s {output.blocklist} ]; then
            seqkit grep -v -f {output.blocklist} {input.fastq} | gzip -c > {output.filt_fastq}
        else
            gzip -c {input.fastq} > {output.filt_fastq}
        fi

        # 6. Stats & Cleanup
        BAD_COUNT=$(wc -l < {output.blocklist})
        echo "SMRTbell Adapter-Contaminated Reads Removed: $BAD_COUNT" > {output.stats}

        rm -rf $LOCAL_BLAST $LOCAL_FASTA $LOCAL_CHUNKS
        echo "Filtration complete." >> {log}
        """
rule run_read_qc_nanoplot:
    input:
        fastq="results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        report_html="results/00_qc/reads/nanoplot_{sample}/NanoPlot-report.html",
        stats_txt="results/00_qc/reads/nanoplot_{sample}/NanoStats.txt"
    log:
        "results/logs/slurm/rule_nanoplot_{sample}.log"
    threads: 8
    resources:
        mem_mb=24000,
        runtime=360
    conda:
        "../envs/read_qc.yaml"
    params:
        outdir="results/00_qc/reads/nanoplot_{sample}"
    shell:
        """
        mkdir -p {params.outdir}
        NanoPlot \
            --fastq {input.fastq} \
            --outdir {params.outdir} \
            --threads {threads} \
            --plots hex dot \
            --N50 \
            --title "{wildcards.sample} HiFi Read QC" > {log} 2>&1
        """

rule calculate_theoretical_coverage:
    input:
        fastq="results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        tsv="results/00_qc/reads/{sample}_theoretical_coverage.tsv",
        json="results/00_qc/reads/{sample}_theoretical_coverage.json"
    log:
        "results/logs/slurm/rule_coverage_{sample}.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/read_qc.yaml"
    params:
        expected_size=config["organism"]["expected_genome_size_bp"]
    shell:
        """
        python workflow/scripts/calculate_theoretical_coverage.py \
            --fastq {input.fastq} \
            --expected-size {params.expected_size} \
            --sample-id {wildcards.sample} \
            --output-tsv {output.tsv} \
            --output-json {output.json} > {log} 2>&1
        """

rule count_kmers_fastk:
    input:
        "results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        hist="results/00_qc/kmer/{sample}.hist",
        ktab="results/00_qc/kmer/{sample}.ktab"
    log:
        "results/logs/slurm/rule_fastk_{sample}.log"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60
    conda:
        "../envs/kmer_qc.yaml"
    shell:
        """
        mkdir -p results/00_qc/kmer
        ROOT_DIR=$(pwd)
        
        cd results/00_qc/kmer
        FastK -k21 -t{threads} -N{wildcards.sample} -p $ROOT_DIR/{input} > $ROOT_DIR/{log} 2>&1
        
        Histex -G {wildcards.sample} > {wildcards.sample}.hist.tmp
        mv {wildcards.sample}.hist.tmp {wildcards.sample}.hist
        
        cd $ROOT_DIR
        """

rule run_genomescope:
    input:
        hist="results/00_qc/kmer/{sample}.hist"
    output:
        model_tsv="results/00_qc/kmer/genomescope_{sample}/model.txt",
        summary_txt="results/00_qc/kmer/genomescope_{sample}/summary.txt",
        linear_plot="results/00_qc/kmer/genomescope_{sample}/linear_plot.png",
        log_plot="results/00_qc/kmer/genomescope_{sample}/log_plot.png"
    log:
        "results/logs/slurm/rule_genomescope_{sample}.log"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60
    conda:
        "../envs/kmer_qc.yaml"
    params:
        kmer_len=config["kmer_profiling"]["kmer_length"],
        ploidy=config["kmer_profiling"]["ploidy"],
        out_dir="results/00_qc/kmer/genomescope_{sample}"
    shell:
        """
        mkdir -p {params.out_dir}
        genomescope2 \
            -i {input.hist} \
            -o {params.out_dir} \
            -k {params.kmer_len} \
            -p {params.ploidy} > {log} 2>&1
        """

rule pre_assembly_qc_gate:
    input:
        cov_json="results/00_qc/reads/{sample}_theoretical_coverage.json",
        gscope_summary="results/00_qc/kmer/genomescope_{sample}/summary.txt",
        nanostats="results/00_qc/reads/nanoplot_{sample}/NanoStats.txt",
    output:
        gate_report="results/00_qc/{sample}_pre_assembly_gate.tsv"
    log:
        "results/logs/slurm/rule_gate_{sample}.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15
    params:
        min_cov=config["organism"]["expected_min_coverage"],
        max_cov=config["organism"]["expected_max_coverage"]
    run:
        import json

        with open(input.cov_json) as f:
            cov_data = json.load(f)

        actual_cov = cov_data["theoretical_coverage"]
        n50_len = cov_data["n50_read_length"]

        status = "PASSED"
        reasons = []

        if actual_cov < params.min_cov:
            status = "FAILED"
            reasons.append(f"Insufficient coverage ({actual_cov:.1f}X < {params.min_cov}X).")
        elif actual_cov > params.max_cov:
            status = "WARNING"
            reasons.append(f"Excessive coverage ({actual_cov:.1f}X > {params.max_cov}X); consider downsampling to save memory.")

        if n50_len < 8000:
            status = "WARNING"
            reasons.append(f"Short N50 read length ({n50_len} bp < 8000 bp); transposable elements may fragment assembly.")

        with open(output.gate_report, "w") as out:
            out.write("sample_id\tstatus\ttheoretical_coverage\tn50_read_length\tnotes\n")
            out.write(f"{wildcards.sample}\t{status}\t{actual_cov}\t{n50_len}\t{'; '.join(reasons) if reasons else 'Dataset cleared for assembly.'}\n")

        if status == "FAILED":
            raise ValueError(f"QC GATE FAILURE for {wildcards.sample}: {'; '.join(reasons)}")
