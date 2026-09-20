# workflow/rules/01_pre_assembly_qc.smk

rule filter_hifi_adapters:
    """
    Screens and removes residual SMRTbell adapters using HiFiAdapterFilt.
    """
    input:
        fastq=lambda wildcards: samples_df.loc[wildcards.sample, "file_path"]
    output:
        filt_fastq="results/00_qc/reads/{sample}.filt.fastq.gz",
        blocklist="results/00_qc/reads/{sample}.contaminant.blocklist",
        stats="results/00_qc/reads/{sample}.adapter_stats.txt"
    log:
        "slurm/logs/read_qc/adapter_filter_{sample}.log"
    threads: 8
    resources:
        mem_mb=16000,
        time="01:00:00"
    conda:
        "../../envs/read_qc.yaml"
    params:
        out_dir="results/00_qc/reads"
    shell:
        """
        mkdir -p {params.out_dir} slurm/logs/read_qc
        pbadapterfilt.sh -p {params.out_dir}/{wildcards.sample} -t {threads} > {log} 2>&1
        
        gzip -c {params.out_dir}/{wildcards.sample}.filt.fastq > {output.filt_fastq}
        rm -f {params.out_dir}/{wildcards.sample}.filt.fastq
        
        mv {params.out_dir}/{wildcards.sample}.contaminant.blocklist {output.blocklist} 2>/dev/null || touch {output.blocklist}
        mv {params.out_dir}/{wildcards.sample}.stats {output.stats} 2>/dev/null || touch {output.stats}
        """

rule run_read_qc_nanoplot:
    """
    Profiles read length distribution, quality metrics (Q-scores), and yield.
    """
    input:
        fastq="results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        report_html="results/00_qc/reads/nanoplot_{sample}/NanoPlot-report.html",
        stats_txt="results/00_qc/reads/nanoplot_{sample}/NanoStats.txt"
    log:
        "slurm/logs/read_qc/nanoplot_{sample}.log"
    threads: 8
    resources:
        mem_mb=24000,
        time="01:00:00"
    conda:
        "../../envs/read_qc.yaml"
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
    """
    Calculates empirical data yield and compares against expected genome size.
    """
    input:
        fastq="results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        tsv="results/00_qc/reads/{sample}_theoretical_coverage.tsv",
        json="results/00_qc/reads/{sample}_theoretical_coverage.json"
    log:
        "slurm/logs/read_qc/coverage_{sample}.log"
    threads: 1
    resources:
        mem_mb=8000,
        time="00:30:00"
    conda:
        "../../envs/read_qc.yaml"
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
    """
    Calculates k-mer frequencies (k=21) using FastK.
    """
    input:
        fastq="results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        hist="results/00_qc/kmer/{sample}.hist",
        ktab="results/00_qc/kmer/{sample}.ktab"
    log:
        "slurm/logs/kmer_qc/fastk_{sample}.log"
    threads: 16
    resources:
        mem_mb=64000,
        time="01:30:00"
    conda:
        "../../envs/kmer_qc.yaml"
    params:
        kmer_len=config["kmer_profiling"]["kmer_length"],
        work_dir="results/00_qc/kmer",
        base_name="{sample}"
    shell:
        """
        mkdir -p {params.work_dir} slurm/logs/kmer_qc
        ROOT_DIR=$(pwd)
        
        cd {params.work_dir}
        FastK \
            -k{params.kmer_len} \
            -t{threads} \
            -N{params.base_name} \
            -p $ROOT_DIR/{input.fastq} > $ROOT_DIR/{log} 2>&1
            
        Histex -G {params.base_name} > {wildcards.sample}.hist
        cd $ROOT_DIR
        """

rule run_genomescope:
    """
    Fits mathematical models to k-mer spectra to determine genome size,
    heterozygosity rate, and repeat volume.
    """
    input:
        hist="results/00_qc/kmer/{sample}.hist"
    output:
        model_tsv="results/00_qc/kmer/genomescope_{sample}/model.txt",
        summary_txt="results/00_qc/kmer/genomescope_{sample}/summary.txt",
        linear_plot="results/00_qc/kmer/genomescope_{sample}/linear_plot.png",
        log_plot="results/00_qc/kmer/genomescope_{sample}/log_plot.png"
    log:
        "slurm/logs/kmer_qc/genomescope_{sample}.log"
    threads: 4
    resources:
        mem_mb=16000,
        time="00:30:00"
    conda:
        "../../envs/kmer_qc.yaml"
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

rule run_smudgeplot:
    """
    Extracts k-mer pairs to deconvolve organismal ploidy (diploid vs triploid).
    """
    input:
        hist="results/00_qc/kmer/{sample}.hist"
    output:
        smudge_png="results/00_qc/kmer/smudgeplot_{sample}/{sample}_smudgeplot.png",
        smudge_tsv="results/00_qc/kmer/smudgeplot_{sample}/{sample}_summary.tsv"
    log:
        "slurm/logs/kmer_qc/smudgeplot_{sample}.log"
    threads: 8
    resources:
        mem_mb=32000,
        time="01:00:00"
    conda:
        "../../envs/kmer_qc.yaml"
    params:
        kmer_dir="results/00_qc/kmer",
        out_prefix="results/00_qc/kmer/smudgeplot_{sample}/{sample}"
    shell:
        """
        mkdir -p results/00_qc/kmer/smudgeplot_{wildcards.sample}
        smudgeplot.py cutoff {input.hist} > {params.out_prefix}_cutoffs.txt 2> {log}
        smudgeplot.py plot \
            -o {params.out_prefix} \
            {params.kmer_dir}/{wildcards.sample}.hist > {log} 2>&1 || touch {output.smudge_png} {output.smudge_tsv}
        """

rule pre_assembly_qc_gate:
    """
    Aggregates metrics and evaluates whether the sequencing run meets 
    the minimum criteria to justify launching hifiasm.
    """
    input:
        cov_json="results/00_qc/reads/{sample}_theoretical_coverage.json",
        gscope_summary="results/00_qc/kmer/genomescope_{sample}/summary.txt",
        nanostats="results/00_qc/reads/nanoplot_{sample}/NanoStats.txt",
        smudge_png="results/00_qc/kmer/smudgeplot_{sample}/{sample}_smudgeplot.png"
    output:
        gate_report="results/00_qc/{sample}_pre_assembly_gate.tsv"
    log:
        "slurm/logs/read_qc/gate_{sample}.log"
    threads: 1
    resources:
        mem_mb=4000,
        time="00:15:00"
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
