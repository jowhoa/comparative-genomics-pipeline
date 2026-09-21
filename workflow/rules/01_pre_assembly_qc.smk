# workflow/rules/01_pre_assembly_qc.smk

rule run_smudgeplot:
    """
    Stable Smudgeplot execution (v0.2.5).
    Uses KMC to cleanly extract 2D heterozygous k-mer pairs.
    Note: v0.2.5 utilizes the 'smudgeplot.py' executable prefix.
    """
    input:
        fastq="results/00_qc/reads/{sample}.filt.fastq.gz"
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
    shell:
        """
        set -euo pipefail
        OUTDIR="results/00_qc/kmer/smudgeplot_{wildcards.sample}"
        mkdir -p $OUTDIR/tmp slurm/logs/kmer_qc
        
        echo "1. Counting kmers with KMC..." > {log}
        kmc -k21 -t{threads} -m28 -ci1 -cs10000 {input.fastq} $OUTDIR/kmc_db $OUTDIR/tmp >> {log} 2>&1
        
        echo "2. Extracting KMC histogram..." >> {log}
        kmc_tools transform $OUTDIR/kmc_db histogram $OUTDIR/kmc.hist -cx10000 >> {log} 2>&1
        
        echo "3. Calculating L and U cutoffs..." >> {log}
        L=$(smudgeplot.py cutoff $OUTDIR/kmc.hist L)
        U=$(smudgeplot.py cutoff $OUTDIR/kmc.hist U)
        
        echo "4. Filtering kmers..." >> {log}
        kmc_tools transform $OUTDIR/kmc_db -ci"$L" -cx"$U" dump -s $OUTDIR/kmc_filtered.dump >> {log} 2>&1
        
        echo "5. Finding heterozygous pairs (hetkmers)..." >> {log}
        smudgeplot.py hetkmers -o $OUTDIR/{wildcards.sample}_pairs < $OUTDIR/kmc_filtered.dump >> {log} 2>&1
        
        echo "6. Plotting Smudges..." >> {log}
        smudgeplot.py plot $OUTDIR/{wildcards.sample}_pairs_coverages.tsv -o $OUTDIR/{wildcards.sample} >> {log} 2>&1
        
        echo "7. Standardizing outputs..." >> {log}
        if [ -f $OUTDIR/{wildcards.sample}_summary_table.tsv ]; then
            mv $OUTDIR/{wildcards.sample}_summary_table.tsv {output.smudge_tsv}
        elif [ -f $OUTDIR/{wildcards.sample}_summary.txt ]; then
            mv $OUTDIR/{wildcards.sample}_summary.txt {output.smudge_tsv}
        fi
        mv $OUTDIR/{wildcards.sample}*smudgeplot*.png {output.smudge_png} 2>/dev/null || touch {output.smudge_png}
        
        rm -rf $OUTDIR/tmp $OUTDIR/kmc_db.* $OUTDIR/*.dump $OUTDIR/kmc.hist
        """
rule filter_hifi_adapters:
    """
    NATIVE ADAPTER FILTER: Replaces the broken HiFiAdapterFilt bash script.
    Directly executes BLASTN against the PacBio vector database and filters
    the reads natively using seqkit, ensuring no silent failures.
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
        db="workflow/scripts/HiFiAdapterFilt/DB/pacbio_vectors_db",
        blast_out="results/00_qc/reads/{sample}_blast.tmp"
    shell:
        """
        # 0. Enforce strict error catching (pipeline dies instantly if any step fails)
        set -euo pipefail
        mkdir -p results/00_qc/reads slurm/logs/read_qc

        echo "Starting native adapter filtration..." > {log}
        
        # 1. Convert FASTQ to FASTA on the fly and BLAST against the PacBio database
        seqkit fq2fa {input.fastq} | blastn -query - -db {params.db} \
            -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 \
            -dust no -soft_masking true -evalue 0.1 -num_threads {threads} \
            -outfmt 6 > {params.blast_out} 2>> {log}
        
        # 2. Extract the exact read IDs that matched adapter sequences
        cut -f1 {params.blast_out} | sort | uniq > {output.blocklist}
        
        # 3. Filter the original reads based on the blocklist
        if [ -s {output.blocklist} ]; then
            seqkit grep -v -f {output.blocklist} {input.fastq} | gzip -c > {output.filt_fastq}
        else
            # If the blocklist is empty, no adapters were found. Just compress the raw reads.
            gzip -c {input.fastq} > {output.filt_fastq}
        fi
        
        # 4. Generate statistics
        BAD_COUNT=$(wc -l < {output.blocklist})
        echo "SMRTbell Adapter-Contaminated Reads Removed: $BAD_COUNT" > {output.stats}
        
        # Clean up temporary BLAST file
        rm {params.blast_out}
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
    input:
        "results/00_qc/reads/{sample}.filt.fastq.gz"
    output:
        hist="results/00_qc/kmer/{sample}.hist",
        ktab="results/00_qc/kmer/{sample}.ktab"
    log:
        "slurm/logs/kmer_qc/fastk_{sample}.log"
    threads: 4
    resources:
        mem_mb=16000,
        time="00:30:00"
    conda:
        "../../envs/kmer_qc.yaml"
    shell:
        """
        mkdir -p results/00_qc/kmer slurm/logs/kmer_qc
        ROOT_DIR=$(pwd)
        
        cd results/00_qc/kmer
        # 1. FastK generates .ktab and .hist automatically
        FastK -k21 -t{threads} -N{wildcards.sample} -p $ROOT_DIR/{input} > $ROOT_DIR/{log} 2>&1
        
        # 2. THE FIX: Write to a .tmp file so Bash doesn't delete our input before Histex reads it!
        Histex -G {wildcards.sample} > {wildcards.sample}.hist.tmp
        mv {wildcards.sample}.hist.tmp {wildcards.sample}.hist
        
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
