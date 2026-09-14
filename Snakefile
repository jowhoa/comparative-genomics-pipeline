import pandas as pd

# 1. Load config and sample sheet
configfile: "config.yaml"
samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample_id", drop=False)

PACBIO_SAMPLES = samples_df[samples_df['data_type'] == 'pacbio']['sample_id'].tolist()

# 2. Target rule (UPDATED)
rule all:
    input:
        expand("results/01_assembly/{sample}.fasta", sample=PACBIO_SAMPLES),
        expand("results/01_assembly/QC/quast_{sample}/report.txt", sample=PACBIO_SAMPLES),
        expand("results/01_assembly/QC/busco_{sample}/short_summary.txt", sample=PACBIO_SAMPLES)

# 3. Input lookup function
def get_pacbio_reads(wildcards):
    return samples_df.loc[wildcards.sample, "file_path"]

# 4. Assembly Rule
rule run_hifiasm:
    input:
        reads=get_pacbio_reads
    output:
        gfa="results/01_assembly/{sample}.bp.p_ctg.gfa",
        fasta="results/01_assembly/{sample}.fasta"
    threads:
        config["assembly_threads"]
    log:
        "logs/hifiasm_{sample}.log"
    benchmark:
        "benchmarks/hifiasm_{sample}.txt"
    conda:
        "envs/assembly.yaml"
    shell:
        """
        hifiasm -o results/01_assembly/{wildcards.sample} -t {threads} {input.reads} 2> {log}
        awk '/^S/{{print ">"$2"\\n"$3}}' {output.gfa} > {output.fasta}
        """

# 5. Physical QC Rule
rule run_quast:
    input:
        "results/01_assembly/{sample}.fasta"
    output:
        "results/01_assembly/QC/quast_{sample}/report.txt"
    conda:
        "envs/qc.yaml"
    shell:
        "quast.py {input} -o results/01_assembly/QC/quast_{wildcards.sample}"

# 6. Biological QC Rule
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
        # BUSCO automatically creates the out_path directory if it doesn't exist
        "busco -i {input} -o busco_{wildcards.sample} --out_path results/01_assembly/QC -l {params.lineage} -m genome --force"
