import pandas as pd

configfile: "config.yaml"
samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample_id", drop=False)

PACBIO_SAMPLES = samples_df[samples_df['data_type'] == 'pacbio']['sample_id'].tolist()
RNA_READS = samples_df[samples_df['data_type'] == 'rnaseq']['file_path'].iloc[0]

# --- CHECKPOINT TARGETS ---

# Target 1: Stop after assembly and QC
rule run_assembly:
    input:
        expand("results/01_assembly/{sample}.fasta", sample=PACBIO_SAMPLES),
        expand("results/01_assembly/QC/quast_{sample}/report.txt", sample=PACBIO_SAMPLES),
        expand("results/01_assembly/QC/busco_{sample}/short_summary.txt", sample=PACBIO_SAMPLES)

# Target 2: Stop after RNA mapping and Gene Prediction
rule run_annotation:
    input:
        expand("results/02_annotation/{sample}_braker/braker.gff3", sample=PACBIO_SAMPLES)

# Target 3: The full pipeline (if you ever want to run it all at once)
rule all:
    input:
        rules.run_assembly.input,
        rules.run_annotation.input
# --- Phase 1: Assembly & QC ---
def get_pacbio_reads(wildcards):
    return samples_df.loc[wildcards.sample, "file_path"]

rule run_hifiasm:
    input:
        reads=get_pacbio_reads
    output:
        gfa="results/01_assembly/{sample}.bp.p_ctg.gfa",
        fasta="results/01_assembly/{sample}.fasta"
    threads:
        config["assembly_threads"]
    conda:
        "envs/assembly.yaml"
    shell:
        """
        hifiasm -o results/01_assembly/{wildcards.sample} -t {threads} {input.reads}
        awk '/^S/{{print ">"$2"\\n"$3}}' {output.gfa} > {output.fasta}
        """

rule run_quast:
    input:
        "results/01_assembly/{sample}.fasta"
    output:
        "results/01_assembly/QC/quast_{sample}/report.txt"
    conda:
        "envs/qc.yaml"
    shell:
        "quast.py {input} -o results/01_assembly/QC/quast_{wildcards.sample}"

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
        "busco -i {input} -o busco_{wildcards.sample} --out_path results/01_assembly/QC -l {params.lineage} -m genome --force"

# --- Phase 2: Annotation ---
rule map_rnaseq:
    input:
        genome="results/01_assembly/{sample}.fasta",
        reads=RNA_READS
    output:
        bam="results/02_annotation/{sample}_rna_mapped.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "minimap2 -ax splice {input.genome} {input.reads} | samtools sort -o {output.bam}"

# NEW RULE: BRAKER3
rule run_braker:
    input:
        genome="results/01_assembly/{sample}.fasta",
        bam="results/02_annotation/{sample}_rna_mapped.bam"
    output:
        # BRAKER creates a directory of files. We specify the main output file here.
        gff3="results/02_annotation/{sample}_braker/braker.gff3"
    params:
        # We pull the species name from config.yaml so we can reuse this pipeline easily
        species=config["species_name"]
    threads:
        8
    conda:
        "envs/annotation.yaml"
    shell:
        """
        # Run the BRAKER3 pipeline
        braker.pl \
            --genome={input.genome} \
            --bam={input.bam} \
            --species={params.species}_{wildcards.sample} \
            --workingdir=results/02_annotation/{wildcards.sample}_braker \
            --threads={threads}
        """
