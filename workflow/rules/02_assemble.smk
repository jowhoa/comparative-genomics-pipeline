import os

# 1. Target rule for the assembly phase
# Snakemake looks at the config to see which assemblies are active,
# then requests the final FASTA files for each one.
rule all_assemblies:
    input:
        expand(
            "results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta",
            sample="Calibration_Isolate",
            assembly_id=list(config["active_assemblies"].keys())
        )

# 2. The core hifiasm execution rule
rule run_hifiasm:
    """
    Executes hifiasm using dynamic parameters pulled from config.yaml.
    """
    input:
        reads="results/00_qc/reads/{sample}.filt.fastq.gz" # Adjust this path if your QC'd reads live elsewhere
    output:
        # We target the primary contig GFA. hifiasm generates several other files, 
        # but Snakemake only strictly needs to track the ones we use downstream.
        primary_gfa="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.gfa",
        alternate_gfa="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.a_ctg.gfa"
    log:
        "slurm/logs/assembly/hifiasm_{sample}_{assembly_id}.log"
    threads: 
        config["assembly_threads"]
    resources:
        mem_mb=120000, # 120GB RAM is standard for a 240Mb genome
        time="24:00:00"
    params:
        # MAGIC HAPPENS HERE:
        # The lambda function looks at the {assembly_id} wildcard (e.g., "l3_default"),
        # goes into the config dictionary, and pulls out the string "-l3".
        flags=lambda wildcards: config["active_assemblies"][wildcards.assembly_id],
        prefix="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}"
    shell:
        """
        # Run hifiasm with the injected flags
        hifiasm {params.flags} -o {params.prefix} -t {threads} {input.reads} 2> {log}
        """

# 3. Graph-to-FASTA conversion
rule gfa2fasta:
    """
    Downstream tools (BUSCO, QUAST) require FASTA files, not GFA graphs.
    This uses awk to extract the sequence lines from the GFA.
    """
    input:
        gfa="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.gfa"
    output:
        fasta="results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta"
    log:
        "slurm/logs/assembly/gfa2fasta_{sample}_{assembly_id}.log"
    threads: 1
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        # Extracts the 'S' (Segment) lines from the GFA and formats them as standard FASTA
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.gfa} > {output.fasta} 2> {log}
        """
