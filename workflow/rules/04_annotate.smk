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
