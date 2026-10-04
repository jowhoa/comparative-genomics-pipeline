import pandas as pd

# 1. LOAD CONFIGURATION
configfile: "config/config.yaml"

samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample_id", drop=False)
wildcard_constraints:
    sample="|".join(samples_df.index)

# 2. TARGET RULE (The Finish Line)
rule all:
    input:
        # 1. Ensure the QC gate passes
        "results/00_qc/Calibration_Isolate_pre_assembly_gate.tsv",

        # 2. Ensure hifiasm generates the FASTAs
        expand(
            "results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta",
            sample="Calibration_Isolate",
            assembly_id=list(config["active_assemblies"].keys())
        ),

        # 3. Ensure Evaluation Reports are generated
        expand("results/03_evaluate/quast/Calibration_Isolate_{assembly_id}/Calibration_Isolate_{assembly_id}_report.html",
               assembly_id=list(config["active_assemblies"].keys())),
        expand("results/03_evaluate/busco/Calibration_Isolate_{assembly_id}/short_summary.txt",
               assembly_id=list(config["active_assemblies"].keys())),
        expand("results/03_evaluate/merqury/Calibration_Isolate_{assembly_id}/Calibration_Isolate_{assembly_id}.qv",
               assembly_id=list(config["active_assemblies"].keys())),
        "results/01_assembly/QC/blast_Calibration_Isolate/Calibration_Isolate_dotplot.png",

        # 4. Ensure Contamination Filtration, Scaffolding, Telomeres, and RBH runs
        expand(
            [  "results/telomeres/{sample}_{assembly_id}/{sample}_{assembly_id}_telomeric_repeat_windows.tsv",
                "results/rbh/{sample}_{assembly_id}_reciprocal_best_hits.tsv",
                "results/rbh/{sample}_{assembly_id}_dotplot.png"
            ],
            sample="Calibration_Isolate",
            assembly_id="l3_default"
        ),

        # 5. Ensure BRAKER3 Gene Annotation runs
        expand(
            "results/annotation/{sample}_{assembly_id}/braker/braker.aa",
            sample="Calibration_Isolate",
            assembly_id="l3_default"
        )

# 3. MODULE INCLUDES
include: "workflow/rules/01_pre_assembly_qc.smk"
include: "workflow/rules/02_assemble.smk"       # Contains hifiasm code
include: "workflow/rules/03_evaluate.smk"       # Contains QUAST/BUSCO/Merqury & BLAST/Dotplot
include: "workflow/rules/04_curate_and_scaffold.smk"
include: "workflow/rules/05_annotate.smk"       # Contains BRAKER rules
