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
        expand("results/03_evaluate/blast/Calibration_Isolate_{assembly_id}/Calibration_Isolate_{assembly_id}_dotplot.png",
               assembly_id=list(config["active_assemblies"].keys())),

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
            "results/annotation/{sample}_{assembly_id}/braker/braker.codingseq",
            sample="Calibration_Isolate",
            assembly_id="l3_default"
        ),
        expand(
            "results/annotation/{sample}_{assembly_id}/braker/braker.gtf",
            sample="Calibration_Isolate",
            assembly_id="l3_default"
        ),
        # 6. Reciprocal Best Hits Stringency Report
        expand(
            "results/rbh/{sample}_{assembly_id}/rbh_stringency_report.txt",
            sample="Calibration_Isolate",
            assembly_id="l3_default"
        )
# 3. MODULE INCLUDES
include: "workflow/rules/01_pre_assembly_qc.smk"
include: "workflow/rules/02_assemble.smk"       # Contains hifiasm code
include: "workflow/rules/03_evaluate.smk"       # Contains QUAST/BUSCO/Merqury & BLAST/Dotplot
include: "workflow/rules/04_curate_and_scaffold.smk"
include: "workflow/rules/05_annotate.smk"       # Contains BRAKER rules
include: "workflow/rules/06_rbh.smk"
# 4. PIPELINE COMPLETION HOOKS
onsuccess:
    print("Pipeline finished successfully. Generating PI job ledger...")
    shell("python3 workflow/scripts/generate_pi_ledger.py")

onerror:
    print("Pipeline crashed. Generating PI job ledger to check for failures...")
    shell("python3 workflow/scripts/generate_pi_ledger.py")
