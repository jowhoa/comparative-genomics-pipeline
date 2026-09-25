# ==================================================
# MASTER SNAKEFILE
# Phytophthora infestans Comparative Genomics
# ==================================================

import pandas as pd

# 1. LOAD CONFIGURATION
configfile: "config/config.yaml"

samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample_id", drop=False)
# 2. TARGET RULE (The Finish Line)
rule all:
    input:
        # 1. Ensure the QC gate passes for the Calibration Isolate
        "results/00_qc/Calibration_Isolate_pre_assembly_gate.tsv",
        
        # 2. Ensure hifiasm generates the FASTAs for all active parameters
        expand(
            "results/assemblies/{sample}/{assembly_id}/{sample}_{assembly_id}.p_ctg.fasta",
            sample="Calibration_Isolate",
            assembly_id=list(config["active_assemblies"].keys())
        )
        
        # (We will add the QUAST/BUSCO/BRAKER outputs back here in the next step
        # once we verify that the hifiasm matrix successfully launches!)

# 3. MODULE INCLUDES
include: "workflow/rules/01_pre_assembly_qc.smk" 
include: "workflow/rules/02_assemble.smk"       # Contains  hifiasm code
include: "workflow/rules/03_evaluate.smk"       # Contains QUAST/BUSCO rules
#include: "workflow/rules/04_annotate.smk"       # Contains BRAKER rules
