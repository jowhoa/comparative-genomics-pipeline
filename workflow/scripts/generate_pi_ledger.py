import subprocess
import re
import glob

internal_to_rule = {}
external_to_rule = {}

# 1. Parse Snakemake's hidden logs to map SLURM UUIDs back to human-readable rule names
for log in glob.glob(".snakemake/log/*.log"):
    with open(log, "r") as f:
        content = f.read()
        
    # Map internal jobid -> rule name
    blocks = content.split("rule ")[1:]
    for block in blocks:
        lines = block.split("\n")
        rule_name = lines[0].strip().strip(":")
        for line in lines:
            m = re.search(r"jobid:\s*(\d+)", line)
            if m:
                internal_to_rule[m.group(1)] = rule_name
                break
                
    # Map external SLURM id -> internal jobid
    for line in content.split("\n"):
        m = re.search(r"Submitted job (\d+) with external jobid '(\d+)'", line)
        if m:
            internal_id = m.group(1)
            external_id = m.group(2)
            if internal_id in internal_to_rule:
                external_to_rule[external_id] = internal_to_rule[internal_id]

# 2. Pull the permanent accounting data from SLURM for your username
sacct_cmd = ["sacct", "-u", "jowhoa26", "--format=JobID,State,Elapsed,MaxRSS", "-n", "-P"]
result = subprocess.run(sacct_cmd, capture_output=True, text=True)

# 3. Merge the data and write a clean spreadsheet for your PI
with open("PI_Job_Ledger.tsv", "w") as out:
    out.write("SLURM_JobID\tRule_Name\tState\tRuntime\tPeak_Memory\n")
    
    for line in result.stdout.strip().split("\n"):
        if not line: continue
        parts = line.split("|")
        job_id = parts[0].split(".")[0] # Clean out the .batch and .extern trailing IDs
        state = parts[1]
        elapsed = parts[2]
        memory = parts[3] if len(parts) > 3 else ""
        
        # Skip the duplicate .batch/.extern lines to keep the ledger clean
        if "." in parts[0]:
            continue
            
        rule_name = external_to_rule.get(job_id, "Unknown_Rule_or_Old_Job")
        out.write(f"{job_id}\t{rule_name}\t{state}\t{elapsed}\t{memory}\n")

print("Success: Generated PI_Job_Ledger.tsv")
