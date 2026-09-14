import pandas as pd

# 1. Load configuration
configfile: "config.yaml"

# 2. Read the sample sheet using Pandas (the industry standard way)
samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample_id", drop=False)

# 3. Target rule (We will add our final expected files here later)
rule all:
    input:
        []
