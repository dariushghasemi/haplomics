
import os
from pathlib import Path
import pandas as pd

path_loci = config.get("path_loci")

# Define input for the rules
data = []
with open(path_loci, "r") as fp:
    lines = fp.readlines()

for line in lines:
    # Split the line by tab and get the relevant parts
    parts = line.strip().split("\t")
    
    # Assuming the locus is the fourth element (index 3) and coordinates are the first three elements
    locus = parts[3].strip()  # Remove any extra whitespace or newline characters
    coordinates = parts[:3]    # First three columns represent the coordinates
    
    # Filter out the line with 'LOCUS'
    if locus != "LOCUS":
        data.append((locus, coordinates))





# def parse_loci(loci_path):
#     """Parse the input config file to extract loci."""
#     loci = []

#     with open(loci_path) as f:
#         next(f)  # Skip header
#         for line in f:
#             parts = line.strip().split("\t")
#             locus = parts[3].strip()
            
#             loci.append(locus)
    
#     return loci


analytes = (
    pd.DataFrame.from_records(data, columns=["locus", "coordinates"])
    .set_index("locus", drop=False)
    .sort_index()
)

# config.py
def parse_phenotype(data):
    """Parse the phenotype file paths for each data type."""
    pheno_paths = {}

    with open(data) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            data_type = parts[1]  # phen, meta, or prot
            path = parts[2]       # corresponding file path
            pheno_paths[data_type] = path
    
    return pheno_paths


# controlling phenotypes list
data_types = config.get("datasets").split(",")
phenotype_files = config.get("phenotype_files").split(",")

# Split the strings to get individual elements
file_paths = [path.strip() for path in phenotype_files]

# Create the DataFrame
df = (
    pd.DataFrame({"data_type": data_types, "data_path": file_paths})
    .set_index("data_type", drop=False)
    .sort_index() )

#print(df)
#print(analytes)

def get_pheno(wildcards):
    return str(Path(df.loc[wildcards, "data_path"]))


def get_region(wildcards):
    chrom, beg, end = analytes.loc[wildcards, "coordinates"]
    # Format as 'chr:beg-end' for bcftools input region
    return f"{chrom}:{beg}-{end}"


# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("path_workspace"), file_path))

# define the functions generating files' path
def full_path(file_path):
   return str(os.path.abspath(Path(config.get("path_workspace"), file_path)))
