import pandas as pd
import re
import os

# Set this as directory for BLASTp NCycDB outputs
directory = '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/manuscript_test/metagenome_analysis/sunbeam_output/annotation/blastp/ncyc/prodigal'

def parse(a):
  # Split commas
  sampleID, chunk, aln_len, e_val = a.split(',')
  # Parse second entry
  ID = chunk.split(' ')[0]
  # Parse last
  gene, ontology, source = re.findall('=(.\w+)', chunk)
  data = dict(zip(['sampleID', 'ID', 'gene', 'ontology', 'source', 'aln_len', 'e_val'], 
                  [sampleID, ID, gene, ontology, source, int(aln_len), float(e_val)]))
  return data
  
data = []

# Loop through all output files and append to one CSV
for entry in os.scandir(directory):
  try:
    print(entry.path)
    file = os.path.basename(entry.path)
    for line in open(file, 'r').readlines():
      try:
        d = parse(line)
        data.append(d)
      except:
        continue

df = pd.DataFrame(data)
df.head()

# Specify output location
df.to_csv('/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/manuscript_test/metagenome_analysis/sunbeam_output/annotation/blastp/blastp_ncyc_allsamples.csv', index = False, header=True)
