# Script to take BLASTxml file outputs, parse hits, and save to CSV for visualization
from Bio.Blast import NCBIXML 
import sys
import csv
import re
import os

E_VALUE_THRESH=0.03
# Set this as directory for BLASTp NCycDB outputs
directory = '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/sunbeam_output/annotation/blastp/ncyc/prodigal'
# Specify output location
with open('/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/sunbeam_output/BLAST/blastp_ncyc.csv','w') as f:
  output = csv.writer(f, delimiter=',')
# Loop through all output files and append to one CSV
  for entry in os.scandir(directory):
    try:	
      print(entry.path)
      base_path = os.path.basename(entry.path)
      result_handle = open(entry.path)
      blast_records = NCBIXML.parse(result_handle)
      blast_records = list(blast_records)
      for blast_record in blast_records:
        for alignment in blast_record.alignments:
          for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
              row = [base_path, alignment.title, alignment.length, hsp.expect]
              output.writerow(row)
    except:
      pass
        
