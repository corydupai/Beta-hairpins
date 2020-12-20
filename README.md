# Beta hairpins
Code to identify/analyze beta hairpin motifs in pdb files

## Key scripts
* scripts/parse_secondary_structures.R finds motifs in pdb files
* main.Rmd does most of the heavy lifting on analysis and figure generation
* logo_generation.Rmd makes the motiff logo for supplemental figure
* make_dssp.sh converts pdb .cif.gz files to .dssp like files with secondary structure info but no headers

## Data files
* data/beta_hairpin_pairs.csv has info on contacting residues between beta strands of analyzed sequences
* data/beta_hairpin_sequences.csv has info on included sequences
* data_mid/pair_data_bigger.csv.gz has all of the pairs output by scripts/parse_secondary_structures.R. Unzip this to read into main.Rmd