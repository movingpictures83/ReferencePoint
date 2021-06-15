# ReferencePoint
# Language: Python
# Input: TXT (keyword-value pairs)
# Output: CSV
# Tested with: PluMA 1.1, Python 3.6
# Dependency:

PluMA plugin that takes as input a TXT file of keyword value pairs.

Keywords:
csvfile: A CSV file with genes as rows, assumed to have columns for start and end.
startpos: New starting position for refgene
refgene: Name of reference gene
refstrand: Direction of reference strand (-=complement, +=original)
genomestart: Starting spot of overall gemone
genomeend: Ending spot of overall genome

Output file will be the same as the original CSV file, with indices updated according to the reference
