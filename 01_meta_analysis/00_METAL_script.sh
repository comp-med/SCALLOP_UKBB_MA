SEPARATOR  TAB
SCHEME STDERR

MARKER  cpaid
ALLELE  EFFECT_ALLELE OTHER_ALLELE
FREQ EAF
EFFECT  BETA
STDERR  SE

AVERAGEFREQ ON
MINMAXFREQ ON

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

ADDFILTER INFO > 0.8
ADDFILTER EAF >= 0.001
ADDFILTER EAF <= 0.999

#Add all the relevant proteins for each protein target
PROCESS file1.txt.gz
PROCESS file2.txt.gz
PROCESS file3.txt.gz
PROCESS file4.txt.gz
PROCESS file5.txt.gz
PROCESS file6.txt.gz
PROCESS file7.txt.gz
PROCESS file8.txt.gz
PROCESS file9.txt.gz
PROCESS file10.txt.gz
PROCESS file11.txt.gz
PROCESS file12.txt.gz
PROCESS file13.txt.gz
PROCESS file14.txt.gz
PROCESS file15.txt.gz
PROCESS file16.txt.gz
PROCESS file17.txt.gz
PROCESS file18.txt.gz
PROCESS file19.txt.gz
PROCESS file20.txt.gz
PROCESS file21.txt.gz
PROCESS file22.txt.gz
PROCESS file23.txt.gz
PROCESS file24.txt.gz
PROCESS file25.txt.gz
PROCESS file26.txt.gz
PROCESS file27.txt.gz
PROCESS file28.txt.gz
PROCESS file29.txt.gz
PROCESS file30.txt.gz

OUTFILE output_SCALLOPMA .tbl
ANALYZE RANDOM

QUIT
