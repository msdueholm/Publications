# First run subset.sh to subset the raw data set to decrease run time and to rename sequence so they are compatibel with medaka 
# Then run canu to make assembly 
# Then run medaka to polish the assemblies 
# then do orientation 
# and in the end make an alignment 
# remove primers and plasmid sequences using trimming based on the alignment. Trimming was done manually in CLC