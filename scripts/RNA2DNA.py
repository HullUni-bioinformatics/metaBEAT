#! /usr/bin/env python

print("\n"+ 'RNA2DNA.py converts an RNA fasta file to uppercase DNA fasta:-')

# Set the input file name
# Run from within the directory containing the data file
InFileName = 'rna.fas' #change this to the name of your input file
 
# Open the input file for reading 'r'
InFile = open(InFileName, 'r')
 
# Initialize a counter used to keep track of how many records are processed
FASTA_records_counter = 0
 
# Set output filename and open this output file for writing 'w'
OutFileName=InFileName + "_dna.fas"
OutFile=open(OutFileName,'w') # You can append to the file instead with 'a'
 
# Loop through each line in the file, do nothing with header lines, make changes to sequence
for line in InFile:
        line = line.strip (); # this tidies up by removing whitespace
        if line.startswith('>'): # ie header lines, they start with the greater-than symbol >
            FASTA_records_counter += 1  #<-- short way of writing 'add 1 to value of FASTA_records_counter'
            output_line = line
        else: # if not a header line its sequence
            output_line = (line.upper()).replace('U','T') #make everything UPPERCASE and replace all U with T
 
        OutFile.write(output_line +"\n") #write each line to the OutFile and add newline \n at the end
 
 
# After the loop is completed, report the actions, and close the files
print('Data was read from file: ' +InFileName)
print(str(FASTA_records_counter) +' records were processed') #str makes this value a 'string' that can be printed
print('Data was written to file: ' +OutFileName +"\n")
InFile.close()
OutFile.close()