'''
This program takes a 2 fastq file(R1 and R2) and output X number fastq file with reads with a the specific adaptors (In R2), special for 3'
adaptor spliting. Output = name1_R1.fastq, name1_R2.fastq,..... nameN_R1.fastq, nameN_R2.fastq.
'''
#Import systems and functions

import sys,gzip,regex
import pandas as pd

'''
Command in UNIX: python arg1: pipelineLocation arg2: fasqfile1Location arg3: fasqfile2Location arg4: Sample sheet arg5: No of mismaches allowed
Be carefull to put a different name to do not erase previous files
'''

#Generate a list with R1 files locations
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

# Define the path to the file 1, in the argument 2 in UNIX commmand
path1 = (sys.argv[1])

# Define the path to the file 2, in the argument 3 in UNIX commmand
path2 = (sys.argv[2])

#Load sample sheet
df = pd.read_csv(sys.argv[3])
print('SampleSheet Ok!')



# Write a list with the names of samples associated with the adaptors (In proper order)
names =list(df['Sample'])
seq =list(df['Adaptor'])
nadaptor=len(seq)
print(names)
print(seq)
#Define mismaches in the unix comand
mis = (sys.argv[4])
print('Mis allowed= '+mis)

# Create output files
files1 = []
for n in names:
	files1.append( open( n + '_R1.fastq', 'w'))
files2 = []
for n in names:
	files2.append( open( n + '_R2.fastq', 'w'))

#Demultiplexing 3' Barcoded libraries Version 3
def fastq_spliter3(path1, path2,nadaptor):
	file1 = gzip.open(path1,'r')
	file2 = gzip.open(path2,'r')
	Total = 0
	determined=0
	undetermined=0
	Goal=1000000
    # Start reading the file
	a1 = file1.readline()
	while a1 != '':
        # Read File 1
		b1 = file1.readline()
		c1 = file1.readline()
		d1 = file1.readline()
        # Read File2
		a2 = file2.readline()
		b2 = file2.readline()
		c2 = file2.readline()
		d2 = file2.readline()
		barcode = b2[ 0 : 7 ]
		for i in range(nadaptor):
			if regex.match('('+seq[i]+')'+'{s<='+mis+'}',barcode):
				barcode1 = b2[:10]
				out_file = files1[ i ]
				out_file.write( a1.strip() + "+" + barcode1[7:] + "/" + barcode1[:7] +"\n")
				out_file.write( b1.strip() +"\n")
				out_file.write( c1.strip()+"\n")
				out_file.write( d1.strip() +"\n")
				out_file2 = files2[ i ]
				out_file2.write( a2.strip() + "+" + barcode1[7:] + "/" + barcode1[:7] +"\n")
				out_file2.write( b2[10:].strip() +"\n")
				out_file2.write( c2.strip()+"\n")
				out_file2.write( d2[10:].strip() +"\n")
				determined=determined + 1
				break

        # Continue reading
		a1 = file1.readline()
		Total = Total + 1
		if Total > Goal:
			Goal+=1000000
			print(str(Total)+' reads demultiplexed')


	print ('Total reads =' + str( Total ) )
	print ('Determined reads =' + str( determined ) )
	print ('Undetermined reads =' + str( Total - determined ) )

# Use the function to subset the fastq file
fastq_spliter3(path1,path2,nadaptor)



for f in files1:
    f.close()

for f in files2:
    f.close()
