#!/usr/bin/python


from Bio import SeqIO
from decimal import *


def matchNucleotides(a,b):
	if b =="-" and a=="-":
		return 0
	elif b==a:
		return 2
	else:
		return -1

def maximum_from_tab_separated_matrix_with_labels_file(inputfile):
	maximum_list = []
	with open(inputfile) as f:
		lines = f.readlines()
		line_list_list = []
		line_list_list2 = []
		x = 0
		for line in lines:	
			x = x + 1
			line = line.rstrip("\n")
			line = line.rstrip("\r")
			line_list = line.split('\t', 1)
			line_list2 = line_list[1].split('\t')
			if x > 1:
				line_list_list2.append(line_list2)
			else:
				line_list_list.append(line_list2)
		y = 0
		max_score = 0
		for d in line_list_list2:
			x = 0
			y = y + 1
			for n in d:
				x = x + 1
				if x != y and Decimal(n) > Decimal(max_score):
					max_score = n
		y = 0		
		for d in line_list_list2:
			x = 0
			y = y + 1
			for n in d:
				x = x + 1
				if Decimal(n) == Decimal(max_score):					
					maximum_list.append('|'.join(sorted([line_list_list[0][y-1],line_list_list[0][x-1]])))
	return max_score,list(set(maximum_list))	


print ('Task 4:\n')
input_file = 'input1.fasta'
my_array = []
my_header = ['    ']
for sequence in SeqIO.parse(open(input_file), "fasta"):
	my_header.append(sequence.id)
my_array.append(my_header)
for sequence in SeqIO.parse(open(input_file), "fasta"):
	row = [sequence.id]
	for sequence2 in SeqIO.parse(open(input_file), "fasta"):
		x = 0
		for count in range(0,len(sequence2.seq.tostring())):
			x = x + matchNucleotides(sequence.seq.tostring()[count],sequence2.seq.tostring()[count])
		row.append(x)
	my_array.append(row)
	
for row in my_array:
    for val in row:
        print ('{:7}'.format(val),)
    print ()
print ('\n*******************************\n')

for input_file in ('data2.txt','data.txt','input1.txt','input2.txt'):
	print ('Task 5 met {0}'.format(input_file))
	max_score,max_list = maximum_from_tab_separated_matrix_with_labels_file(input_file)
	print ('Highest overlapscore={0}'.format(max_score))
	for m in max_list:
		maxseq = m.split('|')
		print ('sequence "{0}" and sequence "{1}"'.format(maxseq[0],maxseq[1]))
	print ('have the highest overlap score in the file\n')
	print ('\n*******************************\n'	)


