#!/usr/bin/python

from Bio import SeqIO


input_file = '/home/gehau/svn/testing/src/main/python/biopython/test.fasta'

amb_nuc_hash = {'R' : ['A','G'],'Y' : ['C','T'],'S' : ['G','C'],'W' : ['A','T'],'K' : ['G','T'],'M' : ['A','C'],'B' : ['C','G','T'],'D' : ['A','G','T'],'H' : ['A','C','T'],'V' : ['A','C','G'],'N' : ['A','C','G','T']} 
ambiguousNucleotides = {'A' : '[ARWMDHVN]','T' : '[TYWKBDHN]','C' : '[CYSMBHVN]','G' : '[GRSKBDVN]','R' : '[RAG]','Y' : '[YCT]','S' : '[SGC]','W' : '[WAT]','K' : '[KGT]','M' : '[MAC]','B' : '[BCGT]','D' : '[DAGT]','H' : '[HACT]','V' : '[VACG]','N' : '[NACGT]'}
unambiguousNucleotides = {'A' : '[ARWMDHVN]','T' : '[TYWKBDHN]','C' : '[CYSMBHVN]','G' : '[GRSKBDVN]'}
ambiguousNucleotidesNo = {'N': 4,"V": 3,"H": 3,"D": 3,"B": 3,"M": 2,"K": 2,"W": 2,"S": 2,"Y": 2,"R": 2}
my_not_unique_list = []


def numberOfPossibleSequences(sequenceString):	
	no = 1
	for key in ambiguousNucleotidesNo:
		no = no * ambiguousNucleotidesNo[key] ** list(sequenceString).count(key)			
	return  no
	
def regexpBuilder(sequenceString,subString):
	regexp = []	
	for char in list(sequenceString):
		regexp.append(subString[char])					
	return  ''.join(regexp)
	
def regexpMatcher(regexp,sequenceString):
	import re
	pattern = re.compile(regexp)
	match = pattern.findall(sequenceString)
	return match
	
def sequenceMultiplier(sequenceString):	
	newSequenceList = []
	for key in amb_nuc_hash:
		if not list(sequenceString).count(key) == 0 :					
			for nuc in amb_nuc_hash[key]:						
				sequenceList = sequenceString.split(key, 1)	
				#print amb_nuc				 
				newSequenceList.append(('%s') % nuc.join(sequenceList))				
	return newSequenceList
	
def sequenceEvaluator(sequenceString):
	newSequenceList = []
	sequenceList = sequenceMultiplier(sequenceString)
	for sequenceString in sequenceList:
		hasAmbiguousNucletides = False
		for key in amb_nuc_hash:
			if not list(sequenceString).count(key) == 0 :	
				hasAmbiguousNucletides = True				
		if not hasAmbiguousNucletides:
			my_not_unique_list.append(sequenceString)
		if hasAmbiguousNucletides:
			newSequenceList.append(sequenceString)
	if len(newSequenceList) > 0 :
		for sequenceString in newSequenceList:
			sequenceEvaluator(sequenceString)
			

def inputSequence():		
	seq = None
	import re
	pattern = re.compile('[ATGCRYSWKMBDHVN]') 
	pattern2 = re.compile('[RYSWKMBDHVN]') 
	while seq is None:
		isValidSequence = True
		isAmbiguousSequence = 'N'
		input_value = raw_input("Insert sequence: ")
		try:	
			for s in list(input_value):				
				match = pattern.match(s)
				match2 = pattern2.match(s)
				if not match :
					isValidSequence = None	
				if match2 :
					isAmbiguousSequence = 'Y'
			if isValidSequence:
				seq = input_value		
		except ValueError:
			print "{input} is not correct".format(input=input_value)


	return seq,isAmbiguousSequence
	
	
	
#A FASTA file contains a number of DNA sequences. Unfortunately, some of the symbols are ambiguous. 
#The encoding is IUPAC (http://www.bioinformatics.org/sms/iupac.html). Write a Python script that, given the name of the FASTA file,
# writes the sequence identifier and the number of possible sequences for each sequence in the file.  		
print 'Task 1:\n'	
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:  
	number = numberOfPossibleSequences(fasta.seq.tostring())
	print 'filename: {0}\nsequenceId: {1}\nthe number of possible sequences: {2}\n'.format(input_file,fasta.id,number)
print '\n*******************************\n'	     
  
  
  
        
#As above, a FASTA file contains sequences with ambiguous symbols (IUPAC). Write a function that, given the name of the FASTA file
# and an unambiguous DNA string, writes out the identifiers of the sequences of which the given sequence could be a subsequence. 
#Try to implement this without generating all possible sequences. 
#Implement the case also when the subsequence can have ambiguous symbols as well as the sequences in the FASTA file.
print 'Task 2:\n'	
seq,yesOrNo=inputSequence()
ambiguous = {'Y' : ambiguousNucleotides,'N' : unambiguousNucleotides}
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
	regExp =  regexpBuilder(seq,ambiguous[yesOrNo])
	match = regexpMatcher(regExp,fasta.seq.tostring())
	print regExp 
	print match
	if len(match) > 0 :
		print 'filename: {0}\n"{2}" could be a subsequence of "{1}"\n'.format(input_file,fasta.id,seq)
print '\n*******************************\n'	




#Write a script that, given the name of a FASTA file and a sequence identifier, writes all possible sequences, 
#but only if there are less than 100. If the given sequence identifier does not occur in the FASTA file, 
#or if there are more than 100 possible sequences, an appropriate error message is printed.
print 'Task 3:\n'
print 'What is the sequence identifier?'
prompt = '> '
f_id = raw_input(prompt)	
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
sequences_found = None
for fasta in fasta_sequences:		
	if fasta.id == f_id :
		print 'the number of possible sequences: {0}\n'.format(numberOfPossibleSequences(fasta.seq.tostring()))
		sequences_found = True
		if numberOfPossibleSequences(fasta.seq.tostring()) < 100 :
			sequenceEvaluator(fasta.seq.tostring())
		else :
			print  'The number of possibilities exeeds 99'
if not sequences_found:
	print  'sequence not found'
	
for unique_sequence in list(set(my_not_unique_list)):
	print unique_sequence
print '\n*******************************\n'	



#Write a script that given a FASTA file that contains aligned sequences (no ambiguous symbols), 
#writes a matrix of the overlap scores of the sequences. A match has a score of 2, a mismatch a score of -1, 
#while a gap scores -1 as well.  If both sequences share a gap, the score is 0. Label the rows and columns of the matrix by the sequence identifiers.

print 'Task 4:\n'
	
