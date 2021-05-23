#!/usr/bin/env python



import os
import sys
import re

# argv[1] = "express_matrix.txt"
# argv[2] = "genes_title.txt"
# argv[3] = "genes_length.txt"
print sys.argv[5]



x = open("../source/seq.txt",'r')

out = open("../source/marker_genes.txt",'w')

for line in x:

	if not re.search("^$", line):
	
	
		array_genes = line.split(",")
		
		for i in xrange(len(array_genes)):
		
			print >> out, re.sub(" ","", array_genes[i])
			
out.close()			

os.system("R --vanilla --no-save --args %s %s %s %s %s %s < scrna_2.r > ../result/temp.log" %(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[7]))

os.system("open ../result/tsne_plot_*.pdf")

os.system("open ../result/*_marker_genes.pdf")

#os.system("rm ../source/seq.txt")
