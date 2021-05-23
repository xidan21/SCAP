#!/usr/bin/env python



import os
import sys
import re

# argv[1] = "express_matrix.txt"
# argv[2] = "genes_title.txt"
# argv[3] = "genes_length.txt"



print ">>>>>>>>>>>>>>>>>>>>>>>"
print sys.argv[1]

print sys.argv[2]

print sys.argv[3]


os.system("R --vanilla --no-save --args %s %s %s < scrna_1.r > ../result/temp.log" %(sys.argv[1], sys.argv[2], sys.argv[3]))

os.system("open ../result/detected_transcripts.pdf ")