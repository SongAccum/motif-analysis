import pandas as pd
from collections import Counter
import os
import sys
import re

#in_file='/home/songyabing/jupyter/motif_analysis/dot.info'
#1547544533723_m1_run1.search.txt---argv[1]
#out_folder:argv[2]
#out folder:sys.argv[3]
if '/' in sys.argv[1]:
    sample_name=re.search('(.*)/(.*)_BEAMready_(m.).*',sys.argv[1]).group(2)+'_'+\
                 re.search('(.*)/(.*)_BEAMready_(m.).*',sys.argv[1]).group(3)
else:
    sample_name=re.search('(.*)_BEAMready_(m.).*',sys.argv[1]).group(1)+'_'+\
                 re.search('(.*)_BEAMready_(m.).*',sys.argv[1]).group(2)
    
dot_file=sys.argv[2]+'/'+sample_name+'.dot'
cmd="sed -n '/^#DB/,/#PSSM$/p' "+ sys.argv[1] \
     +" |awk '{if ($1!~/^#/ &&  length($1)!=0 && $2!~/bg/ ) print;}' | cut -f 1 >"+dot_file
os.system(cmd)
#remove bg motif info 
#os.system("wc -l "+sys.argv[2])
with open(dot_file,'r') as dot:
    dlst=[]
    dot_line=dot.readline().strip()
    while dot_line:
        dot_lst=[i for i in dot_line]
        dlst.append(dot_lst)
        dot_line=dot.readline().strip()

RNA_chr='G'*len(dlst[0])
sym='>'+sample_name+'\n'+RNA_chr+'\n'

out_file=sys.argv[2]+'/'+sample_name+'.dbn'
with open(out_file,'w') as out:
    for i in range(len(dlst[0])):
        result=Counter([j[i] for j in dlst])
        sym+=max(result,key=result.get)
        if i == (len(dlst[0])-1):
            #print (sym)
            out.write(sym)
software="/BioII/lulab_b/songyabing/motif_analysis/software/VARNA"
out_png=sys.argv[2]+'/'+sample_name+'.png'
visual_cmd="java -cp "+software+'/'+\
            "VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -baseName '#FFFFFF' -baseInner '#FFFFFF' -resolution 2. -i "+\
            out_file +" -o "+out_png
#print (visual_cmd)
os.system(visual_cmd)



