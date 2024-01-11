#!/bin/python

import pysam
import sys
import os
import re

# print(sys.argv[1])
mysam = sys.argv[1]
# mysam = "5hmc_1P_clean.tag.bam"
# tmpsam = ''.join([mysam.split(".bam")[0], "_rmDup.tmp.bam"]) 
outsam = ''.join([mysam.split(".bam")[0], "_rmDup.bam"])  ## for expression matrix
sam_rd = pysam.AlignmentFile(mysam, 'rb')

read_dic = {}
read_unique = []
n_total = 0
n_highQuality = 0
n_unique = 0
# with pysam.AlignmentFile(mysam, 'rb') as sam_rd:
for read in sam_rd:
    n_total = n_total + 1
    # if(read.mapping_quality >= 30):
    #     n_highQuality = n_highQuality+1
    # ID = '-'.join([read.reference_name, str(read.pos), read.cigarstring, read.get_tag('CB')])
    map_mat = [int(x) for x in re.split('[A-Z]', read.cigarstring)[0:-1]] ## length of frag
    state = re.findall('[A-Z]', read.cigarstring) ## mapping state of frag
    [map_mat[i] for i in [x.start() for x in re.finditer('M|N|D', ''.join(state))]]
    # read_end = int(read.pos) + sum([map_mat[i] for i in [x.start() for x in re.finditer('M|N|D', ''.join(state))]])
    read_end = read.reference_end
    # print(" ".join([str(read.pos),str(read_end),read.flag,read.cigarstring]))
    if(read.flag == 0):
        ID = '-'.join([read.reference_name, str(read.pos), read.get_tag('CB')])
        if(ID in read_dic):
            if (read_dic[ID].mapping_quality >= read.mapping_quality):
                continue
            else:
                read_dic[ID] = read
        else:
            n_unique = n_unique +1
            read_dic['-'.join([read.reference_name, str(read.pos), read.get_tag('CB')])] = read
    elif(read.flag == 16):
        ID = '-'.join([read.reference_name, str(read_end), read.get_tag('CB')])
        if(ID in read_dic):
            if (read_dic[ID].mapping_quality >= read.mapping_quality):
                continue
            else:
                read_dic[ID] = read
        else:
            n_unique = n_unique +1
            read_dic['-'.join([read.reference_name, str(read_end), read.get_tag('CB')])] = read

for key in read_dic:
    read_unique.append(read_dic[key])

sam_wrt = pysam.AlignmentFile(outsam + ".tmp", 'wb', template=sam_rd)
for r in read_unique:
    sam_wrt.write(r)
sam_wrt.close()

pysam.sort(outsam + ".tmp", '-o', outsam)
os.remove(outsam + ".tmp")
 
sam_rd.close()
sam_wrt.close()
pysam.index(outsam)
print("Total mapped read pairs:" + str(n_total))
# print("Confidently mapped read pairs:" + str(n_highQuality))
print("Unique DNA fragments:" + str(n_unique))