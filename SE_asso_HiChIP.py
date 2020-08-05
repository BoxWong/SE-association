#!/bin/python

#To find associataed genes using HiChIP interaction maps
#Xin Wang, xin.wang@uni-goettingen.de
#20200326

import argparse
import os
import sys
import pybedtools

parser = argparse.ArgumentParser(description='To find associated genes for certain peaks')
parser.add_argument('-bp', type=argparse.FileType('r'), help='bedpe file to determine the associated genes')
parser.add_argument('-a', type=str, help='Gene body annotation file in bed format')
parser.add_argument('-q', type=str, help='Bed file of interests for region-gene association')
parser.add_argument('-e', type=int, default=0, help='Extension length upstream and downstream of the gene body')
parser.add_argument('-o', type=argparse.FileType('w'), help='File name to write output')

args=parser.parse_args()
#
query = pybedtools.BedTool(args.q)
bp_left='' ; bp_right=''
#Make two separate bed files out of bedpe files.
ln=0
bedpe = args.bp.readlines()
for i in bedpe:
    ln+=1
    line = i.strip().split('\t')
    bp_left = bp_left + line[0] + '\t' + line[1] + '\t' +line[2] +'\t' + str(ln) + '\n'
    bp_right = bp_right + line[3] + '\t' + line[4] + '\t' +line[5] +'\t' + str(ln) + '\n'
bp_left_bed = pybedtools.BedTool(bp_left, from_string = True)
bp_right_bed = pybedtools.BedTool(bp_right, from_string = True)
del bp_left; del bp_right

#Intersect the query bed with left and right bed files.
left_bp_query = query.intersect(bp_left_bed, wao = True)
right_bp_query = query.intersect(bp_right_bed, wao = True)

#get the interactions which are overlapped with query regions
query_bedpe_dict = dict()
for i in left_bp_query:
    a = str(i).strip().split()
    query_item = a[0] + '\t' + a[1] + '\t' + a[2]
    if query_item not in query_bedpe_dict.keys():
        query_bedpe_dict[query_item] = {a[6]}
    else:
        query_bedpe_dict[query_item].add(a[6])
for i in right_bp_query:
    a = str(i).strip().split()
    query_item = a[0] + '\t' + a[1] + '\t' + a[2]
    if query_item not in query_bedpe_dict.keys():
        query_bedpe_dict[query_item] = {a[6]}
    else:
        query_bedpe_dict[query_item].add(a[6])

#extend the gene body by the distance provided
Complete_gene_set = set()
gene_body = ''
for i in pybedtools.BedTool(args.a):
    a = str(i).strip().split('\t')
    gene_body = gene_body + a[0] + '\t' + str(max([(int(a[1])-args.e),1])) + '\t' + str(int(a[2]) + args.e) + '\t' + a[3] + '\n'
gene_body_bed = pybedtools.BedTool(gene_body, from_string = True)
To_write=''

for i in query:
    a = str(i).strip().split('\t')
    query_item = a[0] + '\t' + a[1] + '\t' + a[2]
    query_item_bed = pybedtools.BedTool(query_item, from_string = True)
    query_re_gene_bed = query_item_bed.intersect(gene_body_bed, wb = True) #get the genes overlapped with query regions
    set1 = set()
    for j in query_re_gene_bed:
        b = str(j).strip().split('\t')
        set1.add(b[6])
    if query_bedpe_dict[query_item] == {'.'}:#if the query region doesn't overlap with any interactions
        #print(query_item,'doesn\'t overlap with any interaction from your bedpe file...')
        pass
    else:
        #print(query_item, 'has overlap with interactions from your bedpe file, looking for distally regulated genes...')
        set2 = query_bedpe_dict[query_item] - {'.'}
        #print(query_bedpe_dict[query_item], set2)
        interaction = ''
        for k in set2:
            c = bedpe[int(k)-1].strip().split('\t')
            interaction = interaction + c[0] + '\t' + c[1] + '\t' + c[2] + '\n' + c[3] + '\t' + c[4] + '\t' + c[5] + '\n'
        #print(interaction)
        interaction_bed = pybedtools.BedTool(interaction, from_string = True)
        interaction_re_gene_bed = interaction_bed.intersect(gene_body_bed, wb = True) # get genes associataed with interactions
        for k in interaction_re_gene_bed:
            c = str(k).strip().split('\t')
            set1.add(c[6])
    Complete_gene_set = Complete_gene_set.union(set1)
    for l in set1:
        To_write = To_write + query_item + '\t' + l + '\n'

print('Found',len(Complete_gene_set),' genes associated with the query.\n writing to file...')       
args.o.write(To_write)
args.o.close()
