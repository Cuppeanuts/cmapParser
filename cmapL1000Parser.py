#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.write_gctx as wg

gene_path = 'GSE92742/GSE92742_Broad_LINCS_gene_info.txt'
sig_path = 'GSE92742/GSE92742_Broad_LINCS_sig_info.txt'
data_path = 'GSE92742/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'

genesymbol = raw_input("Enter Gene Symbol: ")
geneid = raw_input("Enter Gene ID: ")

#GET gene data from datamatrix
if not os.path.exists("raw_%s.gctx" % genesymbol):
	print "Extract Gene Raw Data... "
	gene_info = pd.read_csv(gene_path, sep='\t', dtype=str)
	landmark_gene_row_ids = gene_info["pr_gene_id"][gene_info["pr_gene_symbol"] == genesymbol]
	landmark_only_gctoo = parse(data_path, rid = landmark_gene_row_ids)
	wg.write(landmark_only_gctoo, "raw_%s.gctx" % genesymbol)
	print "Finished!"

print "Analysing Data... "
#GET sig info
sig_info = pd.read_csv(sig_path, sep='\t', dtype=str)
itemlist = [u'sig_id', u'pert_id', u'pert_iname', u'pert_type', u'cell_id',u'pert_idose', u'pert_itime']
annolist = np.vstack([sig_info[item].values for item in itemlist])
#print sig_info.columns
#print annolist
#print np.unique(sig_info[u'pert_type'].values,return_counts = True)
index = dict([[annolist[0,i], annolist[:,i]] for i in range(len(annolist[0,:]))])

#ADD annotations
data = parse("raw_%s.gctx" % genesymbol)
val = data.data_df.values.T
colname = data.data_df.columns.values
annoinfo = np.array([index[sid] for sid in colname])
context = np.hstack([val, annoinfo])

#REMOVE control
t_group = ['trt_cp', 'trt_lig']
condition = np.array([item in t_group for item in context[:,4]])
context = context[condition]
np.shape(context)
np.set_printoptions(threshold = np.nan)
#print np.unique(context[np.array([item[:4] == 'BRD-' for item in context[:,3]]),3], return_counts = True)

#FILTER pert
filterpert = context[np.where(context[:,0].astype('float16') > 1.644854)]
sortpert = filterpert[np.argsort(filterpert[:,0])][::-1].astype('str')
header = np.array(['z-score','sig_id','pert_id', 'pert_iname', 'pert_type', 'cell_id','pert_idose', 'pert_itime'])
f1 = np.vstack([header, sortpert])
f2 = np.vstack([f1[:,3], f1[:,0],f1[:,4:].T]).T

#WRITE file
with open('%s_cmap_L1000_Activation.txt' % genesymbol, 'w') as fo:
	fo.write('\n'.join(['\t'.join(row) for row in f2.tolist()]))
print "Finished!"

if "y" in raw_input("Find Suppressive molecules? (y/n)").lower():
	#FILTER pert
	filterpert = context[np.where(context[:,0].astype('float16') < -1.644854)]
	sortpert = filterpert[np.argsort(filterpert[:,0])].astype('str')
	header = np.array(['z-score','sig_id','pert_id', 'pert_iname', 'pert_type', 'cell_id','pert_idose', 'pert_itime'])
	f1 = np.vstack([header, sortpert])
	f2 = np.vstack([f1[:,3], f1[:,0],f1[:,4:].T]).T
	
	#WRITE file
	with open('%s_cmap_L1000_Suppression.txt' % genesymbol, 'w') as fo:
		fo.write('\n'.join(['\t'.join(row) for row in f2.tolist()]))
	print "Finished!"
