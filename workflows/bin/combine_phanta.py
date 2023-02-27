#!/opt/conda/bin/python

import sys
import re
import os
import numpy as np
import pandas as pd
from functools import reduce

pd.options.mode.chained_assignment = None

tax_dict ={'species':'S', 'genus':'G','family':'F',
	'order':'O','class':'C','phylum':'P','kingdom':'D'}
tax_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

bracken_out_dir = sys.argv[1] + '/'
kraken_out_dir = sys.argv[2] + '/'

def read_bracken(x):
	sample = x.removesuffix('.out.bracken').removeprefix('phanta_')
	# get the kraken/bracken output
	# check if the file is empty!
	filesize = os.stat(bracken_out_dir + x).st_size
	if filesize==0:
		kraken_file = pd.read_csv(kraken_out_dir + '/phanta_' + sample + '.out.filtered', 
			sep="\t", 
			names=['fraction','fragments', 'assigned', 'rank', 
			'minimizers','uniqminimizers',
			'ncbi_taxID','sci_name'], nrows=2)
		all_reads = kraken_file['fragments'].sum()
		temp = pd.DataFrame({'taxonomy': ['unclassified'], sample: [all_reads]})
	else :
		bracken_file = pd.read_csv(bracken_out_dir + x, sep="\t", 
			names=['fraction','fragments', 'assigned', 'rank', 
			'ncbi_taxID','sci_name'])
		all_bracken_reads = bracken_file.loc[0,'fragments']
		bracken_file['sci_name'] = bracken_file['sci_name'].str.strip()
		# take only the levels we are interested in
		bracken_red = bracken_file.loc[bracken_file['rank'].isin(tax_dict.values())]

		# make nice tax names
		# first add kingdom level
		bracken_red['kingdom'] = np.nan
		bracken_red.loc[bracken_red['rank']=='D','kingdom'] =\
	 	'k_' + bracken_red.loc[bracken_red['rank']=='D','sci_name']
		bracken_red['kingdom'].fillna(method='ffill', inplace=True)

		for tl in range(1, len(tax_levels)):
			lvl = tax_levels[tl]
			bracken_red[lvl] = np.nan
			bracken_red.loc[bracken_red['rank']==tax_dict[lvl],lvl] = \
			lvl[0].lower() + '_' + \
			bracken_red.loc[bracken_red['rank']==tax_dict[lvl],'sci_name']
			bracken_red[lvl] = bracken_red.groupby(tax_levels[:tl])[lvl].fillna(method='ffill')
			bracken_red[lvl].fillna('', inplace=True)

		# nicer names
		bracken_red['nice_name'] = bracken_red[tax_levels].apply(lambda x: ';'.join(x), axis=1)
		bracken_red['nice_name'] = bracken_red['nice_name'].str.replace("\s+", "_", regex=True)
		bracken_red['nice_name'] = bracken_red['nice_name'].str.replace(";+$", "", regex=True)

		# return only what is needed
		temp = bracken_red[['nice_name', 'fragments']]
		
		temp.rename(columns={'nice_name': 'taxonomy', 'fragments': sample}, inplace=True)

		# add non-resolved reads from kraken as unclassified
		kraken_file = pd.read_csv(kraken_out_dir + '/phanta_' + sample + '.out.filtered', 
			sep="\t", 
			names=['fraction','fragments', 'assigned', 'rank', 
			'minimizers','uniqminimizers',
			'ncbi_taxID','sci_name'], nrows=2)
		all_reads = kraken_file['fragments'].sum()
		temp_uncl = pd.DataFrame({'taxonomy': ['unclassified'], 
			sample: [all_reads - all_bracken_reads ]})
		temp = pd.concat([temp, temp_uncl])

	return(temp)


files_bracken = os.listdir(bracken_out_dir)

dfs = [read_bracken(f) for f in files_bracken]
df_merged = reduce(lambda  left,right: pd.merge(left,right, on=['taxonomy'],
	how='outer'), dfs).fillna(0)
df_merged.sort_values('taxonomy', inplace=True)
df_merged.to_csv('phanta_all.tsv', sep="\t", index=False)
