#!/opt/conda/bin/python

import sys
import numpy as np
import pandas as pd

kraken_report = sys.argv[1]
kraken_db = sys.argv[2]
filter_threshold = int(sys.argv[3])

# get the kraken file and merge with the inspect file to get the coverage
kraken_file = pd.read_csv(kraken_report, sep="\t", 
	names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 
	'classification_rank','ncbi_taxa','scientific name'])
inspect_file = pd.read_csv(kraken_db + '/inspect.out', sep="\t", 
	names=['frac','minimizers_clade', 'minimizers_taxa', 'rank',
	'ncbi_taxonomy','sci_name'])
kraken_full = kraken_file.merge(inspect_file, 
	left_on='ncbi_taxa', right_on='ncbi_taxonomy')
kraken_full['cov'] = kraken_full['uniqminimizers']/kraken_full['minimizers_taxa']
kraken_full.loc[kraken_full.minimizers_taxa < 5,'cov']=np.nan
kraken_full['sci_name'] = kraken_full['sci_name'].str.strip()

# make new columns for the domain and species information
kraken_full['domain'] = np.nan
kraken_full.loc[kraken_full['rank']=='D','domain']=kraken_full.loc[kraken_full['rank']=='D','sci_name']
kraken_full['domain'].fillna(method='ffill', inplace=True)
kraken_full['species'] = np.nan
kraken_full.loc[kraken_full['rank']=='S','species']=kraken_full.loc[kraken_full['rank']=='S','sci_name']
kraken_full['species'].fillna(method='ffill', inplace=True)

# get only genomes and filter out sub-species (e.g. S1, when S2 is present)
species_kraken = kraken_full.copy()[kraken_full['rank'].str.startswith('S', na=False)]
# sub-species stuff
# get difference to next value, keep only those with zero (or higher) :D
species_kraken['subspecies'] = pd.to_numeric(species_kraken['rank'].str.replace('S', ''))
species_kraken['subspecies'] = species_kraken['subspecies'].fillna(0)
species_kraken['change'] = species_kraken.groupby(['species'])['subspecies'].diff(periods=-1)
species_kraken['change'] = species_kraken['change'].fillna(0)
species_kraken = species_kraken.loc[species_kraken['change'] >= 0]

# get max coverage and max uniqminimizers
dflookup = species_kraken[['domain', 'species', 'cov', 'uniqminimizers']].groupby(['domain', 'species'], as_index=False).max()
dflookup['pass'] = False
# if the coverage is above the threshold, keep the species
dflookup['cov_pass'] = np.where(dflookup['domain'] == 'Viruses', 0.10, 0.01)
dflookup['uniqminimizers_pass'] = np.where(dflookup['domain'] == 'Viruses', 0, 0)
# these are the default parameters for phanta... For now, we can 
# ignore the number of uniq minimizers, I guess?
dflookup.loc[dflookup['cov'] >= dflookup['cov_pass'], 'pass'] = True
species_keep = np.append(dflookup.loc[dflookup['pass']==True,'species'].values,
	np.array(['not_species_level']))

# new kraken report
kraken_file = kraken_file.merge(kraken_full[['species', 'ncbi_taxa']], on='ncbi_taxa')

# make sure that everything above S/S1 is included :D
kraken_file['species_keep']= 'not_species_level'
kraken_file.loc[kraken_file['classification_rank'].str.startswith('S', na=False),'species_keep'] = kraken_file.loc[kraken_file['classification_rank'].str.startswith('S', na=False),'species']

# filter the kraken file and write it out
kraken_file = kraken_file[kraken_file['species_keep'].isin(species_keep)]
del kraken_file['species']
del kraken_file['species_keep']
kraken_file.to_csv(kraken_report + '.filtered', sep="\t", index=False, header=False)

# check if braken would fail
# if yes, make output files and skip braken
temp = kraken_file.loc[kraken_file['classification_rank'].str.startswith('S', na=False),]
temp = temp.loc[temp['fragments'] >= filter_threshold,]
if temp.shape[0] == 0:
	with open(kraken_report+'.bracken', 'w') as outfile:
		pass
