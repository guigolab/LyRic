#!/usr/bin/env python
import pandas as pd
from pprint import pprint
import json


#########################
#### Process ENTRIES ####
#########################

entriesDf = pd.read_table('entries.tsv', header=0, sep='\t')
entriesDf['entry_id2'] = entriesDf['entry_id']
entriesDf.set_index('entry_id2', inplace=True)

# convert certain columns to JSON:
for column in ['samples', 'contacts', 'experiment_ids', 'library_preps', 'platforms']:
    entriesDf[column] = entriesDf[column].astype("string")
    entriesDf[column] = entriesDf[column].apply(json.loads)

#entriesDf['team_id'] = entriesDf['team_id'].astype("string")
# print one JSON file per row in input TSV
# this line does not output the index columns, this is why we need a copy of entry_id above
entriesDf.apply(lambda x: x.to_json(
    "entries/{}/entry.json".format(x.entry_id)), axis=1)
#    "entries/{}/entry.json".format(x.name)), axis=1)

#############################
#### Process EXPERIMENTS ####
#############################

experimentsDf = pd.read_table('experiments.tsv', header=0, sep='\t')
experimentsDf['experiment_id2'] = experimentsDf['experiment_id']
experimentsDf.set_index('experiment_id2', inplace=True)
experimentsDf = experimentsDf.fillna('')
experimentsDf['description'] = 'LyRic transcript models for ' + \
    experimentsDf['species'] + ' ' + experimentsDf['experiment_id']
experimentsDf = experimentsDf.drop(['file', 'entry_id'], axis=1)
# convert certain columns to JSON:
for column in ['platforms', 'library_preps', 'samples', 'libraries', 'software']:
    experimentsDf[column] = experimentsDf[column].astype("string")
    experimentsDf[column] = experimentsDf[column].apply(json.loads)

experimentsDf['notes'] = experimentsDf['notes'].astype("string")
# the following cannot use 'entry_id' to name output path and skip it in the json file, so we'll have to move the output later:
experimentsDf.apply(lambda x: x.to_json(
    "entries/{}.experiment.json".format(x.experiment_id)), axis=1)
