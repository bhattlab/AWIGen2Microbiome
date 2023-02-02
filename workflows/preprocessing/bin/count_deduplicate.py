#!/usr/bin/env python3
import json
import sys

sample_id = sys.argv[1]
stats_log = sys.argv[2]

with open(stats_log) as user_file:
	file_contents = user_file.read()

parsed_json = json.loads(file_contents)
parsed_json = parsed_json[0]

reads_in = parsed_json.get('Paired_end').get('in') + parsed_json.get('Single_end').get('in')
reads_out = parsed_json.get('Paired_end').get('out') + parsed_json.get('Single_end').get('out')

print('\t'.join((sample_id, 'Deduplication input', str(reads_in))))
print('\t'.join((sample_id, 'Deduplication output', str(reads_out))))
print('\t'.join((sample_id, 'Duplicated reads', str(reads_in-reads_out))))

