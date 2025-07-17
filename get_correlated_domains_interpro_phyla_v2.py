import sys
import operator
from scipy.stats import pearsonr

thresholds = []
subgroups = []
f = open('Phyla_all')
for line in f:
	line = line.strip()
	subgroups.append(line.split()[0])
	thresholds.append(float(line.split()[1]))
f.close()

for i in range(0, len(subgroups)):
	subgroup = subgroups[i]
	threshold = thresholds[i]
	f = open('interpro_files.txt')
	for line in f:
		query = line.strip()
		print(subgroup,query)
		domains = {}
		descriptions = {}
		f = open('interpo_all_reps_domains_clusterid_'+subgroup)
		for line in f:
			if line.split('\t')[0] != 'Domain':
				name = line.split('\t')[0]
				descriptions[name] = line.split('\t')[1]
				values = [float(i) for i in line.split('\t')[2:]]
				if sum(values) == 0 or sum(values) == len(values):
					continue
				domains[name] = values
		f.close()
		if query not in domains:
			print(query+' not found')
			continue
		domain_set = set()
		unsorted_correlations = {}
		query_values = domains[query]
		for subject in domains:
			subject_values = domains[subject]
			coefficient, pvalue = pearsonr(query_values, subject_values)
			if coefficient > threshold:
				unsorted_correlations[subject+'\t'+descriptions[subject]+'\t'+str(coefficient)+'\t'+str(pvalue)] = coefficient
				domain_set.add(subject)
		print('Interpro_new_correlations/correlations_'+query.replace(':','-')+'_'+subgroup+'.txt')
		output = open('Interpro_new_correlations/correlations_'+query.replace(':','-')+'_'+subgroup+'.txt', 'w')
		output.write('domain\tdescription\tcorrelation_coefficient\tpvalue\n')
		count = 0
		for a,b in sorted(unsorted_correlations.items(), key = operator.itemgetter(1), reverse = True):
			output.write(a+'\n')
		output.close()
	f.close()
