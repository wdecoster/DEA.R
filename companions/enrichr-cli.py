# Script to perform EnrichR enrichment analysis
# Based on code available on enrichR website,
# adapted by wdecoster

import json
import requests
import sys
import os
import argparse

def readDBfile():
	'''
	Reads possible databases from a file in the installation dir of this script.
	Dangerous if stuff is moved around!
	Symlinks are okay.
	'''
	db = os.path.join(os.path.dirname(os.path.realpath(__file__)), "enrichr-databases.txt")
	if not os.path.isfile(db):
		sys.exit("ERROR: Could not find my database location. Did you move files around?")
	else:
		with open(db) as databases:
			return([line.strip() for line in databases.readlines() if not line == ""])


def getArgs():
	parser = argparse.ArgumentParser(description="Perform enrichment analysis for a gene list and multiple databases.")
	parser.add_argument("-g", "--genes",
					help="A genelist to be queried, either a list in a file or '-' for a list of genes on stdin")
	parser.add_argument("-d", "--databases",
					help="Databases to query, omit to use a default set.",
					nargs='*')
	parser.add_argument("-w", "--which",
					help="List databases which can be queried and quit.",
					action='store_true')
	parser.add_argument("-p", "--prefix",
					help="Fixed prefix to name the output files with",
					default="")
	parser.add_argument("-o", "--outdir",
					help="Output directory to store files in. Will be created if it doesn't exist.",
					default='.')
	args = parser.parse_args()
	if not args.genes and not args.which:
		sys.exit("Input Error: Required argument is -g/--genes containing list of genes to use for enrichment.")
	return args

def senddata(genes):
	'''
	Send the input gene list to enrichr, return query
	Call function to check how many genes were recognized
	'''
	input = {
		'list': (None, '\n'.join(genes)),
		'description': (None, 'enrichR.py query')
		}
	response = requests.post('http://amp.pharm.mssm.edu/Enrichr/addList', files=input)
	if not response.ok:
		raise Exception('Error uploading gene list')
	queryId = json.loads(response.text)['userListId']
	askgenelist(queryId, genes)
	return(queryId)


def askgenelist(id, inlist):
	'''
	Compare the genes send and received to get succesfully recognized genes
	'''
	response = requests.get('http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s' % id)
	if not response.ok:
		raise Exception('Error getting gene list back')
	returnedL = json.loads(response.text)["genes"]
	returnedN = sum([1 for gene in inlist if gene in returnedL])
	print('{} genes succesfully recognized by Enrichr'.format(returnedN))


def whichdb():
	'''
	Check user-chosen databases or use a default set
	'''
	standarddb = ['KEGG_2015', 'BioCarta_2016', 'WikiPathways_2016', 'Reactome_2016',
					'GO_Biological_Process_2015', 'GO_Cellular_Component_2015', 'GO_Molecular_Function_2015',
					'MSigDB_Computational', 'Panther_2016']
	if args.databases:
		chosendb = [db for db in args.databases if db in databases]
		if chosendb:
			return chosendb
		else:
			print("Could not find a valid database in arguments, falling back on default databases.")
			return standarddb
	else:
		print("Using {} default databases.".format(len(standarddb)))
		return standarddb


def procesinput():
	'''
	check input to the script, either as a file or on stdin
	'''
	if args.which:
		for option in databases:
			print(option)
		sys.exit(0)
	if args.outdir:
		if not os.path.exists(args.outdir):
			os.makedirs(args.outdir)
	if args.genes == '-':
		print("Expecting input on stdin.")
		genes = set([item.strip() for item in sys.stdin.readlines() if not item == ""])
	else:
		if os.path.isfile(args.genes):
			with open(args.genes) as inputgenes:
				genes = set([item.strip() for item in inputgenes.readlines() if not item == ""])
		else:
			sys.exit("Input file not found, is the path correct?")
	print('Input contains {} unique gene names'.format(len(genes)))
	return genes


def getresults(id, gene_set_library):
	'''
	Receive the enrichment for the chosen databases
	write to files with default names
	'''
	filename = gene_set_library + '_enrichment'
	url = 'http://amp.pharm.mssm.edu/Enrichr/export?userListId=%s&filename=%s&backgroundType=%s' % (id, filename, gene_set_library)
	response = requests.get(url, stream=True)
	with open(os.path.join(args.outdir, args.prefix + filename + '.txt'), 'wb') as output:
		for chunk in response.iter_content(chunk_size=1024):
			if chunk:
				output.write(chunk)

if __name__ == '__main__':
	args = getArgs()
	databases = readDBfile()
	genelist = procesinput()
	id = senddata(genelist)
	for db in whichdb():
		getresults(id, db)
