# Script to perform EnrichR enrichment analysis
# Based on code available on enrichR website,
# adapted by wdecoster

import json
import requests
import sys
import os
import argparse


databases = ['Achilles_fitness_decrease', 'Achilles_fitness_increase', 'Aging_Perturbations_from_GEO_down', 'Aging_Perturbations_from_GEO_up',
	'Allen_Brain_Atlas_down', 'Allen_Brain_Atlas_up', 'BioCarta_2013', 'BioCarta_2015', 'BioCarta_2016', 'Cancer_Cell_Line_Encyclopedia', 'ChEA_2013', 'ChEA_2015',
	'Chromosome_Location', 'CORUM', 'dbGaP', 'Disease_Perturbations_from_GEO_down', 'Disease_Perturbations_from_GEO_up', 'Disease_Signatures_from_GEO_down_2014', 'Disease_Signatures_from_GEO_up_2014',
	'Drug_Perturbations_from_GEO_2014', 'Drug_Perturbations_from_GEO_down', 'Drug_Perturbations_from_GEO_up', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
	'ENCODE_Histone_Modifications_2013', 'ENCODE_Histone_Modifications_2015', 'ENCODE_TF_ChIP-seq_2014', 'ENCODE_TF_ChIP-seq_2015', 'Epigenomics_Roadmap_HM_ChIP-seq',
	'ESCAPE', 'Genes_Associated_with_NIH_Grants', 'GeneSigDB', 'Genome_Browser_PWMs', 'GO_Biological_Process_2013', 'GO_Biological_Process_2015', 'GO_Cellular_Component_2013',
	'GO_Cellular_Component_2015', 'GO_Molecular_Function_2013', 'GO_Molecular_Function_2015', 'GTEx_Tissue_Sample_Gene_Expression_Profiles_down', 'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
	'HMDB_Metabolites', 'HomoloGene', 'Human_Gene_Atlas', 'Human_Phenotype_Ontology', 'HumanCyc_2015', 'Humancyc_2016', 'KEA_2013', 'KEA_2015', 'KEGG_2013', 'KEGG_2015', 'KEGG_2016',
	'Kinase_Perturbations_from_GEO_down', 'Kinase_Perturbations_from_GEO_up', 'Ligand_Perturbations_from_GEO_down', 'Ligand_Perturbations_from_GEO_up', 'LINCS_L1000_Chem_Pert_down',
	'LINCS_L1000_Chem_Pert_up', 'LINCS_L1000_Kinase_Perturbations_down', 'LINCS_L1000_Kinase_Perturbations_up', 'LINCS_L1000_Ligand_Perturbations_down', 'LINCS_L1000_Ligand_Perturbations_up',
	'MCF7_Perturbations_from_GEO_down', 'MCF7_Perturbations_from_GEO_up', 'MGI_Mammalian_Phenotype_2013', 'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotype_Level_4',
	'Microbe_Perturbations_from_GEO_down', 'Microbe_Perturbations_from_GEO_up', 'Mouse_Gene_Atlas', 'MSigDB_Computational', 'MSigDB_Oncogenic_Signatures', 'NCI-60_Cancer_Cell_Lines',
	'NCI-Nature_2015', 'NCI-Nature_2016', 'NURSA_Human_Endogenous_Complexome', 'Old_CMAP_down', 'Old_CMAP_up', 'OMIM_Disease', 'OMIM_Expanded', 'Panther_2015', 'Panther_2016', 'Pfam_InterPro_Domains',
	'Phosphatase_Substrates_from_DEPOD', 'PPI_Hub_Proteins', 'Reactome_2013', 'Reactome_2015', 'Reactome_2016', 'SILAC_Phosphoproteomics', 'Single_Gene_Perturbations_from_GEO_down',
	'Single_Gene_Perturbations_from_GEO_up', 'TargetScan_microRNA', 'TF-LOF_Expression_from_GEO', 'Tissue_Protein_Expression_from_Human_Proteome_Map', 'Tissue_Protein_Expression_from_ProteomicsDB',
	'Transcription_Factor_PPIs', 'TRANSFAC_and_JASPAR_PWMs', 'Virus_Perturbations_from_GEO_down', 'Virus_Perturbations_from_GEO_up', 'VirusMINT', 'WikiPathways_2013', 'WikiPathways_2015', 'WikiPathways_2016']


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
	genelist = procesinput()
	id = senddata(genelist)
	for db in whichdb():
		getresults(id, db)
