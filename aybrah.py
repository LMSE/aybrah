
import pandas as pd
from ete3 import Tree



class AYBRAH:
	def __init__(self,path_aybrah,load_sequences=False):
		print('\tLoad AYbRAH')
		self.path_aybrah=path_aybrah
		self.path_file_aybrah=self.path_aybrah+'/aybrah.xlsx'
		self.path_file_phylosignal=self.path_aybrah+'/phylogeny/phylosignal-sequences.faa'
		self.path_file_species_tree=self.path_aybrah+'/phylogeny/species-tree.newick'
		self.path_file_proteome=self.path_aybrah+'/yeast_proteome_aybrah.faa'
		#
		self.df=pd.read_excel(self.path_file_aybrah,sheet_name='aybrah').fillna('').set_index('FOG')
		#self.df=pd.read_csv('latest_aybrah.txt',sep='\t').fillna('').set_index('FOG')
		self.organisms=pd.read_excel(self.path_file_aybrah,sheet_name='organisms').fillna('').set_index('oid')
		try:
			self.lookup_fog_by_entry=pd.read_csv(self.path_aybrah+'/entry-to-fog.txt',sep='\t').fillna('').set_index('entry')
		except:
			print('Need to run "aybrah.update_table()"')
			self.lookup_fog_by_entry=None
		self.organisms=pd.read_excel(self.path_file_aybrah,sheet_name='organisms').set_index('oid').fillna('')
		self.oids=['lma','ddi','mbr','gla',
			'bdb', 'uma', 'cne', 'pgr', 'rgm', 'spo', 'sai', 'ssl', 'tre', 'ncr', 'pno', 'ani', 'afm', 'ang',
			'lst', 'yli', 'arx', 'nfu',
			'aru', 'wan', 'cut', 'hva', 'kla', 'erc', 'ago', 'lkl', 'kwal', 'lth', 'tdl', 'zro', 'tbl', 'tpf', 'vpo', 'ndi', 'ncs', 'kaf', 'knag', 'cgr', 'suva', 'skud', 'smik', 'sce', 'pta', 'kpa', 'ppa', 'kcp', 'car', 'opol', 'opm', 'dbx', 'pku', 'pme', 'bin', 'caur', 'clu', 'mbi', 'cten', 'mgu', 'dha', 'pic', 'spa', 'lel', 'cpa', 'cmt', 'cot', 'ctp', 'cdu', 'calwo', 'cal']
		#
		try:
			self.species_tree=Tree(self.path_file_species_tree)
		except:
			self.species_tree=open(self.path_aybrah+'/phylogeny/species-tree.txt','r').read()
		# manual overrides
		self.organisms.loc['pku']['UniProt Proteome']='UP000249293'
		self.df['arx']['FOG00903']='arx|A0A060TFF3'
		# load sequences creates a dictionary of all sequence records in AYbRAH
		if load_sequences==True:
			self.YP={record.id:record for record in SeqIO.parse(self.path_file_proteome,'fasta')}
		else:
			self.YP=None
	def __repr__(self):
		return "AYbRAH DB"
	def open_organisms(self,application="Sublime Text"):
		if application=='Excel':
			application='Microsoft Excel'
		os.system('open -a "{}" {}'.format(application,self.path_file_organisms))
	def save_organisms(self):
		self.organisms.order=[int(cell) if len(str(cell))>0 else '' for cell in self.organisms.order]
		self.organisms.NCBITaxId=[int(cell) if len(str(cell))>0 else '' for cell in self.organisms.NCBITaxId]
		self.organisms.to_csv(self.path_file_organisms,sep='\t')
	def open_species_tree(self):
		os.system('open -a "Sublime Text" {}'.format(self.path_file_species_tree))
	def open(self):
		os.system('open -a "Microsoft Excel" {}'.format(self.path_aybrah)  )
	def get_hogs_from_fogs(self,fogs):
		hogs=sorted(set(self.df[self.df.index.str.contains('|'.join(fogs))].HOG))
		return hogs
	def get_hog_from_fog(self,fog):
		return self.df.HOG[fog]
	def get_fogs_from_hogs(self,hogs):
		fogs=sorted(set(self.df[self.df.HOG.str.contains('|'.join(hogs))].index.tolist()))
		return fogs
	def get_entries_from_hogs(self,hogs):
		fogs=self.get_fogs_from_hogs(hogs)
		return [entry for fog in fogs for oid in oids for entry in self.get_entries_by_fog_and_oid(fog,oid) ]
	def get_seqids_from_hogs(self,hogs):
		fogs=self.get_fogs_from_hogs(hogs)
		return [seqid for fog in fogs for oid in self.oids for seqid in self.get_seqids(fog,oid)]
	def get_seqids_from_fog(self,fog):
		return [seqid for oid in self.oids for seqid in self.get_seqids(fog,oid)]
	def get_seqids(self,fog,oid,manual_annotation=True):
		seqids=list(filter(None,self.df[oid][fog].split(';')))
		if len(seqids)==0:
			return []
		if manual_annotation==False:
			seqids=[seqid for seqid in seqids if '_AYbRAH_' not in seqid]
		return seqids
	def get_oid_from_entry(self,entry):
		return self.lookup_fog_by_entry['seqid'][entry].split('|')[0]
	def get_entries_by_fog_and_oid(self,fog,oid,manual_annotation=False):
		seqids=self.get_seqids(fog,oid,manual_annotation=False)
		if seqids==None:
			return None
		return [seqid.split('|')[1] for seqid in seqids]
	def create_lookup_table(self):
		"""
		used to create entry-to-fog table, faster to find FOG
		given an protein entry
		"""
		table=pd.DataFrame(columns=['entry','FOG','seqid'])
		for fog,row in self.df.iterrows():
			print(fog)
			seqids=list(filter(None,[self.get_seqids(fog,oid) for oid in self.organisms.index]))
			seqids=[s for seqid in seqids for s in seqid]
			for s in seqids:
				table.loc[len(table)]=[s.split('|')[1],fog,s]
		table.to_csv(path_aybrah+'/entry-to-fog.txt',sep='\t',index=False)
	def get_fogs_from_gene_name(self,gene_name):
		# requires sce to be loaded...
		orf=sce.table[(sce.table.feature=='gene') & (sce.table.symbol==gene_name)].index.item()
		entries=sce.mapping.get_entries_from_orf(orf)
		return [self.get_fog_from_entry(entry) for entry in entries]
	def get_seqids_for_oid(self,oid):
		# returns all seqids for given oid
		return list(filter(None,';'.join(self.df[oid]).split(';')))
	def get_fogs_from_entries(self,entries):
		return [self.get_fog_from_entry (entry) for entry in entries]
	def get_fog_from_entry(self,entry):
		if entry not in self.lookup_fog_by_entry.index:
			return None
		return self.lookup_fog_by_entry.loc[entry].FOG
	def get_fog_from_seqid(self,seqid):
		oid,entry=seqid.split('|')
		return self.get_fog_from_entry(entry)
	def get_root(self,hog):
		for oid in self.oids:
			for fog in self.get_fogs_from_hogs([hog]):
				for seqid in self.get_seqids(fog,oid):
					#raise Exception
					return seqid
	def get_entries_by_fogs_and_oids(self,fogs,oids):
		entries=[entry for fog in fogs for oid in oids for entry in self.get_entries_by_fog_and_oid(fog,oid,True)]
		return entries
		#			
	"""
	Functions below requires additional modules not yet integrated into AYbRAH
	"""
	def get_sce_ortholog_info(self,entry):
		# requires get_gene_info_from_entry function
		if entry==None:
			return None,None
		fog=self.get_fog_from_entry(entry)
		if fog==None:
			return None,None
		entries_oid_sce=self.get_oid_ortholog(fog,'sce')
		entries_sce=[x.split('|')[1] for x in entries_oid_sce]
		orfs_sce=', '.join([get_gene_info_from_entry(entry_sce,'Gene_Name') if len(get_gene_info_from_entry(entry_sce,'Gene_Name'))!=0 else get_gene_info_from_entry(entry_sce,'Gene_OrderedLocusName') for entry_sce in entries_sce])
		return fog,orfs_sce
	def get_entries_by_fog_and_oid(self,fog,oid,return_seqid=False):
		e1=self.df[oid][fog].split(';')
		e2=self.df[self.df.Parent.str.contains('ortholog:{}'.format(fog))][oid].tolist()
		fog_equivalent=self.df.Parent[fog].replace('ortholog:','')
		if fog_equivalent in self.df.index:
			e3=self.df[oid][fog_equivalent].split(';')
		else:
			e3=[]
		entries_oid=list(filter(None,e1+e2+e3))
		if return_seqid==True:
			return entries_oid
		else:
			entries=[e.split('|')[1] for e in entries_oid]
			return entries
	def get_orf_by_locus_identifer(self,locus_identifer,genome_reference):
		# should be a unique lookup
		# first tries gene name like ALD6
		orf_reference=genome_reference.mapping.get_orf_from_gene_name(locus_identifer)
		entry_reference=genome_reference.mapping.get_entry_from_gene_name(locus_identifer)
		# then tries wildcards like YLR108C
		if entry_reference==None:
			entry_reference=genome_reference.mapping.get_entry_from_wildcard(locus_identifer)
		return orf_reference
		#return Gene(entry_reference,genome_reference)
	"""
	def get_genes_via_reference_genome_orthology(self,locus_identifer,genome_subject,genome_reference):
		# mapping index is Entry
		gene_reference_info=self.get_gene_by_locus_identifer(locus_identifer,genome_reference)
		#
		entries_subject=self.get_entries_by_fog_and_oid(gene_reference_info.fog,genome_subject.oid)
		genes=[Gene(entry,genome_subject) for entry in entries_subject if 'AYbRAH' not in entry]
		for gene in genes:
			gene.sce_ortholog=gene_reference_info
		return genes
	"""
	#def get_orfs_from_fogs(self,fogs):
	def get_orfs_from(self,genome_subject,entries):
		orfs=genome_subject.mapping.get_orfs_from_entries(entries)
		return orfs
	def get_orfs_via_reference_genome_orthology(self,locus_identifier,genome_subject,genome_reference):
		# mapping index is Entry
		#entry_reference=genome_reference.mapping.get_entry_from_wildcard(locus_identifer)
		entries_reference=genome_reference.mapping.get_entry_from_gene_name(locus_identifier)
		fogs=[self.get_fog_from_entry(entry) for entry in entries_reference]
		fogs+=[self.df.Parent.loc[fog].split(':')[1].split('?')[0] for fog in fogs if 'ortholog:' in self.df.Parent.loc[fog]]
		fogs+=[f for fog in fogs for f in self.df[self.df.Parent.str.contains('ortholog:{}'.format(fog))].index.tolist()]
		#aybrah.df.loc[fogs[0]]['arx']
		entries_subject=[entry for fog in set(fogs) for entry in self.get_entries_by_fog_and_oid(fog,genome_subject.oid)]
		if len(entries_subject)==0:
			return []
		orfs_to_search=genome_subject.mapping.get_orfs_from_entries(entries_subject)
		orfs_to_search+=[genome_subject.mapping.get_field_name_from_entry(entry,'EnsemblGenome') for entry in entries_subject]
		orfs_to_search_text='|'.join(['ID=gene-{};'.format(orf) for orf in orfs_to_search])
		matched=genome_subject.gff[genome_subject.gff[8].str.contains(orfs_to_search_text)]
		orfs=[field.split('ID=gene-')[1].split(';')[0] for field in matched[8].tolist()]
		return orfs

#path_aybrah='/Users/kcorreia/terminal/Genotecha-example/aybrah-aybrah-expanded'
path_aybrah='/Users/kcorreia/terminal/Genotecha-example/aybrah'

aybrah=AYBRAH(path_aybrah)

