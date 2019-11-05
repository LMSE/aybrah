

"""

BETA
Run from scripts folder

"""

import pandas as pd
import json
import datetime
from Bio import SeqIO

def write(data_json):
    with open('../aybrah_'+version_new+'.json', 'w') as outfile:  
        json.dump(data_json, outfile,sort_keys=True,indent=4)
    open('../version.txt','w').write(version_new)





comment_update='Remove duplicates by reviewing phylo trees'

version=open('../version.txt','r').read()
digit=2
version_new='.'.join([v if i !=digit else str(int(v)+1) for i,v in enumerate(version.split('.'))])


aybrah=pd.read_excel('../aybrah.xlsx',sheet_name='aybrah').fillna('').set_index('FOG')
aybrah_oid=pd.read_excel('../aybrah.xlsx',sheet_name='organisms').fillna('').set_index('oid')


oids=aybrah_oid.index.tolist()
fogs=aybrah.index.tolist()



data={}

data['homologs']={}




while len(fogs) > 0:
    fog=fogs.pop(0)
    print(fog)
    hog=aybrah.HOG[fog]
    current={}
    # get parent group
    parent=aybrah.Parent[fog]
    if len(parent)>0:
        parent_type=parent.split(':')[0]
        parent_fog=parent.split(':')[1].split('|')[0]
    # get relevant information
    current['info']={key:aybrah[key][fog]for key in\
         ['Protein ID','sce Locus','Features','Protein description',\
                'Parent','Neofunc, Subfunct','Suggested Analysis'] if len(aybrah[key][fog])>0}
    # process seqids into paralogs or orthologs
    seqids=filter(None,[aybrah[oid][fog] for oid in oids])
    seqids_orthologs=filter(None,[seqid for seqid in seqids if ';' not in seqid])
    seqids_paralogs=filter(None,';'.join([seqid for seqid in seqids if ';' in seqid]).split(';'))
    seqids_orthologs={seqid:{} for seqid in seqids_orthologs}
    seqids_paralogs={seqid:{} for seqid in seqids_paralogs}
    #
    if len(seqids_orthologs)>0:
        current['orthologs']=seqids_orthologs
    if len(seqids_paralogs)>0:
        current['paralogs']=seqids_paralogs
    if hog not in data['homologs'].keys():
        data['homologs'][hog]={}
        data['homologs'][hog]['ortholog_group']={}
    data['homologs'][hog]['ortholog_group'][fog]=current



data['organisms']={oid:{col:aybrah_oid[col][oid] for col in aybrah_oid.columns} for oid in oids }

data['meta']={}

data['meta']['version']=version_new
data['meta']['modified']=str(datetime.datetime.now())
data['meta']['comment']=comment_update

write(data)


##############################################################################################################################




from collections import OrderedDict

def get_seqids_by_oid(oid):
	return ';'.join(sorted([seqid for seqid in seqids if seqid.split('|')[0]==oid]))


aybrah = json.load(open('../aybrah_'+version_new+'.json'), object_pairs_hook=OrderedDict)



oids=zip(*sorted([(aybrah['organisms'][oid]['order'],oid) for oid in aybrah['organisms'].keys()]))[1]


hogs=aybrah['homologs'].keys()

headers=['FOG','HOG','Parent']+list(oids)



open('../aybrah_'+version_new+'.tsv','w').write('\t'.join(headers)+'\n')

for hog in hogs:
	print hog
	fogs=aybrah['homologs'][hog]['ortholog_group'].keys()
	for fog in fogs:
		df_fog=aybrah['homologs'][hog]['ortholog_group'][fog]
		try:
			orthologs=df_fog['orthologs'].keys()
		except:
			orthologs=[]
		try:
			paralogs=df_fog['paralogs'].keys()
		except:
			paralogs=[]
		seqids=orthologs+paralogs
		try:
			parent=df_fog['info']['Parent']
		except:
			parent=''
		new_row=[   ';'.join(sorted([seqid for seqid in seqids if seqid.split('|')[0]==oid])) if len(get_seqids_by_oid(oid))>0 else ''   for oid in oids]
		new_row=[fog,hog,parent]+new_row
		open('../aybrah_'+version_new+'.tsv','a').write('\t'.join(new_row)+'\n')



###########################################################################################


import xml.etree.ElementTree as ET
from xml.dom import minidom



# process_seqids(ortholog_fog,seqids)

def process_seqids(ortholog,seqids):
	seqids_orthologs=filter(None,[seqid for seqid in seqids if ';' not in seqid])
	seqids_paralogs=filter(None,';'.join([seqid for seqid in seqids if ';' in seqid]).split(';'))
	for seqid in seqids_orthologs:
		xml_geneid=str(lookup[1][seqid.split('|')[1]])
		#xml_geneid=seqid.split('|')[1]
		# assign to current fog
		generef=ET.SubElement(ortholog,'geneRef')
		generef.set('id',xml_geneid)
	if len(seqids_paralogs)>0:
		ortholog_fogp=ET.SubElement(ortholog,'paralogGroup')
		ortholog_fogp.set('id',fog+'_p')
		prop=ET.SubElement(ortholog_fogp,'property')
		prop.set('name','paralog_type')
		prop.set('value','indistinguishable_duplication')
		for seqid in seqids_paralogs:
			# assign to new paralog group (these must be processed)
			#xml_geneid=seqid.split('|')[1]
			xml_geneid=str(lookup[1][seqid.split('|')[1]])
			generef=ET.SubElement(ortholog_fogp,'geneRef')
			generef.set('id',xml_geneid)

def get_descendants(index):
	indicies_temp=[index]
	fogs=[aybrah_tsv['FOG'][index] for index in indicies_temp]
	children=df_temp[df_temp.Parent.str.contains('|'.join(fogs))].index.tolist()
	children=[child for child in children if child not in indicies_temp]
	while len(set(children)-set(indicies_temp))>0:
		indicies_temp=indicies_temp+children
		fogs=[aybrah_tsv['FOG'][index] for index in indicies_temp]
		children=df_temp[df_temp.Parent.str.contains('|'.join(fogs))].index.tolist()
		children=[child for child in children if child not in indicies_temp]
	return indicies_temp




aybrah = json.load(open('../aybrah_'+version_new+'.json'), object_pairs_hook=OrderedDict)


oids=zip(*sorted([(aybrah['organisms'][oid]['order'],oid) for oid in aybrah['organisms'].keys()]))[1]





#taxon=pd.read_csv('updated_taxons.txt',sep='\t').fillna('')

"""
find . -name "YeastProteome*"
"""



print('create header')

root = ET.Element('orthoXML')
root.set('xmlns',"http://orthoXML.org/2011/")
root.set('version','0.3')
root.set('origin','AYbRAH')
root.set('originVersion',version_new)
root.set('xmlns:xsi',"http://www.w3.org/2001/XMLSchema-instance")
root.set('xsi:schemaLocation',"http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd")


# contains all the sequence information from MycoCosm and Uniprot
#YP=pd.read_csv('/Volumes/5TB/Organized/github/aybrah_old/YeastProteome.txt',sep='\t',header=None)
#YP=pd.read_csv('/Volumes/SD250GB/Organized/github/aybrah_old/YeastProteome.txt',sep='\t',header=None)

records=SeqIO.parse('../yeast_proteome_aybrah.faa','fasta')



yeast_proteome={oid:[]for oid in oids}

for record in records:
	seqid=record.id
	oid=seqid.split('|')[0]
	yeast_proteome[oid].append(seqid)




aybrah_tsv=pd.read_csv('../aybrah_'+version_new+'.tsv',sep='\t').fillna('')


# counter for sequences
i=0

# create lookup file to convert between seqid and id for orthoxml
f=open('../lookup_xml_genes.txt','w')
f.write('')
f.close()

print('add species')

for oid,row in aybrah_oid.iterrows():
	#oid=taxon['oid'][index]
	species = ET.SubElement(root, 'species')
	species.set('name',aybrah_oid['name'][oid])
	species.set('NCBITaxId',str(int(aybrah_oid['NCBITaxId'][oid])))
	database=ET.SubElement(species,'database')
	url_genome=aybrah_oid['geneLink'][oid]
	database.set('geneLink',url_genome)
	if 'jgi' in url_genome:
		database.set('name','MycoCosm')
		database.set('version',str(aybrah_oid['version'][oid]))
	else:
		database.set('name','Uniprot')
	genes=ET.SubElement(species,'genes')
	#seqids=YP[YP[0].str.contains(oid+'\|')][0].tolist()
	seqids=yeast_proteome[oid]
	seqids_aybrah=[s for s in ';'.join(aybrah_tsv[oid].tolist()).split(';') if 'AYbRAH' in s]
	seqids=seqids+seqids_aybrah
	for seqid in seqids:
		gene=ET.SubElement(genes,'gene')
		i=i+1
		gene.set('id',str(i))
		gene.set('geneId',seqid.split('|')[1])
		f=open('../lookup_xml_genes.txt','a')
		f.write('\t'.join([seqid.split('|')[1],str(i)])+'\n')
		f.close()



lookup=pd.read_csv('../lookup_xml_genes.txt',sep='\t',header=None).set_index(0)




hogs=filter(None,list(sorted(set(aybrah_tsv.HOG.tolist()))))

groups=ET.SubElement(root,'groups')

print('process orthologs')





for hog in hogs:
	print('\t'+hog)
	ortholog_hog=ET.SubElement(groups,'orthologGroup')
	ortholog_hog.set('id',hog)
	df_temp=aybrah_tsv[aybrah_tsv.HOG==hog]
	indicies_base=[index for index,row in df_temp.iterrows() if len(aybrah_tsv.Parent[index])==0]
	process=[]
	for index in indicies_base:
		process.extend(get_descendants(index))
	for index in process:
	#for index,row in df_temp.iterrows():
		# create ortholog group
		fog=df_temp.FOG[index]
		#yog=df_temp.YOG[index]
		print '\t\t'+fog
		parent=aybrah_tsv.Parent[index]
		seqids=filter(None,[aybrah_tsv[oid][index] for oid in oids if 'missing' not in aybrah_tsv[oid][index]])
		# make nested within another ortholog group
		if 'ortholog' in parent:
			ortholog_fog=[ortho for ortho in root.findall('.//orthologGroup')+root.findall('.//paralogGroup') if parent.split(':')[1].split('|')[0]==ortho.attrib['id']][0]
			ortholog_fog.set('synonym',fog)
			process_seqids(ortholog_fog,seqids)
		# should be ortholog group, not paralog
		elif len(parent)>0:
			test=aybrah_tsv[aybrah_tsv.FOG==parent.split(':')[1].split('|')[0]].Parent.tolist()[0]
			# what is this about???
			if 'ortholog' in test:
				parent_old=parent
				parent=test
			#ortholog_fog=[ortho for ortho in root.findall('.//orthologGroup')+root.findall('.//paralogGroup') if parent.split(':')[1].split('|')[0]==ortho.attrib['id']][0]
			ortholog_fog=[ortho for ortho in root.findall('.//orthologGroup')+root.findall('.//paralogGroup') if parent.split(':')[1].split('|')[0]==ortho.attrib['id']][0]
			sub_ortholog_fog=ET.SubElement(ortholog_fog,'orthologGroup')
			sub_ortholog_fog.set('id',fog)
			prop=ET.SubElement(sub_ortholog_fog,'property')
			prop.set('name','paralog_type')
			if 'ortholog' in test:
				prop.set('value',parent_old.split(':')[0])	
			else:
				prop.set('value',parent.split(':')[0])
			print 'process derived paralogs'
			process_seqids(sub_ortholog_fog,seqids)
		else:
			ortholog_fog=ET.SubElement(ortholog_hog,'orthologGroup')
			ortholog_fog.set('id',fog)
			process_seqids(ortholog_fog,seqids)






xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(indent="\t")
with open('../aybrah_'+version_new+'.xml', "w") as f:
    f.write(xmlstr)






