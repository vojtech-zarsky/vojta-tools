## DEPENDENCIES:
## Taxonomy.py (https://github.com/vojtech-zarsky/vojta-tools/blob/master/Taxonomy.py)
## ete3 (http://etetoolkit.org/)
## USAGE:
## python speciesToTree.py  --inFile <list of taxids of species names> --taxdump <NCBI taxdump directory> [--library <location of Taxonomy.py file>]

import sys
import re
import argparse
from ete3 import Tree, TreeStyle, TextFace, faces


parser = argparse.ArgumentParser(prog='Build a tree from a list of species. DEPENDENCIES: Taxonomy.py (https://github.com/vojtech-zarsky/vojta-tools/blob/master/Taxonomy.py), ete3 (http://etetoolkit.org/)')

parser.add_argument('--inFile', required=True, help='list of taxids or species names')
parser.add_argument('--taxdump', required=True, help='NCBI taxdump directory, can be downloaded here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz')
parser.add_argument('--library', required=False, default='./', help='Location of Vojta\'s Taxonomy.py file. Default: current directory')
args = parser.parse_args()

sInFile = args.inFile
sTaxdumpDir = args.taxdump
sLib = args.library
sLib = sLib.split('Taxonomy.py')[0]

sys.path.append(sLib)
from Taxonomy import Taxonomy

#sInFile = 'speciesToTree_test.txt'
#sTaxdumpDir = './taxdump/'

def taxonomies2Tree(lTaxonomies):
    ## bild tree
    tree = Tree(name='root')
    for lTaxonomy in lTaxonomies:
        if len(lTaxonomy) == 0:
            continue
        node = tree
        iTaxLevel = 0
        for nodeTemp in tree.traverse():
            if nodeTemp.name in lTaxonomy:
                node = nodeTemp
                iTaxLevel = lTaxonomy.index(nodeTemp.name)+1
        for iTaxLevel in range(iTaxLevel, len(lTaxonomy)):
            node.add_child(Tree(name=lTaxonomy[iTaxLevel]))
            for child in node.get_children():
                if child.name == lTaxonomy[iTaxLevel]:
                    node = child
                    break
    ## remove nodes with one child
    for node in tree.traverse():
        if len(node.get_children()) == 1:
            if len(node.get_ancestors()) > 0:
                ancestor = node.get_ancestors()[0]
                ancestor.add_child(node.get_children()[0])
                ancestor.remove_child(node)
            else:
                tree = node.get_children()[0]
    return tree

taxonomy = Taxonomy(sTaxdumpDir=sTaxdumpDir)
lTaxonomies = []
iLine = 0
for sLine in open(sInFile):
    iLine += 1
    sTaxId = sLine.strip()
    if sTaxId == '':
        continue
    lTaxonomy = []
    if sTaxId.isdigit():
        lTaxonomy = taxonomy.getTaxonomyByTaxId(sTaxId)
    else:
        lTaxonomy = taxonomy.getTaxonomy(sTaxId)
    if len(lTaxonomy) < 4:
        print('SHORT/MISSING TAXONOMY line:{}, taxid:{}, taxonomy:{}'.format(iLine, sTaxId, lTaxonomy))
        raise
    lTaxonomy = list(map(lambda x:'{} {}'.format(re.sub('^ a-zA-Z0-9_','',x[1]), x[0]), lTaxonomy))[::-1]
    lTaxonomies.append(lTaxonomy)

## for testing
"""
lTaxonomies = [['cellular organisms', 'eukaryota', 'opisthokonta', 'metazoa', 'eumetazoa', 'bilateria', 'deuterostomia', 'chordata', 'craniata', 'vertebrata', 'homo sapiens'],\
    ['cellular organisms', 'eukaryota', 'opisthokonta', 'metazoa', 'eumetazoa', 'bilateria', 'deuterostomia', 'chordata', 'tunicata', 'ascidiacea', 'phlebobranchia', 'cionidae', 'ciona', 'ciona intestinalis'],\
    ['cellular organisms', 'bacteria', 'terrabacteria group', 'tenericutes', 'mollicutes', 'mycoplasmatales', 'mycoplasmataceae', 'mycoplasma', 'mycoplasma pneumoniae'],\
    ['cellular organisms', 'archaea', 'euryarchaeota', 'methanomada group', 'methanococci', 'methanococcales', 'methanococcaceae', 'methanococcus', 'methanococcus maripaludis']]
#"""

tree = taxonomies2Tree(lTaxonomies)

## write a table with the taxonomies
lTaxonomies = []
for leaf in tree.get_leaves():
    lTaxonomy = [leaf.name]+list(map(lambda x:x.name, leaf.get_ancestors()))
    lTaxonomy = lTaxonomy[::-1]
    if lTaxonomy[0] == 'root':
        del lTaxonomy[0]
    lTaxonomies.append( lTaxonomy )
tableOut = open('{}.tsv'.format(sInFile),'w')
for lTaxonomy in sorted(lTaxonomies):
    tableOut.write('\t'.join(lTaxonomy)+'\n')
tableOut.close()

## print the tree
ts = TreeStyle()
ts.show_leaf_name = False
def my_layout(node):
    F = TextFace(node.name, tight_text=True)
    faces.add_face_to_node(F, node, column=0, position="branch-right")
ts.layout_fn = my_layout
#tree.show(tree_style=ts) ## to show the tree
tree.render('{}.pdf'.format(sInFile), tree_style=ts)

## write the tree
for node in tree.traverse():
    node.name = '"{}"'.format(node.name)
tree.write(format=8, outfile='{}.phb'.format(sInFile))
