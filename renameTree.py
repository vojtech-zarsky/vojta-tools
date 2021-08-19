import argparse
import re
from ete3 import Tree


parser = argparse.ArgumentParser(prog='Rename leaves of a newick tree. (ete3 dependency)')


parser.add_argument('--inTree', required=True)
parser.add_argument('--table', required=True, help='table format: original leaf name<tab>new leaf name')
parser.add_argument('--outTree', required=True)
args = parser.parse_args()

dLeaf2Name = {}
dNames = {}
for sLine in open(args.table):
    lLine = re.split('\t',sLine.strip(),1)
    if len(lLine) > 1:
        (sId, sName) = lLine
        sId = sId.strip()
        sName = sName.strip()
    else:
        print('error in line:', sLine.strip())
        continue
    sName = re.sub('\W+',' ',sName)
    sName = sName.strip().replace(' ','_')
    if len(sId) == 0 or len(sName) == 0:
        print('error in line:', sLine.strip())
        continue
    if sId not in dLeaf2Name:
        dLeaf2Name[sId] = sName
    else:
        print('WARNING:', sId, 'leaf ID duplicated!')
    if sName in dNames:
        print('WARNING:', sName, 'leaf name duplicated!')
    dNames[sName] = None

tree = Tree(args.inTree)

for leaf in tree.get_leaves():
    if leaf.name in dLeaf2Name:
        leaf.name = dLeaf2Name[leaf.name]
    else:
        print('WARNING:', leaf.name, 'not in the table!')

tree.write(outfile=args.outTree,format=0)
