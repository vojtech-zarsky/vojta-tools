import sys

sHelpText = 'LGTanalysis (https://github.com/vojtech-zarsky/vojta-tools)\nWritten for Python 2\n'\
    '\nPrerequisities:\n ete2 or ete3, Taxonomy.py from my GitHub (https://github.com/vojtech-zarsky/vojta-tools) in case you want to get taxonomy from NCBI taxids\n'\
    '\nSyntax:\npython LGTanalysis.py --tree_file --support_cutoff --orthology_cutoff --query_sequence_id --taxonomy_of_selected_groups --taxonomy_mapping --tax_id_mapping --taxdumpdir --outfile\n'\
    '\n--tree_file e.g. IQtree output'\
    '\n--support_cutoff Nodes with support bellow this will be removed. (mind e.g. 90 vs 0.9)'\
    '\n--orthology_cutoff Cutoff of the orthology score. (default 0.5)'\
    '\n--query_sequence_id ID of the query sequence.'\
    '\n--taxonomy_of_selected_groups A tab-separated table of selected eukaryotic and prokaryotic groups with their taxonomy. (e.g. https://github.com/vojtech-zarsky/vojta-tools/blob/master/LGTanalysis.groups.tsv) An arrow indicates inclusion of the taxon on the left site to the group on the right. The domain-level taxon should be first.'\
    '\n--taxonomy_mapping A table of sequence IDs with their taxonomies. The domain-level taxon should be first. (e.g. A2FH21_TRIVA<tab>Eukaryota<tab>Metamonada<tab>Parabasalia<tab>Trichomonadida)'\
    '\n--tax_id_mapping Alternatively a table of sequence IDs with their NCBI taxonomies can be provided. (e.g. A2FH21_TRIVA<tab>5722) In that case the Taxonomy.py file is needed and the --taxdump_dir directory must be specified.'\
    '--outfile File where to write the analysis output.'\
    '\n--taxdump_dir directory to the NCBI taxdump folder which can be downloaded here: https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'\
    '\n'

dArgs = {}
lArgs = []
for iArgIndex in range((len(sys.argv)-1)/2):
    sArg = sys.argv[1+iArgIndex*2].split('--')[1]
    sValue = sys.argv[1+iArgIndex*2+1]
    dArgs[sArg] = sValue
    lArgs.append(sArg)

try:
    sTreeFile = dArgs['tree_file']
except:
    print sHelpText
    raise

try:
    if '.' in dArgs['support_cutoff']:
        fSupportCutoff = float(dArgs['support_cutoff'])
    else:
        fSupportCutoff = int(dArgs['support_cutoff'])
except:
    print sHelpText
    raise

fOrthologyCutoff = 0.5
if 'orthology_cutoff' in dArgs:
    fOrthologyCutoff = float(dArgs['orthology_cutoff'])

try:
    sQueryId = dArgs['query_sequence_id']
except:
    print sHelpText
    raise


try:
    sTaxonomyOfSelectedGroups = dArgs['taxonomy_of_selected_groups']
except:
    print sHelpText
    raise

try:
    sTaxonomyMapping = None
    sTaxIdMapping = None
    sTaxdumpDir = None
    if 'taxonomy_mapping' in dArgs:
        sTaxonomyMapping = dArgs['taxonomy_mapping']
    else:
        sTaxIdMapping = dArgs['tax_id_mapping']
        sTaxdumpDir = dArgs['taxdump_dir']
except:
    print sHelpText
    raise

try:
    sOutfile = dArgs['outfile']
except:
    print sHelpText
    raise

fileOut = open(sOutfile,'w')
fileOut.write('Arguments:\n')
for sArg in lArgs:
    fileOut.write('\t'+sArg+': '+dArgs[sArg]+'\n')
fileOut.write('=====\n')
fileOut.flush()
# we need to know about the taxonomy of the sequences
# either by supplying a table with sequence ID - taxonomy mapping
# or e.g. from the taxids of the Uniprot records
"""
First we read a system of taxonomic groups and their relationships. How we define these groups can influence the results a lot.
The table gives us the phylogenetic structure of those groups. The rightmost taxons, that I later on call "groups", should match an NCBI taxonomy group or our supplied taxonomy mapping.
Sometimes an arrow indicates merger of certain groups. E.g. I wanted to merge apicomplexa+chrompodellids into colpodellida, which however isn't annotated in the NCBI taxonomy.
"""

## read taxonomic groups and the phylogeny of eukaryotic groups

try:
    from ete3 import Tree
except ImportError:
    from ete2 import Tree

def readTaxonomy(file):
    lGroups = []
    dTaxon2Group = {}
    dGroup2Taxonomy = {}
    tree = Tree(name='root')
    for sLine in open(file):
        lTaxonomy = sLine.strip().split('\t')
        if '->' in lTaxonomy[-1]:
            (sTaxon, sGroup) = map(lambda x:x.strip(), lTaxonomy[-1].split('->'))
            dTaxon2Group[sTaxon] = sGroup
            del lTaxonomy[-1]
        else:
            dTaxon2Group[lTaxonomy[-1]] = lTaxonomy[-1]
        sGroup = lTaxonomy[-1]
        if sGroup not in dGroup2Taxonomy:
            lGroups.append(sGroup)
            dGroup2Taxonomy[sGroup] = lTaxonomy
            if lTaxonomy[0] == 'eukaryota':
                node = tree
                for iTaxLevel in range(1, len(lTaxonomy)):
                    bChildCheck = False
                    for child in node.children:
                        if child.name == lTaxonomy[iTaxLevel]:
                            node = child
                            bChildCheck = True
                            break
                    if not bChildCheck:
                        childNew = Tree(name=lTaxonomy[iTaxLevel])
                        node.add_child(childNew)
                        node = childNew
    return (lGroups, dGroup2Taxonomy, dTaxon2Group, tree)

(lGroups, dGroup2Taxonomy, dTaxon2Group, treeEukGroups) = readTaxonomy(sTaxonomyOfSelectedGroups)

fileOut.write('Tree of eukaryotic groups:\n')
fileOut.write( treeEukGroups.get_ascii()+'\n' )
fileOut.write('=====\n')
fileOut.flush()

## get taxonomy for each sequence from Uniprot fasta and assign it to one of the selected groups
dSeqId2Group = {}
dSeqId2Taxonomy = {}

if sTaxonomyMapping != None:
    for sLine in open(sTaxonomyMapping):
        lLine = sLine.strip().split('\t')
        sSeqId = lLine[0]
        lTaxonomy = lLine[1:]
        dSeqId2Taxonomy[sSeqId] = lTaxonomy
        sGroup = 'none'
        for sTaxon in lTaxonomy[::-1]:
            if sTaxon in dTaxon2Group:
                sGroup = dTaxon2Group[sTaxon]
                break
        dSeqId2Group[sSeqId] = sGroup
else:
    from Taxonomy import Taxonomy
    taxonomy = Taxonomy(sTaxdumpDir)
    for sLine in open(sTaxIdMapping):
        (sSeqId, sTaxId) = sLine.strip().split()
        lTaxonomy = taxonomy.getTaxonomyByTaxId(sTaxId)
        if len(lTaxonomy) > 0:
            del lTaxonomy[-1] ## remove the 'Cellular organisms' etc. level
        dSeqId2Taxonomy[sSeqId] = map(lambda x:x[1], lTaxonomy)[::-1]
        sGroup = 'none'
        for (sTaxId, sTaxon, sTaxLevel) in lTaxonomy:
            if sTaxon in dTaxon2Group:
                sGroup = dTaxon2Group[sTaxon]
                break
        dSeqId2Group[sSeqId] = sGroup

## read the gene tree
tree = Tree(sTreeFile)

## root the tree by the query sequence and assign the group and a taxonomy to the leaf
## because we are now interested in the prokaryota-eukaryota LGT, will use the lowest 'eukaryota/bacteria/archaea' taxonomy, we will ignore viruses
## also count the prokaryotic groups
dProkGroups = {}
bQueryCheck = False
for leaf in tree.get_leaves():
    if leaf.name == sQueryId:
        tree.set_outgroup(leaf)
        bQueryCheck = True
    leaf.group = dSeqId2Group[leaf.name]
    ## eukaryota/prokaryota is totally sufficient for us
    ## set frequency to zero if the sequence has taxonomy outside eukaryota/prokaryota
    sTaxon = None
    fTaxFreq = 1
    if dSeqId2Taxonomy[leaf.name][0] == 'eukaryota':
        sTaxon = 'eukaryota'
    elif dSeqId2Taxonomy[leaf.name][0] in ['bacteria', 'archaea']:
        sTaxon = 'prokaryota'
        if leaf.group != 'none':
            dProkGroups[leaf.group] = None
    else:
        sTaxon = 'other'
        fTaxFreq = 0
    leaf.taxonomy = {sTaxon:fTaxFreq}
if not bQueryCheck:
    print 'Query id not found!'
    raise

## you can remove nodes with low support
def removeUnsupported(tree, fNodeSupportCutoff):
    for node in tree.traverse():
        if node == tree or node.is_leaf():
            continue
        if node.support < fNodeSupportCutoff:
            ancestor = node.get_ancestors()[0]
            for child in node.get_children():
                child.dist += node.dist
                ancestor.add_child(child)
            ancestor.remove_child(node)
    return tree

removeUnsupported(tree, fSupportCutoff)

fileOut.write('Input tree with unsupported nodes removed. In ascii format, having problems writing in newick:\n')
fileOut.write( tree.get_ascii()+'\n' )
fileOut.write('=====\n')
fileOut.flush()

## Progressively assign a compound taxonomy going from the leaves to the root.
def setTaxonomy(node):
    if node.is_leaf():
        return node
    node.taxonomy = {}
    lChildren = node.get_children()
    fChildrenTaxFreqSum = 0 ## typically the number of children, but some taxons may be masked
    for child in lChildren:
        setTaxonomy(child)
        for (sTaxon, fFreq) in child.taxonomy.items():
            fChildrenTaxFreqSum += fFreq
    for child in lChildren:
        for (sTaxon, fFreq) in child.taxonomy.items():
            if sTaxon not in node.taxonomy:
                node.taxonomy[sTaxon] = 0
            if fChildrenTaxFreqSum > 0:
                node.taxonomy[sTaxon] += float(fFreq)/fChildrenTaxFreqSum
    return node

setTaxonomy(tree)

## create pseudonodes in multifurcations connecting children with dominant 'eukaryota' or 'prokaryota' taxonomy respectively
def createPseudonodes(node):
    if node.is_leaf():
        return node
    for child in node.get_children():
        createPseudonodes(child)
    if len(node.get_children()) > 2:
        dDominantTaxon2Children = {}
        for child in node.get_children():
            sDominantTaxon = 'prokaryota'
            if 'eukaryota' in child.taxonomy and child.taxonomy['eukaryota'] >= 0.5:
                sDominantTaxon = 'eukaryota'
            if sDominantTaxon not in dDominantTaxon2Children:
                dDominantTaxon2Children[sDominantTaxon] = []
            dDominantTaxon2Children[sDominantTaxon].append(child)
        if len(dDominantTaxon2Children) > 1:
            for (sDominantTaxon, lDominantTaxonChildren) in dDominantTaxon2Children.items():
                if len(lDominantTaxonChildren) > 1:
                    newChild = Tree()
                    newChild.dist = min(map(lambda x:x.dist, node.get_children()))/2.
                    for child in lDominantTaxonChildren:
                        child.dist -= newChild.dist
                        newChild.add_child(child)
                        node.remove_child(child)
                    node.add_child(newChild)
    return node

createPseudonodes(tree)

fileOut.write('Input tree with some multifurcations resolved connecting branches with similar dominant taxonomy. In ascii format, having problems writing in newick:\n')
fileOut.write( tree.get_ascii()+'\n' )
fileOut.write('=====\n')
fileOut.flush()

## Progressively assign a compound taxonomy going from the leaves to the root again.
setTaxonomy(tree)

## For each eukaryotic leaf follow the ancestor and record the value of the 'eukaryota' taxonomy.
## Whenever the 'eukaryota' taxonomy value drops bellow 0.5, I considered it not to be a eukaryotic ortholog of the query sequence. You can however play with it and maybe set it lower.
## Count eukaryotic groups considered to be orthologous to query.
dSeqId2Orthology = {}
dEukGroups = {}
for leaf in tree.get_leaves():
    if 'eukaryota' in leaf.taxonomy:
        fMinEuk = 1
        for ancestor in leaf.get_ancestors():
            fEukFreq = 0
            if 'eukaryota' in ancestor.taxonomy:
                fEukFreq = ancestor.taxonomy['eukaryota']
            fMinEuk = min(fMinEuk, fEukFreq)
        if fMinEuk >= fOrthologyCutoff:
            dEukGroups[leaf.group] = None
        dSeqId2Orthology[leaf.name] = fMinEuk

fileOut.write('Orthology and Selected group:\n')
for (sSeqId, fOrthology) in sorted(dSeqId2Orthology.items(), key=lambda x:(-x[1], dSeqId2Group[x[0]])): ####
    fileOut.write('\t'+sSeqId+' '+str(fOrthology)+' '+dSeqId2Group[sSeqId]+'\n')
fileOut.write('=====\n')
fileOut.flush()

fDirectionalityScore = float(len(dProkGroups))/(len(dEukGroups)+len(dProkGroups))
fileOut.write('Deciding directionality of the transfer. (number of eukaryotic vs prokaryotic groups)\n')
fileOut.write('\tEukaryotic groups orthologous to the query: '+str(len(dEukGroups))+'\n')
fileOut.write('\tProkaryotic groups: '+str(len(dProkGroups))+'\n')
fileOut.write('\tDirectionality score: '+str(fDirectionalityScore)+'\n')
fileOut.write('=====\n')
fileOut.flush()

## Now that we have eukaryotic orthologs, we can map presence of the orthologs on the phylogeny of eukaryotes
for leaf in treeEukGroups.get_leaves():
    leaf.presence = 0
    if leaf.name in dEukGroups:
        leaf.presence = 1

## Count number of losses in case the gene was present in the last common ancestor of eukaryotes
## for each node count number of presences and losses

## count sum presence for each node
def countSumPresence(node):
    node.presenceSum = 0
    if node.is_leaf():
        if node.presence:
            node.presenceSum = 1 
    else:
        for child in node.get_children():
            countSumPresence(child)
            node.presenceSum += child.presenceSum

countSumPresence(treeEukGroups)

## find losses
def findLosses(treeGroups):
    dLeavesUsed = {}
    for leaf in treeGroups.get_leaves():
        if not leaf.presence and leaf.name not in dLeavesUsed:
            nodeLast = leaf
            for ancestor in leaf.get_ancestors():
                if ancestor.presenceSum > 0:
                    break
                nodeLast = ancestor
            for leafTemp in nodeLast.get_leaves():
                dLeavesUsed[leafTemp.name] = None
            nodeLast.loss = 1

findLosses(treeEukGroups)

## count sum losses for each node
def countSumLosses(node):
    node.lossesSum = 0
    if not hasattr(node, 'loss'):
        node.loss = 0
    if node.is_leaf():
        if node.loss == 1:
            node.lossesSum = 1
    else:
        for child in node.get_children():
            countSumLosses(child)
            node.lossesSum += child.lossesSum

countSumLosses(treeEukGroups)

## select possible gains
def selectGains(treeGroups):
    lGains = []
    for node in treeGroups.traverse():
        #node.gainScore = None
        if node.presenceSum > 0:
            if node.is_leaf():
                lGains.append(node)
            else:
                iCountChildrenWithPresence = 0
                for child in node.get_children():
                    if child.presenceSum >= 1:
                        iCountChildrenWithPresence += 1
                if iCountChildrenWithPresence >= 2:
                    lGains.append(node)
    return lGains

lGains = selectGains(treeEukGroups)

## find best combination of gains
def scoreGainCombos(lGains, iPresenceSum, iCutoffScore, lIndices=[], dCombo2Score={}):
    dPresences = {}
    bOverlapCheck = False
    iLossesSum = 0
    for iIndex in lIndices:
        node = lGains[iIndex]
        iLossesSum += node.lossesSum
        for leaf in node.get_leaves():
            if leaf.presence:
                if leaf.name in dPresences:
                    bOverlapCheck = True
                else:
                    dPresences[leaf.name] = None
        if bOverlapCheck:
            break
    ## stop adding gain to the combo if the number of events get bigger than in the case of single gain at the root of eukaryotes
    if not bOverlapCheck and iLossesSum+len(lIndices) <= iCutoffScore:
        if len(dPresences) == iPresenceSum:
            tCombo = tuple(sorted(map(lambda x:lGains[x].name, lIndices)))
            ## the score of gain combination is the number of events - gains and losses
            ## I set the to an equal weight.
            dCombo2Score[tCombo] = iLossesSum+len(lIndices)
        iStartIndex = 0
        if len(lIndices) > 0:
            iStartIndex = lIndices[-1]
        for iIndex in range(len(lGains)):
            if len(lIndices) == 0 or iIndex > lIndices[-1]:
                scoreGainCombos(lGains, iPresenceSum, iCutoffScore, lIndices+[iIndex], dCombo2Score)
    return dCombo2Score

dCombo2Score = scoreGainCombos(lGains, treeEukGroups.presenceSum, treeEukGroups.lossesSum+1)

if ('root',) in dCombo2Score:
    del dCombo2Score[('root',)]
if ('eukaryota',) in dCombo2Score:
    del dCombo2Score[('eukaryota',)]

lLossesAncestral = []
for node in treeEukGroups.traverse():
    if node.loss:
        lLossesAncestral.append( node.name )

tBestCombo = None
if len(dCombo2Score) > 0:
    tBestCombo = max( dCombo2Score.keys(), key=lambda x:(-dCombo2Score[x], len(x)) )

lLossesNew = []
for gain in treeEukGroups.traverse():
    if gain.name in tBestCombo:
        for node in gain.traverse():
            if node.loss:
                lLossesNew.append( node.name )

fNonacestralScore = 0
if tBestCombo != None:
    fNonacestralScore = float(treeEukGroups.lossesSum)/( treeEukGroups.lossesSum+len(lLossesNew)+len(tBestCombo)-1 )

fileOut.write('Comparing presence in the common ancestor of eukaryotes vs. later LGT(s):\n')
fileOut.write('\tCount of losses in case the gene was present in LECA: '+str(treeEukGroups.lossesSum)+'\n')
fileOut.write('\tLosses in case the gene was present in LECA: '+', '.join(sorted(lLossesAncestral))+'\n')
if tBestCombo != None:
    fileOut.write('\tCount of gains in case of LGT: '+str(len(tBestCombo))+'\n')
    fileOut.write('\tGains in case of LGT: '+', '.join(sorted(list(tBestCombo)))+'\n')
    fileOut.write('\tCount of losses in case of LGT: '+str(len(lLossesNew))+'\n')
    fileOut.write('\tLosses in case of LGT: '+', '.join(sorted(lLossesNew))+'\n')
else:
    fileOut.write('\tCount of gains in case of LGT: LGT scenario has always more events than ancestral eukaryotic\n')
    fileOut.write('\tGains in case of LGT: LGT scenario has always more events than ancestral eukaryotic\n')
    fileOut.write('\tCount of losses in case of LGT: LGT scenario has always more events than ancestral eukaryotic\n')
    fileOut.write('\tLosses in case of LGT: LGT scenario has always more events than ancestral eukaryotic\n')
fileOut.write('\tNonancestral score: '+str(fNonacestralScore)+'\n')
fileOut.write('=====\n')
fileOut.flush()
fileOut.write('If both scores are above 0.5, I think it is reasonable to assume that these may be prokaryota to eukaryota lateral gene transfers.\n')
fileOut.write('=====\n')
fileOut.flush()



fileOut.write('DONE\n')
fileOut.close()
