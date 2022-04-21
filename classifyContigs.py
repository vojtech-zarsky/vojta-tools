#!/usr/bin/env python3
# python ../myclasses/classifyContigs.py --fasta calonympha.transcripts.sel.fasta --blastn calonympha.transcripts.sel.fasta.nt.blastn --blastx calonympha.transcripts.sel.fasta.refprot.blastx.gz --taxdump ../myclasses/taxdump --library ../myclasses --blastn_eval 1e-20 --blastn_identity 0.8 --best_score_fraction 0.67 --threads 2 --selected_taxons "eukaryota;parabasalia;bacteroidetes;viruses"

## no weights among regions?


import os
import sys
import argparse
import gzip
import numpy as np
from multiprocessing import Process


parser = argparse.ArgumentParser(prog='Classify contigs based on BLASTX and/or BLASTN. DEPENDENCIES: SeqIO.py, ParseBlast.py, Taxonomy.py (https://github.com/vojtech-zarsky/vojta-tools)')

parser.add_argument('--fasta', required=True, help='Contigs in FASTA format.')
parser.add_argument('--blastn', required=False, help='BLASTN search result in standard tabular format + taxids (-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids")')
parser.add_argument('--blastx', required=False, help='BLASTX search result in standard tabular format + taxids (-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids")')
parser.add_argument('--taxdump', required=True, help='NCBI taxdump directory, can be downloaded here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz')
parser.add_argument('--library', required=False, default='./', help='Location of vojta-tools python library. [./]')
parser.add_argument('--blastn_eval', required=False, default='1e-20', help='BLASTN evalue cutoff. [1e-20]')
parser.add_argument('--blastn_identity', required=False, default='0.8', help='BLASTN identity cutoff. [0.8]')
parser.add_argument('--blastx_eval', required=False, default='1e-10', help='BLASTX evalue cutoff. [1e-10]')
parser.add_argument('--best_score_fraction', required=False, default='0.67', help='Evaluate hits with this fraction of overlapping best hit score. [0.67]')
parser.add_argument('--tax_levels', required=False, default='6', help='How many taxonomy levels to use. [6]')
parser.add_argument('--threads', required=False, default='1', help='How many threads. [1]')
parser.add_argument('--selected_taxons', required=False, help='Taxons to select. ("tax1;tax2;...")')
args = parser.parse_args()

args.blastn_eval = float(args.blastn_eval)
args.blastn_identity = float(args.blastn_identity)
if args.blastn_identity < 0 or args.blastn_identity > 1:
    print('--blastn_identity must be in the [0-1] range!')
    raise
args.blastx_eval = float(args.blastx_eval)
args.best_score_fraction = float(args.best_score_fraction)
if args.best_score_fraction < 0 or args.best_score_fraction > 1:
    print('--best_score_fraction must be in the [0-1] range!')
    raise
args.tax_levels = int(args.tax_levels)
args.threads = int(args.threads)
if args.selected_taxons == None:
    args.selected_taxons = []
else:
    args.selected_taxons = list(map(lambda x:x.strip(), args.selected_taxons.split(';')))

sys.path.append(args.library)
from SeqIO import SeqIO
from ParseBlast import ParseBlast
from Taxonomy import Taxonomy

## get gc-content
def getGC(sSeq):
    iGC = 0
    iAT = 0
    for sNt in seq.seq.upper():
        if sNt in ['G','C']:
            iGC += 1
        elif sNt in ['A','T']:
            iAT += 1
    if iGC+iAT > 0:
        return iGC/(iGC+iAT)
    else:
        return None

def getOverlap(lRange1, lRange2):
    lRange1 = sorted(list(lRange1))
    lRange2 = sorted(list(lRange2))
    return min(lRange1[1], lRange2[1])-max(lRange1[0], lRange2[0])+1

def parseRegionsTable(sTable, sType, dContigs, iTaxLevels, fEvalCutoff, fIdentityCutoff):
    for sContig in dContigs.keys():
        dContigs[sContig][sType] = {'lcaTaxonomies':[], 'lcaTaxonomiesLine':['']*(iTaxLevels*2), 'allTaxonomies':[], 'allTaxonomiesLine':['']*(iTaxLevels*2+1), 'bestHit':['', '','','','','','','']}

    dContigsTemp = {}
    lColumns = None
    for sLine in open(sTable):
        lLine = sLine.strip('\n').split('\t')
        if lColumns == None:
            lColumns = lLine
            continue
        if lLine == lColumns:
            continue
        dLine = dict(zip(lColumns, lLine))
        dLine['BH score'] = float(dLine['BH score'])
        dLine['BH eval'] = float(dLine['BH eval'])
        dLine['BH identity'] = float(dLine['BH identity'])
        if dLine['BH eval'] > fEvalCutoff or dLine['BH identity'] < fIdentityCutoff:
            continue
        sContig = dLine['contig']
        print(sContig)
        if sContig not in dContigsTemp:
            dContigsTemp[sContig] = {'lcaTaxonomies':[]}
            for iTaxLevel in range(iTaxLevels):
                dContigsTemp[sContig]['lcaTaxonomies'].append({})
        if dContigs[sContig][sType]['bestHit'][4] == '' or dLine['BH score'] > dContigs[sContig][sType]['bestHit'][4]:
            dContigs[sContig][sType]['bestHit'] = [dLine['BH id'], dLine['BH taxid'], dLine['BH species'], dLine['BH eval'], dLine['BH score'], dLine['BH identity'], str(dLine['start'])+'-'+str(dLine['end']), dContigs[sContig]['seq'][int(dLine['start'])-1:int(dLine['end'])]]
        ## read LCA taxonomies
        lLcaTaxonomy = dLine['LCA taxonomy'].strip().split(';')
        for iTaxLevel in range(min(iTaxLevels, len(lLcaTaxonomy))):
            sTaxon = lLcaTaxonomy[iTaxLevel]
            if sTaxon not in dContigsTemp[sContig]['lcaTaxonomies'][iTaxLevel]:
                dContigsTemp[sContig]['lcaTaxonomies'][iTaxLevel][sTaxon] = 0
            dContigsTemp[sContig]['lcaTaxonomies'][iTaxLevel][sTaxon] += 1 ##dLine['BH score'] -> don't weight scores among regions
        ## read all taxonomies
        lTaxonFrequencies = []
        for sTaxLevelFrequencies in dLine['taxon frequencies'].strip().split(';'):
            dTemp = {}
            for sTaxonFrequency in sTaxLevelFrequencies.split(','):
                lTaxonFrequency = sTaxonFrequency.rsplit(':',1)
                if len(lTaxonFrequency) == 2:
                    dTemp[lTaxonFrequency[0]] = float(lTaxonFrequency[1])
            if dTemp != {}:
                lTaxonFrequencies.append(dTemp)
        for iTaxLevel in range(len(lTaxonFrequencies)):
            if len(dContigs[sContig][sType]['allTaxonomies']) <= iTaxLevel:
                dContigs[sContig][sType]['allTaxonomies'].append({})
            for sTaxon, fFrequency in lTaxonFrequencies[iTaxLevel].items():
                if sTaxon not in dContigs[sContig][sType]['allTaxonomies'][iTaxLevel]:
                    dContigs[sContig][sType]['allTaxonomies'][iTaxLevel][sTaxon] = 0
                dContigs[sContig][sType]['allTaxonomies'][iTaxLevel][sTaxon] += fFrequency
        ## track species
    ## for each contig summarize LCA taxonomies
    for sContig in dContigsTemp.keys():
        fMinScore = 1
        for iTaxLevel in range(iTaxLevels):
            for sTaxon, fScore in sorted(dContigsTemp[sContig]['lcaTaxonomies'][iTaxLevel].items(), key=lambda x:(-x[1], x[0])):
                if iTaxLevel >= len(dContigs[sContig][sType]['lcaTaxonomies']):
                    dContigs[sContig][sType]['lcaTaxonomies'].append({})
                dContigs[sContig][sType]['lcaTaxonomies'][iTaxLevel][sTaxon] = fScore#/fScoreSum
                if dContigs[sContig][sType]['lcaTaxonomiesLine'][iTaxLevel*2] == '':
                    dContigs[sContig][sType]['lcaTaxonomiesLine'][iTaxLevel*2] = sTaxon
                    fMinScore = min(fMinScore, fScore/sum(dContigsTemp[sContig]['lcaTaxonomies'][iTaxLevel].values()))
                    dContigs[sContig][sType]['lcaTaxonomiesLine'][iTaxLevel*2+1] = fMinScore
                #break
    ## all
    for sContig in dContigs.keys():
        ## set sum to min score
        fMinScore = 1
        lAllTaxonomiesLine = []
        for iTaxLevel in range(len(dContigs[sContig][sType]['allTaxonomies'])):
            lTemp = []
            for sTaxon, fScore in sorted(dContigs[sContig][sType]['allTaxonomies'][iTaxLevel].items(), key=lambda x:(-x[1], x[0])):
                lTemp.append('{}:{}'.format(sTaxon, fScore))
                if dContigs[sContig][sType]['allTaxonomiesLine'][iTaxLevel*2] == '':
                    dContigs[sContig][sType]['allTaxonomiesLine'][iTaxLevel*2] = sTaxon
                    fMinScore = min(fMinScore, fScore/sum(dContigs[sContig][sType]['allTaxonomies'][iTaxLevel].values()))
                    dContigs[sContig][sType]['allTaxonomiesLine'][iTaxLevel*2+1] = fMinScore
            lAllTaxonomiesLine.append(','.join(lTemp))
        dContigs[sContig][sType]['allTaxonomiesLine'][-1] = ';'.join(lAllTaxonomiesLine)

def getN50(lLengths):
    fSum = sum(lLengths)
    fCountSum = 0
    for fLen in sorted(lLengths):
        fCountSum += fLen
        if fCountSum >= fSum/2:
            return fLen
    return 0



if __name__ == "__main__":
    ## read contigs
    dContigs = {}
    lContigs = []
    for seq in SeqIO(args.fasta):
        lContigs.append(seq.id)
        fCoverage = None
        if '_cov_' in seq.id: ## SPADES
            fCoverage = float(seq.id.split('_cov_')[1].split('_')[0])
        elif ' cov_' in seq.desc: ## trinity
            fCoverage = float(seq.desc.split(' cov_')[1].split()[0])
        elif 'multi=' in seq.desc: ## MEGAHIT
            fCoverage = float(seq.desc.split('multi=')[1].split()[0])
        dContigs[seq.id] = {'id':seq.id, 'desc':seq.desc, 'seq':seq.seq, 'len':len(seq.seq), 'gc':getGC(seq.seq), 'coverage':fCoverage, 'blastx':{}, 'blastn':{}}
    #taxonomy = Taxonomy(sTaxdumpDir=args.taxdump, bAddSupergroups=True, bCleanTaxons=True)

    parseRegionsTable('{}.parseBlasts.BlastnRegions.tsv'.format(args.fasta), 'blastn', dContigs, args.tax_levels, args.blastn_eval, args.blastn_identity)
    parseRegionsTable('{}.parseBlasts.BlastxRegions.tsv'.format(args.fasta), 'blastx', dContigs, args.tax_levels, args.blastx_eval, 0)
    ## summarize blastn+blastx
    for sContig in dContigs.keys():
        dContigs[sContig]['combined'] = {'lcaTaxonomies':[], 'lcaTaxonomiesLine':['']*(args.tax_levels*2), 'allTaxonomies':[], 'allTaxonomiesLine':['']*(args.tax_levels*2+1)}
        for sType in ['blastx','blastn']:
            ## lca
            for iTaxLevel in range(len(dContigs[sContig][sType]['lcaTaxonomies'])):
                if iTaxLevel >= len(dContigs[sContig]['combined']['lcaTaxonomies']):
                    dContigs[sContig]['combined']['lcaTaxonomies'].append({})
                for sTaxon, fScore in dContigs[sContig][sType]['lcaTaxonomies'][iTaxLevel].items():
                    if sTaxon not in dContigs[sContig]['combined']['lcaTaxonomies'][iTaxLevel]:
                        dContigs[sContig]['combined']['lcaTaxonomies'][iTaxLevel][sTaxon] = 0
                    dContigs[sContig]['combined']['lcaTaxonomies'][iTaxLevel][sTaxon] += fScore
            ## all
            for iTaxLevel in range(len(dContigs[sContig][sType]['allTaxonomies'])):
                if iTaxLevel >= len(dContigs[sContig]['combined']['allTaxonomies']):
                    dContigs[sContig]['combined']['allTaxonomies'].append({})
                for sTaxon, fScore in dContigs[sContig][sType]['allTaxonomies'][iTaxLevel].items():
                    if sTaxon not in dContigs[sContig]['combined']['allTaxonomies'][iTaxLevel]:
                        dContigs[sContig]['combined']['allTaxonomies'][iTaxLevel][sTaxon] = 0
                    dContigs[sContig]['combined']['allTaxonomies'][iTaxLevel][sTaxon] += fScore
        ## lca
        fMinScore = 1
        for iTaxLevel in range(len(dContigs[sContig]['combined']['lcaTaxonomies'])):
            for sTaxon, fScore in sorted(dContigs[sContig]['combined']['lcaTaxonomies'][iTaxLevel].items(), key=lambda x:(-x[1], x[0])):
                dContigs[sContig]['combined']['lcaTaxonomiesLine'][iTaxLevel*2] = sTaxon
                fMinScore = min(fMinScore, fScore/sum(dContigs[sContig]['combined']['lcaTaxonomies'][iTaxLevel].values()))
                dContigs[sContig]['combined']['lcaTaxonomiesLine'][iTaxLevel*2+1] = fMinScore
                break
        ## all
        fMinScore = 1
        lOutLine = []
        for iTaxLevel in range(len(dContigs[sContig]['combined']['allTaxonomies'])):
            lTemp = []
            for sTaxon, fScore in sorted(dContigs[sContig]['combined']['allTaxonomies'][iTaxLevel].items(), key=lambda x:(-x[1], x[0])):
                lTemp.append( '{}:{}'.format(sTaxon, fScore) )
                if dContigs[sContig]['combined']['allTaxonomiesLine'][iTaxLevel*2] == '':
                    dContigs[sContig]['combined']['allTaxonomiesLine'][iTaxLevel*2] = sTaxon
                    fMinScore = min(fMinScore, fScore/sum(dContigs[sContig]['combined']['allTaxonomies'][iTaxLevel].values()))
                    dContigs[sContig]['combined']['allTaxonomiesLine'][iTaxLevel*2+1] = fMinScore
            lOutLine.append(','.join(lTemp))
        dContigs[sContig]['combined']['allTaxonomiesLine'][-1] = ';'.join(lOutLine)
    #print(dContigs['NODE_95_length_3736_cov_0.237009_g42_i0'])

    ## write tables
    ## headers
    tableOut = open('{}.classifyContigs.table.tsv'.format(args.fasta),'w')
    tableOut.write('contig\tlen\tcoverage\tgc')
    for sSelectedTaxon in args.selected_taxons:
        tableOut.write('\t{} x+n freqs'.format(sSelectedTaxon))
    tableOut.write('\txbh id\txbh taxid\txbh organism\txbh eval\txbh score\txbh identity\txbh qrange\txbh qseq')
    tableOut.write('\tnbh id\tnbh taxid\tnbh organism\tnbh eval\tnbh score\tnbh identity\tnbh qrange\tnbh qseq')
    for iTaxLevel in range(args.tax_levels):
        tableOut.write('\tx+n lca tax{0}\tx+n lca tax{0} score'.format(iTaxLevel))
    for iTaxLevel in range(args.tax_levels):
        tableOut.write('\tx+n tax{0}\tx+n tax{0} score'.format(iTaxLevel))
    tableOut.write('\tx+n taxon frequencies')
    for iTaxLevel in range(args.tax_levels):
        tableOut.write('\tx lca tax{0}\tx lca tax{0} score'.format(iTaxLevel))
    for iTaxLevel in range(args.tax_levels):
        tableOut.write('\tx tax{0}\tx tax{0} score'.format(iTaxLevel))
    tableOut.write('\tx taxon frequencies')
    for iTaxLevel in range(args.tax_levels):
        tableOut.write('\tn lca tax{0}\tn lca tax{0} score'.format(iTaxLevel))
    for iTaxLevel in range(args.tax_levels):
        tableOut.write('\tn tax{0}\tn tax{0} score'.format(iTaxLevel))
    tableOut.write('\tn taxon frequencies')
    tableOut.write('\n')
    ## content
    for sContig in lContigs:
        dContig = dContigs[sContig]
        lLineOut = [sContig, dContig['len'], dContig['coverage'], dContig['gc']]
        for sSelectedTaxon in args.selected_taxons:
            fScoreOut = 0
            for iTaxLevel in range(len(dContigs[sContig]['combined']['allTaxonomies'])):
                for sTaxon, fScore in dContigs[sContig]['combined']['allTaxonomies'][iTaxLevel].items():
                    if sTaxon == sSelectedTaxon:
                        fScoreOut = fScore/sum(dContigs[sContig]['combined']['allTaxonomies'][0].values())
                        break
                if fScoreOut > 0:
                    break
            lLineOut += [fScoreOut]
        lLineOut += dContig['blastx']['bestHit']
        lLineOut += dContig['blastn']['bestHit']
        lLineOut += dContig['combined']['lcaTaxonomiesLine']
        lLineOut += dContig['combined']['allTaxonomiesLine']
        lLineOut += dContig['blastx']['lcaTaxonomiesLine']
        lLineOut += dContig['blastx']['allTaxonomiesLine']
        lLineOut += dContig['blastn']['lcaTaxonomiesLine']
        lLineOut += dContig['blastn']['allTaxonomiesLine']
        tableOut.write('\t'.join(list(map(lambda x:str(x), lLineOut)))+'\n')
    tableOut.close()

    ## write selected taxons
    for sSelectedTaxon in args.selected_taxons:
        fileOut = open('{}.classifyContigs.{}.fasta'.format(args.fasta, sSelectedTaxon),'w')
        fileNotOut = open('{}.classifyContigs.not-{}.fasta'.format(args.fasta, sSelectedTaxon),'w')
        ## select by overall frequency ("all")
        for sContig in lContigs:
            dContig = dContigs[sContig]
            for iTaxLevel in range(len(dContig['combined']['allTaxonomies'])):
                for sTaxon, fScore in dContig['combined']['allTaxonomies'][iTaxLevel].items():
                    if fScore > 0.5 and sTaxon.lower() == sSelectedTaxon.lower():
                        fileOut.write('>{} {}\n{}\n'.format(sContig, dContig['desc'], dContig['seq']))
                    else:
                        fileNotOut.write('>{} {}\n{}\n'.format(sContig, dContig['desc'], dContig['seq']))
        fileOut.close()
        fileNotOut.close()
    
    ## summarize taxons
    lTaxLevel2Taxons = []
    for iTaxLevel in range(args.tax_levels):
        lTaxLevel2Taxons.append({'unclassified':{'contigs':[], 'lens':[], 'coverages':[], 'gcs':[]}})
    lColumns = None
    for sLine in open('{}.classifyContigs.table.tsv'.format(args.fasta)):
        lLine = sLine.strip('\n').split('\t')
        if lColumns == None:
            lColumns = lLine
            continue
        dLine = dict(zip(lColumns, lLine))
        for iTaxLevel in range(args.tax_levels):
            sTaxon = dLine['x+n tax{}'.format(iTaxLevel)]
            if sTaxon == '' or float(dLine['x+n tax{} score'.format(iTaxLevel)]) < 2/3: ##
                sTaxon = 'unclassified'
            if sTaxon not in lTaxLevel2Taxons[iTaxLevel]:
                lTaxLevel2Taxons[iTaxLevel][sTaxon] = {'contigs':[], 'lens':[], 'coverages':[], 'gcs':[]}
            lTaxLevel2Taxons[iTaxLevel][sTaxon]['contigs'].append(dLine['contig'])
            lTaxLevel2Taxons[iTaxLevel][sTaxon]['lens'].append(float(dLine['len']))
            lTaxLevel2Taxons[iTaxLevel][sTaxon]['coverages'].append(float(dLine['coverage']))
            lTaxLevel2Taxons[iTaxLevel][sTaxon]['gcs'].append(float(dLine['gc']))
    tableOut = open('{}.classifyContigs.taxons.tsv'.format(args.fasta),'w')
    tableOut.write('tax level\ttaxon\tcontigs\tsum len\tN50 len\tavg len\t10 len\t90 len\tavg coverage\t10 coverage\t90 coverage\tavg gc\t10 gc\t90 gc\n')
    for iTaxLevel in range(args.tax_levels):
        for sTaxon, dTaxon in sorted(lTaxLevel2Taxons[iTaxLevel].items(), key=lambda x:-sum(x[1]['lens'])):
            lLineOut = [iTaxLevel, sTaxon, len(dTaxon['contigs'])]
            ## lengths
            lLineOut.append(sum(dTaxon['lens']))
            lLineOut.append(getN50(dTaxon['lens']))
            array = np.array(dTaxon['lens'])
            lLineOut.append( np.mean(array) )
            lLineOut.append( np.percentile(array, 10) )
            lLineOut.append( np.percentile(array, 90) )
            ## coverages
            array = np.array(dTaxon['coverages'])
            lLineOut.append( np.mean(array) )
            lLineOut.append( np.percentile(array, 10) )
            lLineOut.append( np.percentile(array, 90) )
            ## gcs
            array = np.array(dTaxon['gcs'])
            lLineOut.append( np.mean(array) )
            lLineOut.append( np.percentile(array, 10) )
            lLineOut.append( np.percentile(array, 90) )
            tableOut.write('\t'.join(list(map(lambda x:str(x), lLineOut)))+'\n')
    tableOut.close()






