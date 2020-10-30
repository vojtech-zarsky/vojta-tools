## for python3
import os
import re
import gzip
import sqlite3
import numpy
from multiprocessing import Process

"""
MAIN FUNCTIONALITY:

- evaluate reads (fastqc)
- trim reads
- assemble
- map reads to the contigs
- search contigs against protein database
- assign taxonomies to the contigs, detect chimeric contigs
- get SSU rRNAs, make phylogeny, evaluate phylogeny -> special module for this?
- evaluate using BlobTools
- evaluate using BUSCO
- assign genetic codes and evaluate introns for contigs and taxons
- predict protein-coding genes and RNA genes
- deal with closely related organisms?

========================================

FOLDER STRUCTRE:



========================================

TAXONOMIC ASSIGNMENT:
(1) Identity based:
- preselect using fixed fraction (1/2 - 2/3) of the score of the best e-val hit as cutoff
- sort by identities
- select those that have identity better than bestHitIdentity**2
- identityNorm is computed as identity/sum(selected identities)
- PROBLEM: (case Coronympha k141_100178)
- the true closest hit may have a longer alignment and could then easily have lower overall identity than others
- GOOD FOR:
- taxonomic assignments of closely related sequences (identity > 70%)

(2) bit-score based:
- preselect using fixed fraction (1/2 - 2/3) of the score of the best e-val hit as cutoff
- sort by bit-score
- select those that have identity better than bestHitIdentity**2
- scoreNorm is computed as bit-score/sum(selected bit-scores)

REMAINING PROBLEMS:
- taxonomic assignment gets weak quickly with more distant hits
- contaminations in the target database

BETTER SOLUTION:
- tree-based taxonomic assignment
"""

class SeqProject:
    def __init__(self, sName, sContigFile, lFRReads, sDb, sDbDir, bUseSql):
        from Genome import Genome
        from ParseCm import ParseCm
        from SeqIO import SeqIO
        from ParseBlast import ParseBlast
        from Taxonomy import Taxonomy
        self.sName = sName
        self.sContigFile = sContigFile
        self.lFRReads = lFRReads
        self.sDb = sDb
        self.sDbDir = sDbDir
        self.bUseSql = bUseSql
        self.taxonomy = False
        print('initializing  seqproj name=%s db=%s'%(self.sName, self.sDb))
        ## initialize taxonomy assignment
        if self.bUseSql:
            self.connectionMain = sqlite3.connect('%s/%s.sqlite3'%(self.sDbDir, self.sDb))
            self.cursorMain = self.connectionMain.cursor()
        else:
            print('reading taxonomic assignment..')
            lPairs = []
            for sLine in open('/Data/vojta/uniprot/%s.taxids'%(self.sDb)):
                lPairs.append( sLine.strip().split('\t') )
            self.dUniprotAccession2TaxId = dict(lPairs)
        ## check/create status file
        ## check/create project folder
    def getTaxId(self, sUniprotId, cursor=None):
        if bUseSql:
            if cursor == None:
                cursor = self.cursorMain
            cursor.execute('SELECT taxid FROM seqs WHERE id=\'%s\''%(sUniprotId))
            return str(cursor.fetchone()[0])
        else:
            return self.dUniprotAccession2TaxId[sUniprotId]
    def getSsurnaModels(self, sRfamDatabase):
        self.dDomain2Rfam = {'bacteria':'SSU_rRNA_bacteria', 'archaea':'SSU_rRNA_archaea', 'eukaryota':'SSU_rRNA_eukarya'}
        fileOut = open('rfam.ssu_rrnas.cm','w')
        for sModel in open('../databases/rfam/Rfam.14_3.cm').read().split('\n//')[:-1]:
            sModel = sModel.strip()
            sName = sModel.split('NAME',1)[1].strip().split()[0]
            if sName in dDomain2Rfam.values():
                fileOut.write('%s\n//\n'%(sModel))
        fileOut.close()
    def scanForRRNA(self):
        # !!!!!!!!!!!!!!
        lHits = []
        #for (sDomain, sRfam) in dDomain2Rfam.items():
        #    dRfam2File[sRfam].close()
        #    sRfamModel = dDomain2Rfam[sDomain]
        #    sFileName = '%s.%s'%(sSequenceData.split('/')[-1], sRfamModel)
        #    if sAdditionalName != None:
        #        sFileName = '%s.%s.%s'%(sAdditionalName, sSequenceData.split('/')[-1], sRfamModel)
        #    os.system('cmsearch --cpu 4 --tblout %s.tblout  -o %s.o --notextw %s.rfam_14_3.cm %s'%(sFileName,  sFileName, sRfamModel, sSequenceData))
        #    lHits += ParseCm(open('%s.tblout'%(sFileName))).next()
        #fileOut.close()
        #os.system('cmpress -F rfam.ssu_rrnas.cm')
        #os.system('cmscan --cpu 4 --notextw rfam.ssu_rrnas.cm %s'%(sSequenceData))
        #dSeqs = SeqIO(open(sSequenceData)).getDict()
        genome = Genome(sSequenceData)
        lHits = sorted(lHits, key=lambda x:(-x.score, x.eval))
        fileOut = open('%s.rrnas.fasta'%(sAdditionalName),'w')
        dSeq2RangesUsed = {}
        for hit in lHits:
            if hit.eval > 1e-3:
                continue
            iStart, iEnd = hit.hitStart-1, hit.hitEnd
            sSeq = genome.seqs[genome.dSeqId2SeqIndex[hit.hitId]].seq[iStart:iEnd]
            if hit.strand == '-':
                iStart, iEnd = hit.hitEnd-1, hit.hitStart
                sSeq = genome.reverseComplement(genome.seqs[genome.dSeqId2SeqIndex[hit.hitId]].seq[iStart:iEnd])
            bOverlapCheck = True
            if hit.hitId in dSeq2RangesUsed:
                for (iStart2, iEnd2) in dSeq2RangesUsed[hit.hitId]:
                    iOverlap = min(iEnd, iEnd2)-max(iStart, iStart2)+1
                    if iOverlap >= 1:
                        bOverlapCheck = False
                        break
            if bOverlapCheck:
                if hit.hitId not in dSeq2RangesUsed:
                    dSeq2RangesUsed[hit.hitId] = []
                dSeq2RangesUsed[hit.hitId].append((iStart, iEnd))
                fileOut.write('>%s.%s.%s.%i.%i %f %s %f\n%s\n'%(sAdditionalName, hit.queryId, hit.hitId, hit.hitStart, hit.hitEnd, hit.eval, hit.strand, hit.gc, sSeq))
        fileOut.close()
    def analyzeRRNA(self, sFile, dSilvaId2Taxonomy=None):
        ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #os.system('blastn -query %s -db /home/zarsky/databases/silva/SILVA_138.1_SSURef_tax_silva.fasta -out %s.silva_138.tab -outfmt 6 -num_threads 4 -task megablast'%(sFile, sFile))
        dSeqs = SeqIO(open('%s.cdhit100.fasta'%(sFile.rsplit('.',1)[0]))).getDict()
        fastaOut = open('%s.taxonomy.fasta'%(sFile.rsplit('.')[0]),'w')
        for lHits in ParseBlast(open('%s.silva_138.tab'%(sFile))):
            sQueryId = lHits[0].queryId
            if sQueryId not in dSeqs:
                continue
            #del dSeqs[sQueryId]
            lHits = sorted(lHits, key=lambda x:(-x.identity, -x.score))
            print(sQueryId, lHits[0].hitId, lHits[0].identity)
            lTaxonomy = dSilvaId2Taxonomy[lHits[0].hitId].strip().replace(' ','_').split(';')
            for iTax in range(len(lTaxonomy)):
                lTaxonomy[iTax] = re.sub('\W','',lTaxonomy[iTax])
            #or iTax in range(len(lTaxonomy)):
            #    lTaxonomy[iTax] = re.split(lTaxonomy[iTax]
            sId = sQueryId+'.'+str(int(lHits[0].identity))+'.'+'_'.join(lTaxonomy)
            #sId = re.sub('\W','',sId)
            fastaOut.write('>'+sId+'\n'+dSeqs[sQueryId].seq+'\n')
        fastaOut.close()
        #print dSeqs.keys()
    def mapReads(self, sGenomeFileDir=None, lReadFileDirs=None, sAdditionalName=None):
        #lRoots = map(lambda x:x.rsplit('.',2)[0], lReadFileDirs)
        """
        print('indexing...')
        os.system('bwa index %s'%(sGenomeFileDir))
        ## paired
        os.system('bwa mem -t 8 %s %s.paired.fastq.gz %s.paired.fastq.gz > %s.paired.sam'%(sGenomeFileDir, lRoots[0], lRoots[1], sGenomeFileDir))
        os.system('samtools view -S -b -o %s.paired.bam %s.paired.sam'%(sGenomeFileDir, sGenomeFileDir))
        os.system('rm -f %s.paired.sam'%(sGenomeFileDir))
        ## unpaired 1
        os.system('bwa mem -t 8 %s %s.unpaired.fastq.gz > %s.unpaired1.sam'%(sGenomeFileDir, lRoots[0], sGenomeFileDir))
        os.system('samtools view -S -b -o %s.unpaired1.bam %s.unpaired1.sam'%(sGenomeFileDir, sGenomeFileDir))
        os.system('rm -f %s.unpaired1.sam'%(sGenomeFileDir))
        ## unpaired 2
        os.system('bwa mem -t 8 %s %s.unpaired.fastq.gz > %s.unpaired2.sam'%(sGenomeFileDir, lRoots[1], sGenomeFileDir))
        os.system('samtools view -S -b -o %s.unpaired2.bam %s.unpaired2.sam'%(sGenomeFileDir, sGenomeFileDir))
        os.system('rm -f %s.unpaired2.sam'%(sGenomeFileDir))
        #"""
        #os.system('samtools sort -@ 8 -O bam -T temp -o %s.paired_sorted.bam %s.paired.bam'%(sGenomeFileDir, sGenomeFileDir))
        #os.system('samtools sort -@ 8 -O bam -T temp -o %s.unpaired1_sorted.bam %s.unpaired1.bam'%(sGenomeFileDir, sGenomeFileDir))
        #os.system('samtools sort -@ 8 -O bam -T temp -o %s.unpaired2_sorted.bam %s.unpaired2.bam'%(sGenomeFileDir, sGenomeFileDir))
        #os.system('samtools index %s.paired_sorted.bam'%(sGenomeFileDir))
        #os.system('samtools index %s.unpaired1_sorted.bam'%(sGenomeFileDir))
        #os.system('samtools index %s.unpaired2_sorted.bam'%(sGenomeFileDir))
        #lReadTypes = ['paired','unpaired1','unpaired2']
        #for sReadType in lReadTypes:
        #    os.system('rm %s.%s.bam'%(sGenomeFileDir, sReadType))
        ## get depth
        os.system('samtools depth %s.paired_sorted.bam %s.unpaired1_sorted.bam %s.unpaired2_sorted.bam > %s.depth'%(self.sContigFile, self.sContigFile, self.sContigFile, self.sContigFile))
    def blobTools(self, sGenomeFileDir, sOrgn):
        ## fro BlobTools ver. 0.9.8
        ## create taxfile
        """
        ## PROBLEM WHEN RUNNING WITH PYTHON 3!!!!
        tableOut = open('%s.%s.fast.blobtoolstab'%(sGenomeFileDir, sDb),'w')
        for lHits in ParseBlast(open('%s.%s.fast.tab'%(sGenomeFileDir, sDb))):
            lHits = sorted(lHits, key=lambda x:(-x.identity, x.eval, -x.score))
            for hit in lHits:
                if hit.eval <= 1e-5: ##
                    sTaxId = dUniprotAccession2TaxId[hit.hitId.split('|')[2]]
                    lTaxonomy = taxonomy.getTaxonomyByTaxId(sTaxId)[::-1]
                    if len(lTaxonomy) < 3 or (lTaxonomy[0][1] != 'viruses' and lTaxonomy[1][1] not in ('bacteria','archaea','eukaryota')):
                        print(hit.hitId.split('|')[2], sTaxId, lTaxonomy)
                        continue
                    tableOut.write('%s\t%s\t%f\n'%(hit.queryId, sTaxId, hit.score)) ## 4932
                    break
        tableOut.close()
        #"""
        #print('blobtools create -i %s -t %s.%s.fast.blobtoolstab --nodes /home/zarsky/databases/ncbi/taxdump.201002/nodes.dmp --names /home/zarsky/databases/ncbi/taxdump.201002/names.dmp -b %s.paired_sorted.bam -b %s.unpaired1_sorted.bam -b %s.unpaired2_sorted.bam -o %s'%(sGenomeFileDir, sGenomeFileDir, sDb, sGenomeFileDir, sGenomeFileDir, sGenomeFileDir, sOrgn))
        sTaxRank = 'species'
        for sTaxRank in ['superkingdom', 'phylum', 'order', 'family', 'genus', 'species']:
            os.system('blobtools view --rank %s -i %s/%s.blobDB.json --out %s'%(sTaxRank, sGenomeFileDir.rsplit('/',1)[0], sOrgn, sTaxRank))
        for sTaxRank in ['superkingdom', 'phylum', 'order']:
            os.system('blobtools plot --format png --lib covsum --rank %s --plotgroups 8 -i %s/%s.blobDB.json --out ./'%(sTaxRank, sGenomeFileDir.rsplit('/',1)[0], sOrgn)) ## species, genus, family, order, phylum, superkingdom
    def sortContigsFast(sGenomeFileDir, sOrgn, lTargetTaxons):
        dTargetTaxon2Contigs = {}
        for (sTargetTaxon, iTargetTaxonLevel) in lTargetTaxons:
            dTargetTaxon2Contigs[sTargetTaxon] = []
        dContigs = SeqIO(open(sGenomeFileDir)).getDict()
        for sLine in open('%s.%s.fast.blobtoolstab'%(sGenomeFileDir, sDb)):
            (sContig, sTaxId, sScore) = sLine.strip().split()
            lTaxonomy = map(lambda x:x[1], taxonomy.getTaxonomyByTaxId(sTaxId))[::-1]
            for (sTargetTaxon, iTargetTaxonLevel) in lTargetTaxons:
                if len(lTaxonomy) > iTargetTaxonLevel and lTaxonomy[iTargetTaxonLevel] == sTargetTaxon:
                    dTargetTaxon2Contigs[sTargetTaxon].append(dContigs[sContig])
                    del dContigs[sContig]
                    break
        if not os.path.isdir(sGenomeFileDir.rsplit('/',1)[0]+'/target_taxons/'):
            os.system('mkdir '+sGenomeFileDir.rsplit('/',1)[0]+'/target_taxons/')
        for (sTargetTaxon, iTargetTaxonLevel) in lTargetTaxons:
            fileOut = open('%s/target_taxons/%s.fasta'%(sGenomeFileDir.rsplit('/',1)[0], sTargetTaxon),'w')
            for seq in dTargetTaxon2Contigs[sTargetTaxon]:
                fileOut.write('>'+seq.id+' '+seq.desc+'\n'+seq.seq+'\n')
            fileOut.close()
    def parseSearchProcess(self, iThread, iThreads, sGenomeFileDir, sOrgn, fEvalCutoff, fBestHitScoreCutoff, iMaxHits):
        print(iThread)
        fileRanges = open('%s.%s.long.tab.ranges.thread%i'%(sGenomeFileDir, sDb, iThread),'w')
        fileSelectedHits = open('%s.%s.long.tab.selected.thread%i'%(sGenomeFileDir, sDb, iThread),'w')
        iHitsIndex = -1
        for lHits in ParseBlast(open('%s.%s.long.tab'%(sGenomeFileDir, sDb)), mergeHsps=True):
            iHitsIndex += 1
            if (iHitsIndex-iThread)%iThreads != 0:
                continue
            sContig = lHits[0][0].queryId
            dRanges2Hits = {}
            print(iThread, sContig, len(lHits))
            for lHsps in lHits:
                iStart = min(list(map(lambda x:x.queryStart, lHsps))+list(map(lambda x:x.queryEnd, lHsps)))
                iEnd = max(list(map(lambda x:x.queryStart, lHsps))+list(map(lambda x:x.queryEnd, lHsps)))
                iBestOverlap = 0
                tBestOverlapRange = None
                for tRange in dRanges2Hits.keys():
                    iOverlap = min(iEnd, tRange[1])-max(iStart, tRange[0])+1
                    if iOverlap > iBestOverlap:
                        iBestOverlap = iOverlap
                        tBestOverlapRange = tRange
                if iBestOverlap == 0:
                    tRange = (iStart, iEnd)
                    dRanges2Hits[tRange] = [lHsps]
                    fileRanges.write('%s\t%i\t%i\n'%(sContig, iStart, iEnd))
                    for hsp in lHsps:
                        fileSelectedHits.write('\t'.join(list(map(lambda x:str(x), hsp.show())))+'\n')
                elif len(dRanges2Hits[tBestOverlapRange]) < iMaxHits and lHsps[0].score >= dRanges2Hits[tBestOverlapRange][0][0].score*fBestHitScoreCutoff:
                    dRanges2Hits[tBestOverlapRange].append(lHsps)
                    for hsp in lHsps:
                        fileSelectedHits.write('\t'.join(list(map(lambda x:str(x), hsp.show())))+'\n')
            #break #### for test purpose
        fileRanges.close()
        fileSelectedHits.close()

    def parseSearch(self, sGenomeFileDir, sOrgn, iThreads=8, fEvalCutoff=1e-3, fBestHitScoreCutoff=1/2., iMaxHits=100):
        lProcesses = []
        for iThread in range(iThreads):
            p = Process(target=parseSearchProcess, args=(iThread, iThreads, sGenomeFileDir, sOrgn, fEvalCutoff, fBestHitScoreCutoff, iMaxHits))
            p.start()
            lProcesses.append(p)
        for p in lProcesses:
            p.join()
        fileRanges = open('%s.%s.long.tab.ranges'%(sGenomeFileDir, sDb),'w')
        fileSelectedHits = open('%s.%s.long.tab.selected'%(sGenomeFileDir, sDb),'w')
        for iThread in range(iThreads):
            for sLine in open('%s.%s.long.tab.ranges.thread%i'%(sGenomeFileDir, sDb, iThread)):
                fileRanges.write(sLine)
            for sLine in open('%s.%s.long.tab.selected.thread%i'%(sGenomeFileDir, sDb, iThread)):
                fileSelectedHits.write(sLine)
            os.system('rm -f %s.%s.long.tab.ranges.thread%i'%(sGenomeFileDir, sDb, iThread))
            os.system('rm -f %s.%s.long.tab.selected.thread%i'%(sGenomeFileDir, sDb, iThread))
        fileRanges.close()
        fileSelectedHits.close()

    def assignTaxonomiesToRangesProcess(self, iThread, iThreads, sGenomeFileDir, sOrgn, fEvalCutoff, fIdentityCutoff, fBestHitScoreCutoff, iMaxHits):
        print('thread:%i threads:%i genomeFileDir:%s orgn:%s evalcutoff:%f identitycutoff:%f besthitscorecutoff:%f maxhits:%i'%(iThread, iThreads, sGenomeFileDir, sOrgn, fEvalCutoff, fIdentityCutoff, fBestHitScoreCutoff, iMaxHits))
        if bUseSql:
            connectionThread = sqlite3.connect('%s/%s.sqlite3'%(sDbDir, sDb))
            cursorThread = connectionThread.cursor()
        ## parse taxonomies
        tableOut = open('%s.%s.long.tab.rangeTaxonomies.thread%i'%(sGenomeFileDir, sDb, iThread),'w')
        iHitsIndex = -1
        #iCountRanges = 0 ######
        #iCountHits = 0
        for lHits in ParseBlast(open('%s.%s.long.tab.selected'%(sGenomeFileDir, sDb)), mergeHsps=True):
            iHitsIndex += 1
            if iHitsIndex%iThreads!=iThread:
                continue
            sContig = lHits[0][0].queryId
            dRanges2Hits = {}
            #print(iThread, sContig, len(lHits))
            #if lHits[0][0].eval != min(map(lambda x:x[0].eval, lHits)):
            #    print('asdf')
            #    raise
            #continue
            #if lHits[0][0].eval < fEvalCutoff and lHits[0][0].identity > fIdentityCutoff:
            #    iCountHits += 1
            for lHsps in lHits:
                if lHsps[0].eval > fEvalCutoff:
                    continue
                if lHsps[0].identity < fIdentityCutoff:
                    continue
                ## take only the best hsp
                lHsps = lHsps[:1] ####
                ## skips hits with unresolved taxonomy
                sTaxId = getTaxId( lHsps[0].hitId.split('|')[1], cursor=cursorThread )
                lTaxonomy = taxonomy.getTaxonomyByTaxId(sTaxId)
                lTaxonomy = list(map(lambda x:x[1], lTaxonomy))[::-1]
                if len(lTaxonomy) < 4 or (lTaxonomy[0] != 'viruses' and lTaxonomy[1] not in ('bacteria','archaea','eukaryota')):
                    #print('skipping..', lTaxonomy)
                    continue
                #lTaxonomy[-1] += '.'+sTaxId
                ## create hit
                hit = lHsps[0]
                hit.start = min(list(map(lambda x:x.queryStart, lHsps))+list(map(lambda x:x.queryEnd, lHsps)))
                hit.end = max(list(map(lambda x:x.queryStart, lHsps))+list(map(lambda x:x.queryEnd, lHsps)))
                hit.taxonomy = lTaxonomy
                hit.taxId = sTaxId
                hit.identities = list(map(lambda x:x.identity, lHsps))
                hit.identity = hit.identity/100
                hit.scoreSum = sum(map(lambda x:x.score, lHsps))
                iBestOverlap = 0
                tBestOverlapRange = None
                for tRange in dRanges2Hits.keys():
                    iOverlap = min(hit.end, tRange[1])-max(hit.start, tRange[0])+1
                    if iOverlap > iBestOverlap:
                        iBestOverlap = iOverlap
                        tBestOverlapRange = tRange
                if iBestOverlap == 0:
                    tRange = (hit.start, hit.end)
                    dRanges2Hits[tRange] = [hit]
                ## add hits that overlap well with the range
                elif len(dRanges2Hits[tBestOverlapRange]) < iMaxHits \
                    and hit.score >= dRanges2Hits[tBestOverlapRange][0].score*fBestHitScoreCutoff \
                    and iBestOverlap > (min((tBestOverlapRange[1]-tBestOverlapRange[0]),(hit.end-hit.start))+1)/2: ##
                    dRanges2Hits[tBestOverlapRange].append(hit)
                else:
                    #print(tBestOverlapRange, iBestOverlap, dRanges2Hits[tBestOverlapRange][0].show(), hit.show())
                    pass
            #print('dRanges2Hits', dRanges2Hits.keys())
            #dTaxonomy2Scores = {}
            #iCountRanges += len(dRanges2Hits)
            for (tRange, lHits) in dRanges2Hits.items():
                tableOut.write('%s\t%i-%i'%(sContig, tRange[0], tRange[1]))
                #dTaxonomy2Scores = {}
                fMaxIdentity = max(map(lambda x:x.identity, lHits))
                fIdentityCutoff = fMaxIdentity**2
                #print(tRange, fMaxIdentity, fIdentityCutoff)
                for iHitIndex in range(len(lHits))[::-1]:
                    if hit.identity < fIdentityCutoff:
                        del lHits[iHitIndex]
                fIdentitiesSum = sum(map(lambda x:x.identity, lHits))
                for hit in lHits:
                    tableOut.write('\t%s %s %f %f %f %f'%(hit.hitId, hit.taxId, hit.identity/fIdentitiesSum, hit.identity, hit.eval, hit.score))
                tableOut.write('\n')
        tableOut.close()
    def assignTaxonomiesToRanges(self, sGenomeFileDir, sOrgn, iThreads=8, fIdentityCutoff=0, fEvalCutoff=1e-3, fBestHitScoreCutoff=1/2., iMaxHits=100):
        lProcesses = []
        for iThread in range(iThreads):
            p = Process(target=assignTaxonomiesToRangesProcess, args=(iThread, iThreads, sGenomeFileDir, sOrgn, fEvalCutoff, fIdentityCutoff, fBestHitScoreCutoff, iMaxHits))
            p.start()
            lProcesses.append(p)
        for p in lProcesses:
            p.join()
        tableOut = open('%s.%s.long.tab.rangeTaxonomies'%(sGenomeFileDir, sDb),'w')
        tableOut.write( '#contig\\trange\\ttaxid identityNorm identity eval score\\t...\n' )
        for iThread in range(iThreads):
            for sLine in open('%s.%s.long.tab.rangeTaxonomies.thread%i'%(sGenomeFileDir, sDb, iThread)):
                tableOut.write(sLine)
            os.system('rm -f %s.%s.long.tab.rangeTaxonomies.thread%i'%(sGenomeFileDir, sDb, iThread))
        tableOut.close()

    def appendTaxLevel2Taxon2Scores(self, lTaxLevel2Taxon2ScoresToReturn, lTaxLevel2Taxon2ScoresToAppend, fIdentityNormCutoff=0):
        for iTaxLevel in range(len(lTaxLevel2Taxon2ScoresToAppend)):
            dTaxon2ScoresTemp = {}
            if len(lTaxLevel2Taxon2ScoresToReturn) > iTaxLevel:
                dTaxon2ScoresTemp = lTaxLevel2Taxon2ScoresToReturn[iTaxLevel]
            for (sTaxon, (lIdentityNorms, lMaxIdentities, lRanges, lHitIds)) in lTaxLevel2Taxon2ScoresToAppend[iTaxLevel].items():
                if max(lIdentityNorms) < fIdentityNormCutoff:
                    continue
                if sTaxon not in dTaxon2ScoresTemp:
                    dTaxon2ScoresTemp[sTaxon] = [[],[],[],[]]
                dTaxon2ScoresTemp[sTaxon][0] += lIdentityNorms
                dTaxon2ScoresTemp[sTaxon][1] += lMaxIdentities
                dTaxon2ScoresTemp[sTaxon][2] += lRanges
                dTaxon2ScoresTemp[sTaxon][3] += lHitIds
            if len(dTaxon2ScoresTemp) > 0:
                if len(lTaxLevel2Taxon2ScoresToReturn) <= iTaxLevel:
                    lTaxLevel2Taxon2ScoresToReturn.append({})
                lTaxLevel2Taxon2ScoresToReturn[iTaxLevel] = dTaxon2ScoresTemp
            else:
                break
        return lTaxLevel2Taxon2ScoresToReturn

    def assignTaxonomiesToContigs(self, sGenomeFileDir, sOrgn, fIdentityNormCutoff=1/2):
        ## assign source taxons to contigs
        dContig2TaxLevel2Taxon2Scores = {}
        for sLine in open('%s.%s.long.tab.rangeTaxonomies'%(sGenomeFileDir, sDb)):
            if sLine[0] == '#':
                continue
            lLine = sLine.strip().split('\t')
            if len(lLine) <= 2:
                continue
            (sContig, sRange) = lLine[:2]
            if sContig not in dContig2TaxLevel2Taxon2Scores:
                dContig2TaxLevel2Taxon2Scores[sContig] = []
            tRange = tuple(list(map(lambda x:int(x), sRange.split('-'))))
            lTaxLevel2Taxon2ScoresLoc = []
            for sTaxonScore in lLine[2:]:
                lTaxonScore = sTaxonScore.split()
                sHitId, sTaxId = lTaxonScore[:2]
                sHitId = sHitId.split('|')[1] ##
                (fIdentityNorm, fIdentity, fEval, fScore) = list(map(lambda x:float(x), lTaxonScore[2:]))
                lTaxonomy = taxonomy.getTaxonomyByTaxId(sTaxId)[::-1] ## complete taxonomy
                for iTaxLevel in range(len(lTaxonomy)):
                    if len(lTaxLevel2Taxon2ScoresLoc) <= iTaxLevel:
                        lTaxLevel2Taxon2ScoresLoc.append({})
                    sTaxon = tuple(lTaxonomy[iTaxLevel][:2])
                    #print(sTaxon)
                    if sTaxon not in lTaxLevel2Taxon2ScoresLoc[iTaxLevel]:
                        lTaxLevel2Taxon2ScoresLoc[iTaxLevel][sTaxon] = [[0], [0], [tRange], [sHitId]] ## identityNorm, max identity, range, hit
                    lTaxLevel2Taxon2ScoresLoc[iTaxLevel][sTaxon][0][0] += fIdentityNorm
                    lTaxLevel2Taxon2ScoresLoc[iTaxLevel][sTaxon][1][0] = max(lTaxLevel2Taxon2ScoresLoc[iTaxLevel][sTaxon][1][0], fIdentity)
            dContig2TaxLevel2Taxon2Scores[sContig] = appendTaxLevel2Taxon2Scores(dContig2TaxLevel2Taxon2Scores[sContig], lTaxLevel2Taxon2ScoresLoc, 0)
        tableOut = open('%s.%s.long.tab.contigTaxonomies'%(sGenomeFileDir, sDb),'w')
        tableOut.write('#contig\tbest\tranges\thitIds\ttaxLevel\ttaxId\ttaxon\ttaxonomy\tnRangesOverall\tnRanges\tsumIdentityNorms\tmaxIdentityNorms\tidentityNorms\tmaxIdentity\tidentities\n')
        for (sContig, lTaxLevel2Taxon2Scores) in sorted(dContig2TaxLevel2Taxon2Scores.items()):
            iRanges = sum( map(lambda x:len(x[1][3]), lTaxLevel2Taxon2Scores[0].items()) )
            #print(iRanges, lTaxLevel2Taxon2Scores[0])
            bBreak = False
            for iTaxLevel in range(len(lTaxLevel2Taxon2Scores))[::-1]: ##
                for ((sTaxId, sTaxon), (lIdentityNorms, lMaxIdentities, lRanges, lHitIds)) in sorted(lTaxLevel2Taxon2Scores[iTaxLevel].items(), key=lambda x:(-sum(x[1][0]), -sum(x[1][1]))):
                    iCountIdentityNormsAbove50 = 0 ##
                    for fIdentityNorm in lIdentityNorms:
                        if fIdentityNorm >= 0.5:
                            iCountIdentityNormsAbove50 += 1
                    iBest = 0
                    if iCountIdentityNormsAbove50 >= iRanges/2: ##
                        if not bBreak:
                            iBest = 1
                    if iBest:
                        lTaxonomy = list(map(lambda x:x[1], taxonomy.getTaxonomyByTaxId(sTaxId)))[::-1]
                        sRanges = ' '.join(map(lambda x:str(x[0])+'-'+str(x[1]), lRanges))
                        tableOut.write(sContig+'\t'+\
                            str(iBest)+'\t'+\
                            sRanges+'\t'+\
                            ' '.join(lHitIds)+'\t'+\
                            str(iTaxLevel)+'\t'+\
                            sTaxId+'\t'+\
                            sTaxon+'\t'+\
                            '; '.join(lTaxonomy)+'\t'+\
                            str(iRanges)+'\t'+\
                            str(len(lIdentityNorms))+'\t'+\
                            str(sum(lIdentityNorms))+'\t'+\
                            str(max(lIdentityNorms))+'\t'+\
                            ' '.join(list(map(lambda x:str(x), lIdentityNorms)))+'\t'+\
                            str(max(lMaxIdentities))+'\t'+\
                            ' '.join(list(map(lambda x:str(x), lMaxIdentities)))+'\n')
                        bBreak = True
                        break
                if bBreak:
                    break
        tableOut.close()

    def summarizeTaxonomies(self, sGenomeFileDir, sOrgn, fIdentityNormCutoff=1/2, fEvalCutoff=1e-5, fIdentityCutoff=0):
        lTaxLevel2Taxon2Scores = []
        dTaxon2Scores = {}
        for sLine in open('%s.%s.long.tab.rangeTaxonomies'%(sGenomeFileDir, sDb)):
            if sLine[0] == '#':
                continue
            lLine = sLine.strip().split('\t')
            if len(lLine) <= 2:
                continue
            (sContig, sRange) = lLine[:2]
            tRange = tuple(list(map(lambda x:int(x), sRange.split('-'))))
            for sTaxonScore in lLine[2:]:
                lTaxonScore = sTaxonScore.split()
                sHitId, sTaxId = lTaxonScore[:2]
                sHitId = sHitId.split('|')[1] ##
                (fIdentityNorm, fIdentity, fEval, fScore) = list(map(lambda x:float(x), lTaxonScore[2:]))
                if fIdentityNorm < fIdentityNormCutoff:
                    continue
                if fIdentity < fIdentityCutoff:
                    continue
                if fEval > fEvalCutoff:
                    continue
                lTaxonomy = taxonomy.getTaxonomyByTaxId(sTaxId)[::-1] ## complete taxonomy
                for iTaxLevel in range(len(lTaxonomy)):
                    if len(lTaxLevel2Taxon2Scores) <= iTaxLevel:
                        lTaxLevel2Taxon2Scores.append({})
                    sTaxon = tuple(list(map(lambda x:x[1], lTaxonomy[:iTaxLevel+1])))
                    if sTaxon not in lTaxLevel2Taxon2Scores[iTaxLevel]:
                        lTaxLevel2Taxon2Scores[iTaxLevel][sTaxon] = [0, 0]
                    lTaxLevel2Taxon2Scores[iTaxLevel][sTaxon][0] += fIdentityNorm
                    lTaxLevel2Taxon2Scores[iTaxLevel][sTaxon][1] += fIdentity
                ## for krona tools
                tTaxonomy = tuple(list(map(lambda x:x[1], lTaxonomy)))
                if tTaxonomy not in dTaxon2Scores:
                    dTaxon2Scores[tTaxonomy] = []
                dTaxon2Scores[tTaxonomy].append(fIdentityNorm)
        tableOut = open('%s.%s.long.tab.taxons'%(sGenomeFileDir, sDb),'w')
        tableOut.write('#taxLevel\ttaxonomy\ttaxon\tidentityNormsSum\tidentitiesSum\n')
        for iTaxLevel in range(len(lTaxLevel2Taxon2Scores)):
            for (sTaxon, (fIdentityNormSum, fIdentitySum)) in sorted(lTaxLevel2Taxon2Scores[iTaxLevel].items(), key=lambda x:(-x[1][0], -x[1][1])):
                #print(sTaxon)
                #print(iTaxLevel, ';'.join(sTaxon), sTaxon[iTaxLevel], fIdentityNormSum, fIdentitySum)
                tableOut.write('%i\t%s\t%s\t%f\t%f\n'%(iTaxLevel, ';'.join(sTaxon), sTaxon[iTaxLevel], fIdentityNormSum, fIdentitySum))
                #tableOut.write('\t'+'\t'.join(map(lambda x:str(x), sorted(lIdentities, reverse=True))))
                #tableOut.write('\n')
        tableOut.close()
        tableKronaOut = open('%s.%s.long.tab.taxons.krona'%(sGenomeFileDir, sDb),'w')
        for (tTaxonomy, lScores) in sorted(dTaxon2Scores.items(), key=lambda x:-sum(x[1])):
            fScore = sum(lScores)#/len(lScores)
            tableKronaOut.write('%f\t%s\n'%(fScore, '\t'.join(tTaxonomy)))
        tableKronaOut.close()

## read assemblies
lSeqProjects = []
for sLine in open('parabasalids_data.tsv'):
    lLine = sLine.strip('\n').split('\t')
    (sUse, sName, sReads, sAssembly) = lLine
    if sUse != '1':
        continue
    lSeqProjects.append((sName, sAssembly, sReads.split()))
#lAssemblies = [('test','test.fna')]

#sDB = 'nr.200929'
#sDB = 'swissprot.200929'
sDb = 'uniprot_trembl.201002'
sDbDir = '/Data/vojta/databases/uniprot/'
bUseSql = True


#dSilvaId2Taxonomy = {}
#print('reading SILVA taxonomy...')
#for seq in SeqIO(open('/home/zarsky/databases/silva/SILVA_138.1_SSURef_tax_silva.fasta')):
#    dSilvaId2Taxonomy[seq.id] = seq.desc
#print('done.')

#taxonomy = Taxonomy(sTaxdumpDir='/Data/vojta/databases/ncbi/taxdump.201002')

for (sOrgn, sGenomeFileDir, lReadFileDirs) in lSeqProjects:
    seqproj = SeqProject(sOrgn, sGenomeFileDir, [[lReadFileDirs[0]],[lReadFileDirs[1]]], sDb, sDbDir, bUseSql)
    seqproj.mapReads()
    #seqproj.getSsurnaModels('/Data/vojta/databases/rfam/Rfam.14.cm')
    #genome = Genome(sGenomeFileDir)
    #iNSeqs, iSumLen, iN50, lLengths, lGCFractions, fGC, iNs = genome.getStats(bGC=True)
    #lPrint = ['Orgn: %s, Nseqs: %i, Sum: %i, N50: %i, Ns: %f, GC: %f\t'%(sOrgn, iNSeqs, iSumLen, iN50, float(iNs)/iSumLen, fGC)]
    #for i in sorted(range(len(lLengths)), key=lambda x:-lLengths[x]):
    #    lPrint.append('\t'.join(map(lambda x:str(x), [lLengths[i], lGCFractions[i]])))
    #lPrints.append(lPrint)

    ## map reads
    #mapReads(sGenomeFileDir, lReadFileDirs)

    #os.system('busco -c 12 -i %s -l eukaryota_odb10 -o busco.%s -m genome'%(sGenomeFileDir, sOrgn))
    #os.system('/Data/vojta/software/diamond blastx -b8 -c1 --threads 12 --db /Data/vojta/databases/ncbi/nr.200929 --query %s --out %s.nr_200929_fast.daa --outfmt 100'%(sGenomeFileDir, sGenomeFileDir)) ## for quick search

    ## Diamond
    #os.system('/home/zarsky/software/diamond blastx -b8 -c1 --threads 24 --db /home/zarsky/databases/uniprot/%s --query %s --out %s.%s.fast.tab --outfmt 6 --evalue 1e-3 --max-target-seqs 10'%(sDb, sGenomeFileDir, sGenomeFileDir, sDb)) ## fast
    #os.system('/home/zarsky/software/diamond blastx -b8 -c1 --threads 24 --db /home/zarsky/databases/uniprot/%s --query %s --out %s.%s.long.tab --outfmt 6 --evalue 1e-1 --max-target-seqs 1000000 --max-hsps 1000'%(sDb, sGenomeFileDir, sGenomeFileDir, sDb)) ## for long output
    #os.system('pigz -p 8 --best %s.%s.long.tab'%(sGenomeFileDir, sDb))

    #os.system('/home/zarsky/software/diamond blastx -b8 -c1 --threads 24 --db /home/zarsky/databases/uniprot/%s --query %s --out %s.%s.sensitive.tab --outfmt 6 --evalue 1e-1 --max-target-seqs 1000000 --max-hsps 1000 --sensitive'%(sDb, sGenomeFileDir, sGenomeFileDir, sDb)) ## for deeper research

    ## BlobTools
    #blobTools(sGenomeFileDir, sOrgn)

    ## Find ORFs
    #os.system('/home/zarsky/software/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t %s'%(sGenomeFileDir))

    #lTargetTaxons = [('metamonada', 2),]
    ## Parse DIAMOND search
    #parseSearch(sGenomeFileDir, sOrgn, iThreads=24)
    ## Assign taxonomies to contig regions
    #dContigs = {}
    #for lHits in ParseBlast(open('%s.%s.long.tab.selected'%(sGenomeFileDir, sDb)), mergeHsps=True):
    #    if lHits[0][0].eval > 1e-5:
    #        continue
    #    dContigs[lHits[0][0].queryId] = lHits[0][0].eval
    #print(len(dContigs))
    #assignTaxonomiesToRanges(sGenomeFileDir, sOrgn, iThreads=24, fIdentityCutoff=0, fEvalCutoff=1e-5, fBestHitScoreCutoff=1/2.) ##
    #assignTaxonomiesToContigs(sGenomeFileDir, sOrgn, fIdentityNormCutoff=1/2)
    #summarizeTaxonomies(sGenomeFileDir, sOrgn, fIdentityNormCutoff=1/2, fEvalCutoff=1e-10, fIdentityCutoff=0)
    ## BUSCO on target taxons
    #for (sTargetTaxon, iTargetTaxonLevel) in lTargetTaxons:
    #    os.system('busco --config /home/zarsky/software/busco_v4/config/config.ini -c 12 -i %s/target_taxons/%s.fasta -l eukaryota_odb10 -o busco.%s.%s -m genome'%(sGenomeFileDir.rsplit('/',1)[0], sTargetTaxon, sOrgn, sTargetTaxon))
    ## get assemby stats
    #for (sTargetTaxon, iTargetTaxonLevel) in lTargetTaxons:
    #    print(sOrgn, sTargetTaxon)
    #    genome = Genome('%s/target_taxons/%s.fasta'%(sGenomeFileDir.rsplit('/',1)[0], sTargetTaxon))
    #    iNSeqs, iSumLen, iN50, lLengths, lGCFractions, fGC, iNs = genome.getStats(bGC=True)
    #    print('Orgn: %s, Nseqs: %i, Sum: %i, N50: %i, Ns: %f, GC: %f\t'%(sOrgn, iNSeqs, iSumLen, iN50, float(iNs)/iSumLen, fGC))
    ##
    #sortReads(sGenomeFileDir, sOrgn, lTargetTaxons)    
    ## rRNA
    #scanForRRNA(sGenomeFileDir, sOrgn)
    #analyzeRRNA('%s.rrnas.fasta'%sOrgn, dSilvaId2Taxonomy)
    #break ## for test purpose
    #pass

"""
sKrontaToolsCommand = 'ktImportText -o parabasalids.krona0.html'
for (sOrgn, sGenomeFileDir, lReadFileDirs) in lAssemblies:
    sKrontaToolsCommand += ' %s.%s.long.tab.taxons.krona,%s'%(sGenomeFileDir, sDb, sOrgn)
print(sKrontaToolsCommand)
os.system(sKrontaToolsCommand)
"""
#summarizeTaxonomiesKrona(lAssemblies)

