import os
import sys
import re
sys.path.append('/Data/vojta/myclasses/')
from SeqIO import SeqIO
from ParseBlast import ParseBlast

## blastn -task blastn -query final.contigs.fa.ssu.fna -db /Data/databases/silva/SILVA_138.1_SSURef_NR99_tax_silva.fasta -out final.contigs.fa.ssu.fna.silva.blast -num_threads 8 -max_target_seqs 1000 -outfmt 6


try:
    sQueryFile = sys.argv[1]
    if not os.path.isfile(sQueryFile):
        raise
    sBlastFile = sys.argv[2]
    if not os.path.isfile(sBlastFile):
        raise
    sDbFile = sys.argv[3]
    if not os.path.isfile(sDbFile):
        raise
    iBestHits = int(sys.argv[4])
    #sOutFile = sys.argv[3]
    #if not os.path.isfile(sOutFile):
    #    raise
except:
    print('Usage:\npython parseSilva.py <queryfile.fna> <blastfile.blast> <dbfile> <nbesthits>')

fEvalCutoff = 1e-10
#iLenCutoff = 100
#bMegahit = False
#if 'megahit' in sQueryFile:
#    bMegahit = True

tableOut = open('{}.parseSilva.tsv'.format(sQueryFile.rsplit('.',1)[0]),'w')
tableOut.write('id\tseq\tlen\tcov\tbestHit\tscore\teval\tidentity\tspecies\ttaxonomy\n')
if not os.path.isdir('parseSilva_dir'):
    os.system('mkdir parseSilva_dir')

## SILVA !!!!!
dDbSeqs = SeqIO(sDbFile).getDict()

dQuerySeqs = SeqIO(sQueryFile).getDict()
for lHits in ParseBlast([sBlastFile], mergeHsps=True):
    print(lHits[0][0].queryId, len(lHits))
    sQueryId = lHits[0][0].queryId
    seqQuery = dQuerySeqs[sQueryId]
    fileOut = open('./parseSilva_dir/{}.fasta'.format(sQueryId),'w')
    fileOut.write('>{} {}\n{}\n'.format(seqQuery.id, seqQuery.desc, seqQuery.seq))
    lBestHit = ['',-1,0,'','']
    for lHsps in lHits[:iBestHits]:
        hsp = lHsps[0]
        if hsp.eval <= fEvalCutoff:
            sHitId = hsp.hitId.strip()
            ## large DB !!!!!!
            #sRead = os.popen('blastdbcmd -db {} -entry {} -outfmt %f'.format(sDbFile, sHitId)).read()
            #sId = sRead.split()[0][1:]
            #sTaxonomy = sRead.split(' ',1)[1].split('\n')[0]
            #lTaxonomy = sTaxonomy.split(';')
            #sSeq = sRead.split('\n',1)[1].replace('\n','')
            ## SILVA !!!!!!
            seqDb = dDbSeqs[sHitId]
            sTaxonomy = seqDb.desc.strip()
            lTaxonomy = sTaxonomy.split(';')
            sSeq = seqDb.seq
            #print(lTaxonomy)
            if len(lTaxonomy) == 2:
                lTaxonomy = [lTaxonomy[0],'None',lTaxonomy[1]]
            fileOut.write('>{}.{}.{}.{} {}\n{}\n'.format(sHitId, re.sub('\W','',lTaxonomy[2]).replace(' ','_'), re.sub('\W','','_'.join(lTaxonomy[-1].split()[:2])), int(hsp.identity), sTaxonomy, sSeq))
            if lBestHit[0] == '':
                lBestHit = [sHitId, hsp.score, hsp.eval, hsp.identity, lTaxonomy[-1]]+lTaxonomy
    print(lBestHit)
    fileOut.close()
    fCoverage = ''
    if '_cov_' in seqQuery.id: ## spades
        lTemp = seqQuery.id.split('_cov_')[1].split('_')[0].split('.')
        fCoverage = float(lTemp[0])
        if lTemp[1].isdigit():
            fCoverage = float('.'.join(lTemp[:2]))
    elif ' cov_' in seqQuery.desc: ## trinity
        fCoverage = float(seqQuery.desc.split(' cov_')[1].split()[0])
    elif 'multi=' in seqQuery.desc: ## megahit
        fCoverage = float(seqQuery.desc.split('multi=')[1].split()[0])
    tableOut.write('\t'.join(list(map(lambda x:str(x), [sQueryId, seqQuery.seq, len(seqQuery.seq), fCoverage ]+lBestHit)))+'\n')
    tableOut.flush()
    #break
tableOut.close()


