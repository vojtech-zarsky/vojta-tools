import os
import re
import pickle
from multiprocessing import Process
import numpy
#import sklearn
from sklearn.svm import SVC

from SeqIO import SeqIO
from ParseHmmer import ParseHmmer
#from Taxonomy import Taxonomy
#from ParseBlast import ParseBlast


class Annotation:
    def __init__(self, infilename=None, dir='./', moltype='aa', cpus=1):
        self.sInfilename = infilename
        if self.sInfilename != None:
            self.sDir = dir+'/'
            self.sMoltype = moltype
            self.iCpus = cpus
            self.dSeqId2Annotation = {}
            self.lColumns = ['id']
            ## partition by threads
            if self.iCpus > 0:
                lThread2File = []
                for iThread in range(self.iCpus):
                    lThread2File.append( open(self.sDir+self.sInfilename+'.thread'+str(iThread), 'w') )
                iThread = 0
                for seq in SeqIO(open(self.sDir+self.sInfilename)):
                    lThread2File[iThread].write('>'+seq.id+' '+seq.desc+'\n'+seq.seq+'\n')
                    iThread += 1
                    if iThread >= self.iCpus:
                        iThread = 0
                for file in lThread2File:
                    file.close()
    def parseSeqProperties(self, parseIdFunction=None, parseDescFunction=None):
        for seq in SeqIO(open(self.sDir+self.sInfilename)):
            if seq.id not in self.dSeqId2Annotation:
                self.dSeqId2Annotation[seq.id] = {'id':seq.id}
            if parseIdFunction != None:
                self.dSeqId2Annotation[seq.id]['id'] = parseIdFunction(self.dSeqId2Annotation[seq.id]['id'])
            self.dSeqId2Annotation[seq.id]['desc'] = seq.desc
            if parseDescFunction != None:
                self.dSeqId2Annotation[seq.id]['desc'] = parseDescFunction(self.dSeqId2Annotation[seq.id]['desc'])
            self.dSeqId2Annotation[seq.id]['len'] = len(seq.seq)
            ## transdecoder/selenidium specific !!!
            #self.dSeqId2Annotation[seq.id]['id'] = self.dSeqId2Annotation[seq.id]['id'].split(':')[0]
        self.lColumns += ['desc', 'len']
    def runEggnog(self, emapperdir='/home/zarsky/eggnog-mapper', eggnogdb='euk'):
        self.sEggnogdb = eggnogdb
        os.system('python '+emapperdir+'/emapper.py -i '+self.sDir+self.sInfilename+' --output '+self.sDir+self.sInfilename+'.'+self.sEggnogdb+' -d '+self.sEggnogdb+' --usemem --cpu '+str(self.iCpus))
    def parseEggnog(self, eggnogfilename=None, eggnogdb='euk'):
        self.sEggnogdb = eggnogdb
        if eggnogfilename == None:
            self.sEggnogfilename = self.sInfilename+'.'+self.sEggnogdb+'.emapper.annotations'
        else:
            self.sEggnogfilename = eggnogfilename
        #lEggnogColumns = None
        for sLine in open(self.sDir+self.sEggnogfilename):
            lLine = sLine.strip('\n').split('\t')
            if lLine[0] == '#query_name':
                lEggnogColumns = sLine[1:].strip('\n').split('\t')
            if sLine[0] == '#':
                continue
            sSeqId = lLine[0]
            sGoTerms = lLine[5]
            (sOg, sEval, sScore) = lLine[10].split('|')
            lCogs = lLine[11].split(', ')
            sDesc = lLine[12]
            if sSeqId not in self.dSeqId2Annotation:
                self.dSeqId2Annotation[sSeqId] = {'id':sSeqId}
            self.dSeqId2Annotation[sSeqId]['gos'] = sGoTerms
            self.dSeqId2Annotation[sSeqId]['eggnogOg'] = sOg
            self.dSeqId2Annotation[sSeqId]['eggnogEval'] = sEval
            self.dSeqId2Annotation[sSeqId]['cogs'] = ' '.join(lCogs)
            self.dSeqId2Annotation[sSeqId]['eggnogDesc'] = sDesc
            self.dSeqId2Annotation[sSeqId]['keggOg'] = lLine[6]
        self.lColumns += ['eggnogOg', 'eggnogDesc', 'eggnogEval', 'cogs', 'gos']
    def readKegg(self, sDir = '/home/zarsky/databases/kegg/'):
        ## read kegg
        dKegg = {}
        lKeggCategories = ['pathwayFamilies','pathways','modules','reactions','orthologs','compounds']
        for sKeggCategory in lKeggCategories:
            dKegg[sKeggCategory] = {}
            lColumns = None
            dColumns2Process = {}
            iLine = 0
            for sLine in open(sDir+'kegg.'+sKeggCategory+'.tsv'):
                lLine = sLine.strip('\n').split('\t')
                sPrimaryId = lLine[0]
                if lColumns == None:
                    lColumns = map(lambda x:x.rsplit('.',2), lLine)
                    continue
                dKegg[sKeggCategory][sPrimaryId] = {'line':iLine}
                iLine += 1
                for iColumn in range(len(lColumns)):
                    (sColumn, sProcess, sType) = lColumns[iColumn]
                    if sProcess == 's':
                        dKegg[sKeggCategory][sPrimaryId][sColumn] = lLine[iColumn]
                    elif sProcess == 'l':
                        dKegg[sKeggCategory][sPrimaryId][sColumn] = lLine[iColumn].split(';')
                        if dKegg[sKeggCategory][sPrimaryId][sColumn] == ['']:
                            dKegg[sKeggCategory][sPrimaryId][sColumn] = []
                    elif sProcess == 'll':
                        dKegg[sKeggCategory][sPrimaryId][sColumn] =  map(lambda x:x.split(), lLine[iColumn].split(';'))
                        if dKegg[sKeggCategory][sPrimaryId][sColumn] == [[]]:
                            dKegg[sKeggCategory][sPrimaryId][sColumn] = []
                    elif sProcess == 'lll':
                        dKegg[sKeggCategory][sPrimaryId][sColumn] = map(lambda x:map(lambda y:y.split(), x.split(',')), lLine[iColumn].split(';'))
                    else:
                        print(sKeggCategory)
                        print(lLine)
                        print(sColumn, sProcess, sType)
                        raise
        for (sOrtholog, dOrtholog) in dKegg['orthologs'].items():
            dKegg['orthologs'][sOrtholog]['pathways'] = []
        for (sPathway, dPathway) in sorted(dKegg['pathways'].items()):
            for sOrtholog in dPathway['orthologs']:
                dKegg['orthologs'][sOrtholog]['pathways'].append(sPathway)
        self.dKegg = dKegg
        return dKegg
    def runKofamProcess(self, iThread, kofamdb):
        os.system('hmmscan --cpu 1 --domtblout '+self.sDir+self.sInfilename+'.thread'+str(iThread)+'.kofam.domtblout '+kofamdb+' '+self.sDir+self.sInfilename+'.thread'+str(iThread)+' > dump')
    def runKofam(self, kofamdb='/home/zarsky/databases/kegg/kofam/kofam.hmm'):
        lProcesses = []
        for iThread in range(self.iCpus):
            p = Process(target=self.runKofamProcess, args=(iThread,kofamdb))
            p.start()
            lProcesses.append(p)
        for p in lProcesses:
            p.join()
        fileOut = open(self.sDir+self.sInfilename+'.kofam.domtblout','w')
        for iThread in range(self.iCpus):
            fileOut.write(open(self.sDir+self.sInfilename+'.thread'+str(iThread)+'.kofam.domtblout').read())
            os.system('rm '+self.sDir+self.sInfilename+'.thread'+str(iThread)+'.kofam.domtblout')
        fileOut.close()
    def parseKofamFromEggnog(self):
        #lKeggCategories = ['pathwayFamilies','pathways','modules','reactions','orthologs','compounds']
        for (sSeqId, dAnnotation) in self.dSeqId2Annotation.items():
            dPathways = {}
            dPathwayFamilies = {}
            dPathwaySuperfamilies = {}
            sKeggName = ''
            if 'keggOg' in dAnnotation and dAnnotation['keggOg'] != '':
                lKeggOgs = dAnnotation['keggOg'].strip().split(',')
                while len(lKeggOgs) > 0 and lKeggOgs[0] not in self.dKegg['orthologs']:
                    del lKeggOgs[0]
                sKeggOg = None
                if len(lKeggOgs) > 0:
                    sKeggOg = lKeggOgs[0]
                if sKeggOg != None:
                    sKeggName = self.dKegg['orthologs'][sKeggOg]['name']
                    for dPathway in self.dKegg['pathways'].values():
                        if sKeggOg in dPathway['orthologs']:
                            dPathways[dPathway['name']] = None
                            dPathwayFamilies[dPathway['pathwayFamily']] = None
                            dPathwaySuperfamilies[dPathway['pathwaySuperfamily']] = None
            self.dSeqId2Annotation[sSeqId]['kegg.name'] = sKeggName
            self.dSeqId2Annotation[sSeqId]['kegg.pathways'] = sorted( dPathways.keys() )
            self.dSeqId2Annotation[sSeqId]['kegg.pathwayFamilies'] = sorted( dPathwayFamilies.keys() )
            self.dSeqId2Annotation[sSeqId]['kegg.pathwaySuperfamilies'] = sorted( dPathwaySuperfamilies.keys() )
        self.lColumns += ['kegg.name', 'kegg.pathways', 'kegg.pathwayFamilies', 'kegg.pathwaySuperfamilies']
    def parseKofam(self, kofamfilename=None, evalcutoff=1e-10):
        #lKeggCategories = ['pathwayFamilies','pathways','modules','reactions','orthologs','compounds']
        if kofamfilename == None:
            kofamfilename = self.sInfilename+'.kofam.tsv'
        for lHits in ParseHmmer(open(self.sDir+kofamfilename), nooverlap=True):
            if len(lHits) == 0:
                continue
            sSeqId = lHits[0].queryId
            self.dSeqId2Annotation[sSeqId]['kegg.ogs'] = []
            self.dSeqId2Annotation[sSeqId]['kegg.names'] = []
            self.dSeqId2Annotation[sSeqId]['kegg.pathways'] = []
            self.dSeqId2Annotation[sSeqId]['kegg.pathwayFamilies'] = []
            self.dSeqId2Annotation[sSeqId]['kegg.pathwaySuperfamilies'] = []
            for hit in lHits:
                if hit.eval > evalcutoff:
                    continue
                sKeggOg = hit.hitId
                self.dSeqId2Annotation[sSeqId]['kegg.ogs'].append(sKeggOg)
                self.dSeqId2Annotation[sSeqId]['kegg.names'].append( self.dKegg['orthologs'][sKeggOg]['name'] )
                dPathways = {}
                dPathwayFamilies = {}
                dPathwaySuperfamilies = {}
                for dPathway in self.dKegg['pathways'].values():
                    if sKeggOg in dPathway['orthologs']:
                        dPathways[dPathway['name']] = None
                        dPathwayFamilies[dPathway['pathwayFamily']] = None
                        dPathwaySuperfamilies[dPathway['pathwaySuperfamily']] = None
                self.dSeqId2Annotation[sSeqId]['kegg.pathways'].append( sorted( dPathways.keys() ) )
                self.dSeqId2Annotation[sSeqId]['kegg.pathwayFamilies'].append( sorted( dPathwayFamilies.keys() ) )
                self.dSeqId2Annotation[sSeqId]['kegg.pathwaySuperfamilies'].append( sorted( dPathwaySuperfamilies.keys() ) )
        self.lColumns += ['kegg.ogs', 'kegg.names', 'kegg.pathways', 'kegg.pathwayFamilies', 'kegg.pathwaySuperfamilies']
    def runPfam(self, pfamdb='/home/zarsky/databases/pfam/Pfam-A.32.hmm'):
        self.pfamDbName = pfamdb.rsplit('/',1)[-1].rsplit('.',1)[0]
        os.system('hmmscan --cpu '+str(self.iCpus)+' --domtblout '+self.sDir+self.sInfilename+'.pfam.domtblout '+pfamdb+' '+self.sDir+self.sInfilename)
    def parsePfam(self, pfamfilename, evalcutoff=1e-5):
        for lHits in ParseHmmer(open(self.sDir+pfamfilename), nooverlap=True):
            if len(lHits) == 0:
                break
            lPfams = []
            for hit in sorted(lHits, key=lambda x:x.queryStart):
                if hit.eval <= 1e-5:
                    lPfams.append(hit.hitId)
            if len(lPfams) > 0:
                sSeqId = lHits[0].queryId
                if sSeqId not in self.dSeqId2Annotation:
                    self.dSeqId2Annotation[sSeqId] = {'id':sSeqId}
                self.dSeqId2Annotation[sSeqId]['pfam'] = ' '.join(lPfams)
        self.lColumns += ['pfam']
    def predictPtsRes(self):
        dPtsRe = {'1':{}, '2':{}}
        dPtsRe['1']['general'] ='(S|A|C)(K|R|H|Q)(L|M)'
        dPtsRe['1']['kinetoplastids'] ='(S|A|G|C|N|P)(R|H|K|N|Q)(L|I|V|F|A|M|Y)'
        dPtsRe['2']['general'] ='R(L|I|V|Q)\w{2}(L|I|V|Q|H)(L|S|G|A)\w{1}(H|Q)(L|A)'
        dPtsRe['2']['kinetoplastids'] ='(R|K)(L|V|I)\w{5}(H|K|Q|R)(L|A|I|V|F|Y)'
        for seq in SeqIO(open(self.sDir+self.sInfilename)):
            if seq.id not in self.dSeqId2Annotation:
                self.dSeqId2Annotation[seq.id] = {'id':seq.id}
            for sPtsClass in ['1', '2']:
                sSeq = seq.seq
                if sPtsClass == '1':
                    sSeq = sSeq[-3:]
                elif sPtsClass == '2':
                    sSeq = sSeq[:30]
                for sPtsModel in ['general', 'kinetoplastids']:
                    sMatch = ''
                    match = re.search(dPtsRe[sPtsClass][sPtsModel], sSeq)
                    if match != None:
                        sMatch = match.string[match.start():match.end()]
                    self.dSeqId2Annotation[seq.id]['pts'+sPtsClass+'.'+sPtsModel] = sMatch
        for sPtsClass in ['1', '2']:
            for sPtsModel in ['general', 'kinetoplastids']:
                self.lColumns += ['pts'+sPtsClass+'.'+sPtsModel]
    def predictPts1(self):
        lPosition2Residue2Score = [{'S':1., 'A':1., 'C':1., 'H':.5, 'K':.5, 'N':.5, 'P':.5, 'T':.5}, {'K':1., 'R':1., 'H':1., 'N':.5, 'Q':.5, 'S':.5}, {'L':1., 'I':.5, 'M':.5, 'F':.5, 'A':.5, 'V':.5, 'Y':.5}]
        for seq in SeqIO(open(self.sDir+self.sInfilename)):
            if seq.id not in self.dSeqId2Annotation:
                self.dSeqId2Annotation[seq.id] = {'id':seq.id}
            sPts1 = seq.seq.upper()[-3:]
            if len(sPts1) < 3:
                self.dSeqId2Annotation[seq.id]['pts1seq'] = ''
                self.dSeqId2Annotation[seq.id]['pts1matches'] = 0
                self.dSeqId2Annotation[seq.id]['pts1score'] = 0
                continue
            fScore = 0
            iMatches = 0
            for iPos in range(3):
                if sPts1[iPos] in lPosition2Residue2Score[iPos]:
                    fScore += lPosition2Residue2Score[iPos][sPts1[iPos]]
                    iMatches += 1
            self.dSeqId2Annotation[seq.id]['pts1seq'] = sPts1
            self.dSeqId2Annotation[seq.id]['pts1matches'] = iMatches
            self.dSeqId2Annotation[seq.id]['pts1score'] = fScore
        self.lColumns += ['pts1score', 'pts1matches', 'pts1seq']
    def seq2array(self, sSeq):
        sAas = 'ARNDCQEGHILKMFPSTWYV'
        aArray = numpy.zeros((len(sSeq)*20,))
        for i in range(len(sSeq)):
            try:
                aArray[i*20+sAas.index(sSeq[i])] = 1
            except:
                pass
        return aArray
    def trainPts1ML(self, sPositiveFile, sNegativeFile, iCLen=10):
        sAas = 'ARNDCQEGHILKMFPSTWYV'
        self.iCLen = 10
        ## read positive set
        dPts1Used = {}
        lPts1 = []
        for seq in SeqIO(open(sPositiveFile)):
            sPts1 = seq.seq.strip('*').upper().replace('U','S')[-iCLen:]
            #sPts1 = seq.seq.strip('*').upper().replace('U','S')[-iCLen-1:-1] ## for test purpose !!!
            if sPts1 in dPts1Used:
                continue
            dPts1Used[sPts1] = None
            if len(sPts1) >= iCLen and not re.search('[^'+sAas+']', sPts1):
                lPts1.append( sPts1 )
        ## read negative set
        dPts1Used = {}
        lNonPts1 = []
        lNegativeFiles = [(sNegativeFile, 0)]
        for (sNegativeFile, bPeroxisomes) in lNegativeFiles:
            for seq in SeqIO(open(sNegativeFile)):
                if bPeroxisomes:
                    seq.seq = seq.seq.strip('*')[:-1]
                sPts1 = seq.seq.strip('*').upper().replace('U','S')[-iCLen:]
                if sPts1 in dPts1Used:
                    continue
                dPts1Used[sPts1] = None
                if len(sPts1) >= iCLen and not re.search('[^'+sAas+']', sPts1):
                    lNonPts1.append( sPts1 )
        ## sequences to arrays
        x_pts = []
        y_pts = []
        for i in range(len(lPts1)):
            sPts1 = lPts1[i]
            aX = self.seq2array(sPts1)
            x_pts.append(aX)
            y_pts.append(1)
        x_pts = numpy.array(x_pts)
        y_pts = numpy.array(y_pts)
        print('x_pts', x_pts.shape)
        x_nonpts = []
        y_nonpts = []
        for i in range(len(lNonPts1)):
            sPts1 = lNonPts1[i]
            aX = self.seq2array(sPts1)
            x_nonpts.append(aX)
            y_nonpts.append(0)
        x_nonpts = numpy.array(x_nonpts)
        y_nonpts = numpy.array(y_nonpts)
        x = numpy.concatenate((x_pts, x_nonpts))
        y = numpy.concatenate((y_pts, y_nonpts))
        print('x_nonpts', x_nonpts.shape)
        #model = sklearn.svm.SVC(probability=True)
        model = SVC(probability=True)
        model.fit(x, y)
        pickle.dump(model, open('trainPts1ML.model.pickle','wb'))
    def predictPts1ML(self, sModelFile='trainPts1ML.model.pickle', iCLen=10):
        model = pickle.load(open(sModelFile,'rb'))
        for seq in SeqIO(open(self.sDir+self.sInfilename)):
            if seq.id not in self.dSeqId2Annotation:
                self.dSeqId2Annotation[seq.id] = {'id':seq.id}
            sPts1 = seq.seq.strip('*').upper().replace('U','S')[-iCLen:]
            if len(sPts1) < iCLen:
                self.dSeqId2Annotation[seq.id]['pts1MLseq'] = ''
                self.dSeqId2Annotation[seq.id]['pts1MLprob'] = 0
                continue
            aSeq = [self.seq2array(sPts1)]
            fProb = model.predict_proba(aSeq)[0][1]
            self.dSeqId2Annotation[seq.id]['pts1MLseq'] = sPts1
            self.dSeqId2Annotation[seq.id]['pts1MLprob'] = fProb
        self.lColumns += ['pts1MLprob', 'pts1MLseq']
    def predictPts2(self):
        for seq in SeqIO(open(self.sDir+self.sInfilename)):
            if seq.id not in self.dSeqId2Annotation:
                self.dSeqId2Annotation[seq.id] = {'id':seq.id}
            sSeq = seq.seq.strip('*').upper()
            iMotifStart = None
            fPts2Score = 0
            if re.search('[RK][LVIQ].....[HKQR][LAIVFY]', sSeq[:50]):
                iMotifStart = re.search('[RK][LVIQ].....[HKQR][LAIVFY]', sSeq[:50]).start(0)
                fPts2Score = 1
                if re.search('[RK][LI].....[HQ]L', sSeq[:50]):
                    iMotifStart = re.search('[RK][LI].....[HQ]L', sSeq[:50]).start(0)
                    fPts2Score = 2
            sMotif = ''
            if iMotifStart != None:
                sMotif = sSeq[iMotifStart:iMotifStart+9]
            self.dSeqId2Annotation[seq.id]['pts2seq'] = sMotif
            self.dSeqId2Annotation[seq.id]['pts2score'] = fPts2Score
        self.lColumns += ['pts2score', 'pts2seq']
    def runTmhmmProcess(self, iThread):
        fileOut = open(self.sDir+self.sInfilename+'.tmhmm.thread'+str(iThread),'w')
        for seq in SeqIO(open(self.sDir+self.sInfilename+'.thread'+str(iThread))):
            open(self.sDir+'runTmhmm.temp01.'+self.sInfilename+'.thread'+str(iThread),'w').write('>'+seq.id+'\n'+seq.seq+'\n')
            #os.system('cat '+self.sDir+'runTmhmm.temp01.'+self.sInfilename+'.thread'+str(iThread)+' | /home/zarsky/software/tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64 -optionfile /home/zarsky/software/tmhmm-2.0c/lib/TMHMM2.0.options -modelfile /home/zarsky/software/tmhmm-2.0c/lib/TMHMM2.0.model | perl /home/zarsky/software/tmhmm-2.0c/bin/tmhmmformat.pl > '+self.sDir+'runTmhmm.temp02.'+self.sInfilename+'.thread'+str(iThread))
            sNTopology = None
            lTmds = []
            for sLine in open(self.sDir+'runTmhmm.temp02.'+self.sInfilename+'.thread'+str(iThread)):
                if sLine[0] == '#':
                    continue
                (sSeqId, sProgram, sFeature, sStart, sEnd) = sLine.strip().split()
                if sNTopology == None and sFeature in ['inside','outside']:
                    sNTopology = sFeature
                if sFeature == 'TMhelix':
                    lTmds.append( sStart+'-'+sEnd )
            if seq.id != sSeqId:
                raise
            fileOut.write(seq.id+'\t'+str(len(lTmds))+'\t'+sNTopology+'\t'+','.join(lTmds)+'\n')
        fileOut.close()
        os.system('rm -f '+self.sDir+'runTmhmm.temp01.'+self.sInfilename+'.thread'+str(iThread))
        os.system('rm -f '+self.sDir+'runTmhmm.temp02.'+self.sInfilename+'.thread'+str(iThread))
    def runTmhmmArch(self):
        lProcesses = []
        for iThread in range(self.iCpus):
            p = Process(target=self.runTmhmmProcess, args=(iThread,))
            p.start()
            lProcesses.append(p)
        for p in lProcesses:
            p.join()
        fileOut = open(self.sDir+self.sInfilename+'.tmhmm.tsv','w')
        for iThread in range(self.iCpus):
            fileOut.write(open(self.sDir+self.sInfilename+'.tmhmm.thread'+str(iThread)).read())
            os.system('rm '+self.sDir+self.sInfilename+'.tmhmm.thread'+str(iThread))
        fileOut.close()
    def runTmhmm(self, sTmhmmDir='/home/zarsky/software/tmhmm-2.0c/bin/'):
        sMyDir = os.getcwd()
        os.system('cp '+self.sDir+self.sInfilename+' '+sTmhmmDir+self.sInfilename)
        os.chdir(sTmhmmDir)
        os.system('tmhmm '+self.sInfilename+' > '+self.sInfilename+'.tmhmm')
        os.chdir(sMyDir)
        os.system('mv '+sTmhmmDir+self.sInfilename+'.tmhmm '+self.sDir+self.sInfilename+'.tmhmm')
        os.system('rm -f '+sTmhmmDir+self.sInfilename)
        os.system('rm -r -f '+sTmhmmDir+'TMHMM_*')
    def parseTmhmmArch(self, tmhmmfilename=None):
        if tmhmmfilename == None:
            tmhmmfilename = self.sInfilename+'.tmhmm.tsv'
        for sLine in open(self.sDir+tmhmmfilename):
            (sSeqId, sNTmds, sTopology, sTmds) = sLine.strip('\n').split('\t')
            lTmds = sTmds.split(',')
            if lTmds == ['']:
                lTmds = []
            lTmds = map(lambda x:(int(x.split('-')[0]), int(x.split('-')[1])), lTmds)
            self.dSeqId2Annotation[sSeqId]['tmhmmNTmds'] = int(sNTmds)
            self.dSeqId2Annotation[sSeqId]['tmhmmTopology'] = sTopology
            self.dSeqId2Annotation[sSeqId]['tmhmmTmds'] = lTmds
        self.lColumns += ['tmhmmNTmds','tmhmmTmds','tmhmmTopology']
    def parseTmhmm(self, tmhmmfilename=None):
        if tmhmmfilename == None:
            tmhmmfilename = self.sInfilename+'.tmhmm'
        sNTopology = None
        lTmds = []
        lLines = open(self.sDir+tmhmmfilename).read().strip().split('\n') ## not great for super-big files
        for iLine in range(len(lLines)):
            sLine = lLines[iLine]
            if sLine[0] == '#':
                continue
            (sSeqId, sProgram, sFeature, sStart, sEnd) = sLine.strip().split()
            if sNTopology == None and sFeature in ['inside','outside']:
                sNTopology = sFeature
            if sFeature == 'TMhelix':
                lTmds.append( (int(sStart), int(sEnd)) )
            if iLine+1 == len(lLines) or lLines[iLine+1][0] == '#':
                self.dSeqId2Annotation[sSeqId]['tmhmmNTmds'] = len(lTmds)
                self.dSeqId2Annotation[sSeqId]['tmhmmTopology'] = sNTopology
                self.dSeqId2Annotation[sSeqId]['tmhmmTmds'] = lTmds
                sNTopology = None
                lTmds = []
        self.lColumns += ['tmhmmNTmds','tmhmmTmds','tmhmmTopology']
    def runTargetpProcess(self, iThread, sPlant='False'):
        sOrg = 'non-pl'
        if sPlant:
            sOrg = 'pl'
        os.system('/home/zarsky/software/targetp-2.0/bin/targetp -verbose=true -format=short -stdout -batch 1 -org='+sOrg+' -fasta='+self.sDir+self.sInfilename+'.thread'+str(iThread)+' > '+self.sDir+'runTargetp.temp01.'+self.sInfilename+'.thread'+str(iThread))
    def runTargetp(self, sPlant=False):
        lProcesses = []
        for iThread in range(self.iCpus):
            p = Process(target=self.runTargetpProcess, args=(iThread,sPlant))
            p.start()
            lProcesses.append(p)
        for p in lProcesses:
            p.join()
        fileOut = open(self.sDir+self.sInfilename+'.targetp.tsv','w')
        for iThread in range(self.iCpus):
            fileOut.write(open(self.sDir+'runTargetp.temp01.'+self.sInfilename+'.thread'+str(iThread)).read())
            os.system('rm '+self.sDir+'runTargetp.temp01.'+self.sInfilename+'.thread'+str(iThread))
        fileOut.close()
    def parseTargetp(self, targetpfilename=None):
        if targetpfilename == None:
            targetpfilename = self.sInfilename+'.targetp.tsv'
        #dSeqIds = {}
        #for seq in SeqIO(open(self.sInfilename)):
        #    dSeqIds[seq.id] = None
        for sLine in open(self.sDir+targetpfilename):
            if sLine[0] == '#':
                continue
            (sId, sPrediction, sNoTp, sSp, sMtp, sCs) = sLine.strip('\n').split('\t')
            iCs = 0
            fCsProb = 0
            if sCs != '' and '?' not in sCs:
                iCs = int(sCs.split('CS pos:')[1].split('-')[0].strip())
                fCsProb = float(sCs.split('Pr:')[1].strip())
            self.dSeqId2Annotation[sId]['targetp.prediction'] = sPrediction
            self.dSeqId2Annotation[sId]['targetp.notTP'] = float(sNoTp)
            self.dSeqId2Annotation[sId]['targetp.SP'] = float(sSp)
            self.dSeqId2Annotation[sId]['targetp.mTP'] = float(sMtp)
            self.dSeqId2Annotation[sId]['targetp.CS'] = iCs
            self.dSeqId2Annotation[sId]['targetp.CSprob'] = fCsProb
            #del dSeqIds[sId]
        #for sId in dSeqIds.keys():
        #    self.dSeqId2Annotation[sId]['targetp.prediction'] = ''
        #    self.dSeqId2Annotation[sId]['targetp.notTP'] = 0
        #    self.dSeqId2Annotation[sId]['targetp.SP'] = 0
        #    self.dSeqId2Annotation[sId]['targetp.mTP'] = 0
        #    self.dSeqId2Annotation[sId]['targetp.CS'] = 0
        #    self.dSeqId2Annotation[sId]['targetp.CSprob'] = 0
        self.lColumns += ['targetp.prediction', 'targetp.notTP', 'targetp.SP', 'targetp.mTP', 'targetp.CS', 'targetp.CSprob']
    def parseTargetp1(self, targetpfilename=None):
        if targetpfilename == None:
            targetpfilename = self.sInfilename+'.targetp.tsv'
        for sLine in open(self.sDir+targetpfilename):
            (sId, sMtp, sSp, sNoTp, sPrediction, sRs) = sLine.strip('\n').split('\t')
            self.dSeqId2Annotation[sId]['targetp.prediction'] = sPrediction
            self.dSeqId2Annotation[sId]['targetp.notTP'] = float(sNoTp)
            self.dSeqId2Annotation[sId]['targetp.SP'] = float(sSp)
            self.dSeqId2Annotation[sId]['targetp.mTP'] = float(sMtp)
        self.lColumns += ['targetp.prediction', 'targetp.notTP', 'targetp.SP', 'targetp.mTP']
    def runPsortProcess(self, iThread, sPsortDir):
        fileOut = open(self.sDir+'runPsort.temp03.'+self.sInfilename+'.thread'+str(iThread),'w')
        for seq in SeqIO(open(self.sDir+self.sInfilename+'.thread'+str(iThread))):
            if len(seq.seq.strip('*')) <= 20:
                continue
            open(self.sDir+'runPsort.temp01.'+self.sInfilename+'.thread'+str(iThread),'w').write('>'+seq.id+'\n'+seq.seq+'\n')
            os.system('perl '+sPsortDir+'/psort2.'+str(iThread)+'/psort '+self.sDir+'runPsort.temp01.'+self.sInfilename+'.thread'+str(iThread)+' > '+self.sDir+'runPsort.temp02.'+self.sInfilename+'.thread'+str(iThread))
            lFile = open(self.sDir+'runPsort.temp02.'+self.sInfilename+'.thread'+str(iThread)).read().strip('-').strip().split('\n\n')
            sId = lFile[0].strip()
            if sId != seq.id:
                raise
            lAttributes = []
            lTemp = lFile[1].replace(':','').strip().split()
            for iAttribute in range(0, len(lTemp), 2):
                lAttributes.append( (lTemp[iAttribute], float(lTemp[iAttribute+1])) )
            lPredictions = []
            lTemp = lFile[2].strip().split('\n')
            for sPair in lTemp:
                lPair = map(lambda x:x.strip(), sPair.split('%:'))
                lPredictions.append( (lPair[1], float(lPair[0])) )
            fileOut.write( seq.id+'\t'+'\t'.join(map(lambda x:x[0]+':'+str(x[1]), lAttributes+lPredictions))+'\n' )
        fileOut.close()
    def runPsort(self, sPsortDir='/home/zarsky/software/'):
        lProcesses = []
        for iThread in range(self.iCpus):
            os.system('cp -r '+sPsortDir+'/psort2 '+sPsortDir+'/psort2.'+str(iThread))
            p = Process(target=self.runPsortProcess, args=(iThread,sPsortDir))
            p.start()
            lProcesses.append(p)
        for p in lProcesses:
            p.join()
        os.system('rm -r '+sPsortDir+'/psort2.*')
        lAttributes = ['psg', 'gvh', 'alm', 'top', 'tms', 'mit', 'mip', 'nuc', 'erl', 'erm', 'pox', 'px2', 'vac', 'rnp', 'act', 'caa', 'yqr', 'tyr', 'leu', 'gpi', 'myr', 'dna', 'rib', 'bac', 'm1a', 'm1b', 'm2', 'mNt', 'm3a', 'm3b', 'm_', 'ncn', 'lps', 'len']
        lPredictions = ['Golgi', 'cytoplasmic', 'cytoskeletal', 'endoplasmic reticulum', 'extracellular, including cell wall', 'mitochondrial', 'nuclear', 'peroxisomal', 'plasma membrane', 'vacuolar', 'vesicles of secretory system']
        fileOut = open(self.sDir+self.sInfilename+'.psort.tsv','w')
        fileOut.write('id\tprediction\t'+'\t'.join(lPredictions+lAttributes)+'\n')
        for iThread in range(self.iCpus):
            for sLine in open(self.sDir+'runPsort.temp03.'+self.sInfilename+'.thread'+str(iThread)):
                lLine = sLine.strip('\n').split('\t')
                fileOut.write(lLine[0])
                sBestPrediction = ''
                fBestPrediction = 0
                dPredictionsAndAttributes = {}
                for sTemp in lLine[1:]:
                    (sAttr, sValue) = map(lambda x:x.strip(), sTemp.split(':'))
                    dPredictionsAndAttributes[sAttr] = sValue
                for sPrediction in lPredictions:
                    if sPrediction in dPredictionsAndAttributes and float(dPredictionsAndAttributes[sPrediction]) > fBestPrediction:
                        fBestPrediction = float(dPredictionsAndAttributes[sPrediction])
                        sBestPrediction = sPrediction
                fileOut.write('\t'+sBestPrediction)
                for sAttr in lPredictions+lAttributes:
                    fileOut.write('\t')
                    if sAttr in dPredictionsAndAttributes:
                        fileOut.write(dPredictionsAndAttributes[sAttr])
                        del dPredictionsAndAttributes[sAttr]
                if len(dPredictionsAndAttributes) > 0:
                    print(dPredictionsAndAttributes)
                    raise
                fileOut.write('\n')
            os.system('rm -f '+self.sDir+'runPsort.temp*.'+self.sInfilename+'.thread'+str(iThread))
        fileOut.close()
    def parsePsort(self, psortfilename=None):
        if psortfilename == None:
            psortfilename = self.sInfilename+'.psort.tsv'
        lColumns = None
        for sLine in open(self.sDir+psortfilename):
            lLine = sLine.strip('\n').split('\t')
            if lColumns == None:
                lColumns = lLine
                continue
            dLine = dict(zip(lColumns, lLine))
            sId = dLine['id']
            self.dSeqId2Annotation[sId]['psort.prediction'] = dLine['prediction']
            for sAttr in lColumns[2:]:
                if dLine[sAttr] == '':
                    dLine[sAttr] = 0
                self.dSeqId2Annotation[sId]['psort.'+sAttr] = float(dLine[sAttr])
        for sAttr in lColumns[1:]:
            self.lColumns.append('psort.'+sAttr)
    def writeAnnotations(self, lSelectedIdsOnly=None):
        if lSelectedIdsOnly == None:
            tableOut = open(self.sDir+self.sInfilename+'.annotations.tsv','w')
        else:
            tableOut = open(self.sDir+self.sInfilename+'.sel.annotations.tsv','w')
        tableOut.write( '\t'.join( map(lambda x:str(x), self.lColumns) )+'\n' )
        for seq in SeqIO(open(self.sDir+self.sInfilename)):
            if lSelectedIdsOnly != None and seq.id.split('-')[0].strip('A') not in lSelectedIdsOnly: ####!!!!
                continue
            lLineOut = []
            for sColumn in self.lColumns:
                if sColumn in self.dSeqId2Annotation[seq.id]:
                    if type(self.dSeqId2Annotation[seq.id][sColumn]) == list:
                        if len(self.dSeqId2Annotation[seq.id][sColumn]) > 0 and type(self.dSeqId2Annotation[seq.id][sColumn][0]) in (list, tuple):
                            lLineOut.append( '; '.join( map(lambda x:', '.join(map(lambda y:str(y), x)), self.dSeqId2Annotation[seq.id][sColumn]) ) )
                        else:
                            lLineOut.append( '; '.join( map(lambda x:str(x), self.dSeqId2Annotation[seq.id][sColumn]) ) )
                    else:
                        lLineOut.append( str(self.dSeqId2Annotation[seq.id][sColumn]) )
                else:
                    lLineOut.append( '' )
            tableOut.write( '\t'.join(lLineOut)+'\n' )
        tableOut.close()
    def writeFastaAnnotations(self):
        fileOut = open(self.sDir+self.sInfilename.rsplit('.',1)[0]+'.annotations.fasta','w')
        ## adjust !!!
        self.lColumns.remove( 'id' )
        self.lColumns.remove( 'desc' ) ##
        #self.lColumns.remove( 'gos' )
        for seq in SeqIO(open(self.sInfilename)):
            lLineOut = []
            for sColumn in self.lColumns:
                if sColumn in self.dSeqId2Annotation[seq.id]:
                    if type(self.dSeqId2Annotation[seq.id][sColumn]) == list:
                        lLineOut.append( sColumn+'='+'; '.join( map(lambda x:str(x), self.dSeqId2Annotation[seq.id][sColumn]) ) )
                    else:
                        lLineOut.append( sColumn+'='+str(self.dSeqId2Annotation[seq.id][sColumn]) )
            ## specific for transdecoder/selenidium !!!
            #seq.id = seq.id.split(':')[0]
            #sSource = seq.desc.strip().split()[-1]
            #seq.desc = ' '.join( lLineOut )+' source='+sSource
            ## general:
            seq.desc = ' '.join( lLineOut )
            fileOut.write('>'+seq.id+' '+seq.desc+'\n'+seq.seq+'\n')
        fileOut.close()
    def clean(self):
        for iThread in range(self.iCpus):
            os.system('rm -f '+self.sDir+self.sInfilename+'.thread'+str(iThread))


"""

#lEntamoebas = ['eh','en','ed','em','ei']
lEntamoebas = ['eh','en','ed','em','ei']
for sEntamoeba in lEntamoebas:
    print sEntamoeba
    annotation = Annotation(sEntamoeba+'.fa', cpus=1)
    annotation.parseSeqProperties(parseIdFunction=lambda x:x.split('-')[0], parseDescFunction=lambda x:x.split(' gene_product=')[1].split('|')[0].strip())
    #annotation.runEggnog() #
    #annotation.runKofam()
    #annotation.runTmhmm() #
    #annotation.runTargetp() #
    #annotation.runPsort() #
    annotation.readKegg()
    annotation.parseEggnog()
    annotation.parseKofamFromEggnog()
    annotation.parseTmhmm()
    annotation.predictPts1()
    annotation.parseTargetp()
    annotation.parsePsort()
    #print annotation.dSeqId2Annotation['EHI7A_110630-t26_1-p1']
    annotation.writeAnnotations()
    #annotation.clean()
"""

