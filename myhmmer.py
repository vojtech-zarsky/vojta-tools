import os
import re
from tkinter import *
import tkinter.filedialog
#import tkFileDialog

#DB = 'mastiga_aa.fasta'
import re

class Seq:
    def __init__(self, sId, sDesc, sSeq):
        self.id = sId
        self.desc = ''
        if len(re.split('\s+', sDesc, 1)) > 1:
            self.desc = re.split('\s+', sDesc, 1)[1]
        self.seq = sSeq
    def getString(self):
        return '>'+self.id+' '+self.desc+'\n'+self.seq+'\n'

class SeqIO:
    def __init__(self, file):
        self.file = file
        self.end = 1
        sLast = ''
        for sLine in self.file:
            sLine = sLine.strip('\n')
            if re.search('^>', sLast) and len(sLine) > 0 and sLine[0] != '\s':
                self.header = sLast.lstrip('>').strip('\n')
                self.seq = sLine.strip()
                self.end = 0
                break
            sLast = sLine
    def __iter__(self):
        return self
    def __next__(self):
        return self.next()
    def end(self):
        return self.end
    def next(self):
        if self.end:
            raise StopIteration()
        sLine = ''
        while True:
            try:
                sLine = next(self.file).strip('\n')
                if re.search('^>',sLine):
                    seq = Seq(self.header.split()[0], self.header, self.seq)
                    self.header = sLine.lstrip('>').strip('\n')
                    self.seq = ''
                    return seq
                self.seq = self.seq+sLine.strip()
            except StopIteration:
                seq = Seq(self.header.split()[0], self.header, self.seq)
                self.header = None
                self.seq = None
                self.end = 1
                return seq
                break
        return 0
    def getDict(self):
        dSeqs = {}
        while not self.end:
            seq = self.next()
            dSeqs[seq.id] = seq
        return dSeqs
    def writeIndex(self, sDb=None):
        handle = open(self.file.name)
        fileOut = open(handle.name+'.index','w')
        iIndex = handle.tell()
        sLine = handle.readline()
        if sDb == 'uniprot':
            while sLine:
                if len(sLine) > 1 and sLine[0] == '>':
                    fileOut.write(sLine.split()[0].split('|')[2]+'\t'+str(iIndex)+'\n')
                iIndex = handle.tell()
                sLine = handle.readline()
        else:
            while sLine:
                if len(sLine) > 1 and sLine[0] == '>':
                    fileOut.write(sLine.split()[0][1:]+'\t'+str(iIndex)+'\n')
                iIndex = handle.tell()
                sLine = handle.readline()
        return 0
    def readIndex(self):
        self.dId2Index = {}
        for sLine in open(self.file.name+'.index'):
            (sId, sIndex) = sLine.strip().split('\t',1)
            self.dId2Index[sId] = int(sIndex)
        return 0
    def getSeq(self, sId, iIndex=None):
        if iIndex == None:
            iIndex = self.dId2Index[sId]
        self.file.seek( iIndex )
        sHeader = self.file.readline()[1:].strip()
        sSeq = ''
        for sLine in self.file:
            sLine = sLine.strip()
            if sLine == '' or sLine[0] == '>':
                break
            sSeq += sLine
        return Seq(sHeader.split()[0], sHeader, sSeq)


def sortE(dX,dY):
    if dX['e']<=dY['e']:
        return -1
    else:
        return 1
        
class ParseHMM:
    def __init__(self,fileIn):
        self.dResult = {'hits':[]}
        sFile = fileIn.read()
        self.dResult['file'] = sFile
        self.dResult['queryHMM'] = sFile.split('# query HMM file:')[1].split('\n')[0].strip()
        self.dResult['database'] = sFile.split('# target sequence database:')[1].split('\n')[0].strip()
        for sHit in sFile.split('\n>>')[1:]:
            dHit = {}
            dHit['id'] = sHit.strip().split()[0]
            dHit['desc'] = sHit.strip().split('\n')[0].split(' ',1)[-1]
            lStats = list(map( lambda x:x.strip().split(), sHit.split('\n\n')[0].split('\n')[3:] ))
            dHit['hsps'] = []
            iHsp = 0
            for sHsp in sHit.split('\n  ==')[1:]:
                dHsp = {}
                lStatsLine = lStats[iHsp]
                dHsp['score'] = lStatsLine[2]
                dHsp['e'] = lStatsLine[5]
                dHsp['queryStart'] = lStatsLine[6]
                dHsp['queryEnd'] = lStatsLine[7]
                dHsp['envStart'] = lStatsLine[12]
                dHsp['envEnd'] = lStatsLine[13]
                dHsp['hsp'] = ''
                lHsp = sHsp.split('\n\n\n')[0].split('\n')
                iLineIndex = 3
                while iLineIndex < len(lHsp):
                    try:
                        dHsp['hsp'] += lHsp[iLineIndex].strip().split()[2].replace('-','').upper()
                    except:
                        print(lHsp[iLineIndex])
                    iLineIndex += 5
                dHit['hsps'].append( dHsp )
                iHsp += 1
            if len(dHit['hsps']) == 0:
                continue
            #dHit['e'] = sorted( map(lambda x:x['e'], dHit['hsps']), key=lambda x:float(x) )[0]
            dHit['bestHsp'] = dHit['hsps'].index( sorted( dHit['hsps'], key=lambda x:float(x['e']) )[0] )
            dHit['e'] = dHit['hsps'][dHit['bestHsp']]['e']
            self.dResult['hits'].append(dHit)
        #lFile = sFile.split('\n')
        """
        i = -1
        for sLine in lFile:
            i = i+1
            if re.search('^# query HMM file:',sLine):
                self.dResult['queryHMM'] = sLine.split(' ')[-1].strip()
            if re.search('^# target sequence database:',sLine):
                self.dResult['database'] = sLine.split('target sequence database:')[1].strip()
            if re.search('^//',sLine) and 'id' in dHit:
                self.dResult['hits'].append(dHit)
                dHit = {'hsps':[]}
            if re.search('^>>',sLine) and lFile[i+2] != '':
                iDomains = 0
                if 'id' in dHit: 
                    self.dResult['hits'].append(dHit)
                    dHit = {'hsps':[]}
                dHit['id'] = sLine.lstrip('> ').split(' ')[0].strip()
                dHit['desc'] = sLine.lstrip('> ').lstrip(dHit['id']).strip()
                for j in range(i+3,len(lFile)):
                    if lFile[j] == '' or lFile[j] == '^>>':
                        break
                    lLine = lFile[j].strip()
                    lLine = re.split(' +',lLine)
                    dHit['hsps'].append({'score':lLine[2],'e':lLine[5],'queryStart':lLine[6],'queryEnd':lLine[7],'aliStart':lLine[9],'aliEnd':lLine[10],'envStart':lLine[12],'envEnd':lLine[13],'hsp':''})
                    if 'bestHsp' in dHit:
                        if lLine[5] < dHit['hsps'][dHit['bestHsp']]['e']:
                            dHit['bestHsp'] = iDomains
                    else:
                        dHit['bestHsp'] = iDomains
                    if 'e' in dHit:
                        dHit['e'] = min(dHit['e'],lLine[5])
                    else:
                        dHit['e'] = lLine[5]
                    iDomains = iDomains+1
                k = 0
                for j in range(i+iDomains+3,len(lFile)):
                    if re.search('^>>',lFile[j]):
                        break
                    if re.search('== domain '+str(k+1),lFile[j]):
                        if k >= iDomains:
                            break
                        for l in range(j+1,len(lFile)):
                            if re.search('== domain ',lFile[l]) or re.search('^>>',lFile[l]):
                                break
                            if re.search(dHit['id'].replace('|','\|'),lFile[l]):
                                sTemp = lFile[l].strip()
                                lTemp = re.split(' +',sTemp)       
                                dHit['hsps'][k]['hsp'] = dHit['hsps'][k]['hsp']+lTemp[2].strip()
                        k = k+1
        """
    def getHits(self):
        return self.dResult

def runHMM(sIn,dParams):
    os.system('rm -f temp_runHMM*')
    if sIn.strip()[:6] == 'HMMER3':
        open('temp_runHMM02','w').write(sIn)
    else:
        if not re.search('^\s*>', sIn):
            sIn = '>name\n'+sIn
        open('temp_runHMM01','w').write(str(sIn))
        os.system('hmmbuild --informat afa temp_runHMM02 temp_runHMM01')
    sParams = ''
    if 'E' in dParams:
        sParams = sParams+' -E '+str(dParams['E'])
    if 'T' in dParams:
        sParams = sParams+' -T '+str(dParams['T'])
    if 'Z' in dParams:
        sParams = sParams+' -Z '+str(dParams['Z'])
    os.system('hmmsearch -E 10000 --cpu 2 -o temp_runHMM03 '+str(sParams)+' temp_runHMM02 "'+str(dParams['db'])+'"')
    
    return ParseHMM(open('temp_runHMM03')).getHits()

def getSeq(sFile, lIds):
        lSeqs = []
        seqio = SeqIO(open(sFile))
        while not seqio.end:
            seq = seqio.next()
            if seq.id in lIds:
                lSeqs.append({'id':seq.id, 'desc':seq.desc, 'seq':seq.seq})
                lIds.remove(seq.id)
                if len(lIds) == 0:
                    break
        if len(lIds) > 0:
            print('!! ids not found', lIds)
        return lSeqs


class SeqWindow:
    def __init__(self, lSeq, sName):
        seqWindow = Toplevel()
        seqWindow.title(sName)
        mainFrame = Frame(seqWindow)
        mainFrame.pack(fill=BOTH, expand=1)
        scrollbar = Scrollbar(mainFrame)
        scrollbar.pack(side=RIGHT, fill=Y)
        seqText = Text(mainFrame, yscrollcommand=scrollbar.set)
        seqText.pack(fill=BOTH, expand=1)
        for dSeq in lSeq:
            seqText.insert(END, '>'+str(dSeq['id'])+' '+dSeq['desc']+'\n'+str(dSeq['seq'])+'\n')
            scrollbar.config(command=seqText.yview)
    
class ResultWindow:
    def __init__(self, dResult):
        self.dResult = dResult
        resultWindow = Toplevel()
        resultWindow.geometry("1050x700")
        resultWindow.title('Results')
        mainFrame = Frame(resultWindow)
        mainFrame.pack(fill=BOTH, expand=1) 
        scrollbar = Scrollbar(mainFrame)
        scrollbar.pack(side=RIGHT, fill=Y) 
        self.resultText = Text(mainFrame, yscrollcommand=scrollbar.set)
        self.resultText.pack(fill=BOTH, expand=1)
        self.resultText.insert(END, dResult['file'])
        
        scrollbar.config(command=self.resultText.yview)
        hitsFrame = Frame(resultWindow)
        hitsFrame.pack(fill=BOTH, expand=1)
        scrollbar = Scrollbar(hitsFrame)
        scrollbar.pack(side=RIGHT, fill=Y)
        self.listboxHits = Listbox(hitsFrame, selectmode=EXTENDED, yscrollcommand=scrollbar.set)
        for dHit in dResult['hits']:
            self.listboxHits.insert(END, dHit['id']+' ('+str(dHit['e'])+')')
        self.listboxHits.pack(fill=BOTH, expand=1)
        self.listboxHits.select_set(0)
        scrollbar.config(command=self.listboxHits.yview)
        items = self.listboxHits.curselection()
        
        
        buttonFrame = Frame(resultWindow)
        buttonFrame.pack(side=LEFT)
        button = Button(buttonFrame, text='  Hits  ',command=self.getHit)
        button.pack(side=LEFT)
        button = Button(buttonFrame, text='  Hsps  ',command=self.getHsp)
        button.pack(side=LEFT)
        button = Button(buttonFrame, text='Select All',command=self.selectAll)
        button.pack(side=LEFT)
    def getHit(self):
        lIds = []
        lHits = self.listboxHits.curselection()
        for i in lHits:
            lIds.append(self.dResult['hits'][int(i)]['id'])
        lSeqs = getSeq(self.dResult['database'],lIds)
        lSeqsOut = []
        for dHit in self.dResult['hits']:
            for dSeq in lSeqs:
                if dHit['id'] == dSeq['id']:
                    dSeq['seq'] = dSeq['seq'].replace('*','')
                    lSeqsOut.append(dSeq)
        SeqWindow(lSeqsOut,'Hits')
    def getHsp(self):
        lIds = []
        lHits = self.listboxHits.curselection()
        for i in lHits:
            lIds.append(self.dResult['hits'][int(i)]['id'])
        lSeqs = []
        for dHit in self.dResult['hits']:
            if dHit['id'] in lIds:
                lSeqs.append({'id':dHit['id'],'desc':dHit['desc'],'seq':dHit['hsps'][dHit['bestHsp']]['hsp'].replace('-','')})
        SeqWindow(lSeqs,'Hsps')
    def selectAll(self):
        self.listboxHits.select_set(0,self.listboxHits.size())
        
class InputWindow:
    def __init__(self, master):
        mainFrame = Frame(master)
        mainFrame.pack(fill=BOTH, expand=1)
        scrollbar = Scrollbar(mainFrame)
        scrollbar.pack(side=RIGHT, fill=Y)
        self.inputText = Text(mainFrame, yscrollcommand=scrollbar.set)
        self.inputText.pack(fill=BOTH, expand=1)
        #self.inputText.insert(END, '>cpn\nQEATVIAVGPGRF')
        scrollbar.config(command=self.inputText.yview)

        buttonFrame = Frame(mainFrame)
        buttonFrame.pack(side=LEFT)
        self.searchButton = Button(buttonFrame, text='Search',command=self.search)
        self.searchButton.pack(side=LEFT)
        self.clearButton = Button(buttonFrame, text='Clear',command=self.clear)
        self.clearButton.pack(side=LEFT)
        self.quitButton = Button(buttonFrame, text='Quit',command=mainFrame.quit)
        self.quitButton.pack(side=LEFT)
        paramFrame = Frame(mainFrame)
        paramFrame.pack(side=LEFT)
        labelE = Label(paramFrame, text=' E')
        labelE.pack(side=LEFT)
        self.entryE = Entry(paramFrame,width=5)
        self.entryE.pack(side=LEFT)
        labelT = Label(paramFrame, text=' T')
        labelT.pack(side=LEFT)
        self.entryT = Entry(paramFrame,width=5)
        self.entryT.pack(side=LEFT)
        labelZ = Label(paramFrame, text=' Z')
        labelZ.pack(side=LEFT)
        self.entryZ = Entry(paramFrame,width=5)
        self.entryZ.pack(side=LEFT)
        self.database = ''
        self.dbButton = Button(paramFrame, text='Database', command=self.selectDatabase)
        self.dbButton.pack(side=LEFT)
        
    def clear(self):
        self.inputText.delete(1.0, END)
    def search(self):
        sInput = self.inputText.get(1.0, END)
        dParams = {}
        sE = self.entryE.get()
        sT = self.entryT.get()
        sZ = self.entryZ.get()
        if sE:
            dParams['E'] = sE
        if sT:
            dParams['T'] = sT
        if sZ:
            dParams['Z'] = sZ
        dParams['db'] = self.database
        #dParams['db'] = DB
        if dParams['db'] != '':
                dResult = runHMM(sInput,dParams)
                ResultWindow(dResult)
    def selectDatabase(self):
        self.database = tkinter.filedialog.askopenfilename()
        
root = Tk()
root.title('MyHMMER')
inputWindow = InputWindow(root)
root.mainloop()
