from operator import itemgetter
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
import pandas as pd
from Bio import pairwise2
import numpy as np
import time
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalwCommandline
import os
import datetime
import statistics as stat
import seaborn as sns
import matplotlib.pyplot as plt
import math

class SequenceHelper:
    def __init__(self):
        self.AccessionIDs = []
        self.sequenceData = []
        self.combined = []
        self.alignedSequences = []
        pass

    def sequenceReader(self, filename):
        sequences = list(SeqIO.parse(filename, "fasta"))
        for i in sequences:
            self.AccessionIDs.append(i.id)
            self.sequenceData.append(i.seq)
            self.combined.append([i.id, i.seq])
    
    def MetaClass(self, metaDF):
        patientStatus = np.array(metaDF['Patient status'].values)
        patientStatus1 = np.array(metaDF['Patient status'].values)
        for i in range(0, len(patientStatus)):
            if not(isinstance(patientStatus[i], str)):
                if math.isnan(patientStatus[i]):
                    patientStatus[i] = 'unknown'

        for i in range(0, len(patientStatus1)):
            if not(isinstance(patientStatus1[i], str)):
                if math.isnan(patientStatus1[i]):
                    patientStatus1[i] = 'unknown'

        patientStatus = list(np.intersect1d(patientStatus, patientStatus))
        SymptomClassified = [[], [], []]
        OutcomeClassified = [[], [], []]
        for i in patientStatus:
            print('\n\nSymptom Classifier')
            print('1 == Unknown', '2 == Symptomatic', '3 == Asymptomatic')
            classify = int(input('"' + i + '"')) - 1
            SymptomClassified[classify].append(i)
            print('\nOutcome Classifier')
            print('1 == Unknown', '2 == Recovered', '3 == Deceased')
            classify = int(input('"' + i + '"')) - 1
            OutcomeClassified[classify].append(i)
        

        FinalClassified = []
        FinalClassified2 = []
        for i in patientStatus1:
            if i in SymptomClassified[0]:
                FinalClassified.append('Unknown')
            elif i in SymptomClassified[1]:
                FinalClassified.append('Symptomatic')
            elif i in SymptomClassified[2]:
                FinalClassified.append('Asymptomatic')
            else:
                print('BUG')

            if i in OutcomeClassified[0]:
                FinalClassified2.append('Unknown')
            elif i in OutcomeClassified[1]:
                FinalClassified2.append('Recovered')
            elif i in OutcomeClassified[2]:
                FinalClassified2.append('Deceased')
            else:
                print('BUG')
        
        metaDF['Symptoms'] = FinalClassified
        metaDF['Outcome'] = FinalClassified2

        return metaDF


    def metricCalc(self, refSeq = None):
        if refSeq == None:
    
            referenceSequence = np.array(self.sequenceData[0])
            refLen = len(referenceSequence)
            allSequenceMetrics = []

            refACount = len(list(np.argwhere(referenceSequence == 'A').flatten()))
            refCCount = len(list(np.argwhere(referenceSequence == 'C').flatten()))
            refGCount = len(list(np.argwhere(referenceSequence == 'G').flatten()))
            refTCount = len(list(np.argwhere(referenceSequence == 'T').flatten()))

            for i in self.sequenceData:
                i = np.array(i)
                flag = True
                polyA = 0
                tmp = -1
                while flag == True:
                    if i[tmp:tmp+1] == 'A':
                        polyA+=1
                        tmp -=1
                    else:
                        flag = False

                correctA = len(list(np.argwhere(i == 'A').flatten()))
                correctC = len(list(np.argwhere(i == 'C').flatten()))
                correctG = len(list(np.argwhere(i == 'G').flatten()))
                correctT = len(list(np.argwhere(i == 'T').flatten()))

                numberDash = len(list(np.argwhere(i == '-').flatten()))
                numberN = len(list(np.argwhere(i == 'N').flatten()))
                numberY = len(list(np.argwhere(i == 'Y').flatten()))
                numberM = len(list(np.argwhere(i == 'M').flatten()))
                numberR = len(list(np.argwhere(i == 'R').flatten()))

                total = correctA+correctC+correctG+correctT
                # print(total)

                totalCorrect = len(list(np.intersect1d(referenceSequence, i).flatten()))
                seqlength = len(i)
                percentSimilar = totalCorrect / seqlength

                if total != 0:
                    allSequenceMetrics.append([correctA, (correctA / total), correctC, (correctC / total), correctG, (correctG / total), correctT, (correctT / total), 
                                            total, (total / refLen), numberDash, (numberDash / refLen), numberN, (numberN / refLen), (numberN + numberDash), 
                                            ((numberN + numberDash) / refLen), numberY, numberM, numberR, polyA, percentSimilar, total, correctA - refACount,
                                            correctC - refCCount, correctG - refGCount, correctT - refTCount])
                else:
                    allSequenceMetrics.append([correctA, 0, correctC, 0, correctG, 0, correctT, 0, 
                                            total, (total / refLen), numberDash, (numberDash / refLen), numberN, (numberN / refLen), (numberN + numberDash), 
                                            ((numberN + numberDash) / refLen), numberY, numberM, numberR, polyA, percentSimilar, total, correctA - refACount,
                                            correctC - refCCount, correctG - refGCount, correctT - refTCount])

            return allSequenceMetrics

        else:
            referenceSequence = np.array(refSeq)
            refLen = len(referenceSequence)
            allSequenceMetrics = []

            refACount = len(list(np.argwhere(referenceSequence == 'A').flatten()))
            refCCount = len(list(np.argwhere(referenceSequence == 'C').flatten()))
            refGCount = len(list(np.argwhere(referenceSequence == 'G').flatten()))
            refTCount = len(list(np.argwhere(referenceSequence == 'T').flatten()))
            
            for i in self.sequenceData:
                i = np.array(i)
                flag = True
                polyA = 0
                tmp = -1
                while flag == True:
                    if i[tmp:tmp+1] == 'A':
                        polyA+=1
                        tmp -=1
                    else:
                        flag = False

                correctA = len(list(np.argwhere(i == 'A').flatten()))
                correctC = len(list(np.argwhere(i == 'C').flatten()))
                correctG = len(list(np.argwhere(i == 'G').flatten()))
                correctT = len(list(np.argwhere(i == 'T').flatten()))

                numberDash = len(list(np.argwhere(i == '-').flatten()))
                numberN = len(list(np.argwhere(i == 'N').flatten()))
                numberY = len(list(np.argwhere(i == 'Y').flatten()))
                numberM = len(list(np.argwhere(i == 'M').flatten()))
                numberR = len(list(np.argwhere(i == 'R').flatten()))

                total = correctA+correctC+correctG+correctT

                totalCorrect = len(list(np.intersect1d(referenceSequence, i).flatten()))
                percentSimilar = totalCorrect / refLen

                # print(total)
                
                if total != 0:
                    allSequenceMetrics.append([correctA, (correctA / total), correctC, (correctC / total), correctG, (correctG / total), correctT, (correctT / total), 
                                            total, (total / refLen), numberDash, (numberDash / refLen), numberN, (numberN / refLen), (numberN + numberDash), 
                                            ((numberN + numberDash) / refLen), numberY, numberM, numberR, polyA, percentSimilar, total, correctA - refACount,
                                            correctC - refCCount, correctG - refGCount, correctT - refTCount])
                else:
                    allSequenceMetrics.append([correctA, 0, correctC, 0, correctG, 0, correctT, 0, 
                                            total, (total / refLen), numberDash, (numberDash / refLen), numberN, (numberN / refLen), (numberN + numberDash), 
                                            ((numberN + numberDash) / refLen), numberY, numberM, numberR, polyA, percentSimilar, total, correctA - refACount,
                                            correctC - refCCount, correctG - refGCount, correctT - refTCount])

            return allSequenceMetrics

    
    def align2Sequences(self, refseq, seq2):
        tmp = str('>Reference_Sequence\n' + refseq + '\n')
        tmp1 = str('>Sequence2\n' + seq2)
        if os.path.exists('2seqs.fasta'):
            os.remove('2seqs.fasta')
        f = open("2seqs.fasta", "a")
        f.write(tmp)
        f.write(tmp1)
        f.close()
        start = time.time()
        alignment = MuscleCommandline(cmd='muscle', input='2seqs.fasta', out='2seqsout.fasta')
        alignment('2seqs.fasta')
        end = time.time()
        print(end - start, 'seconds to align')
        print('Aligned\n')
        sequences = list(SeqIO.parse("2seqsout.fasta", "fasta"))
        os.remove('2seqsout.fasta')
        os.remove('2seqs.fasta')

        alignedRef = str(sequences[0].seq)
        alignedSeq2 = list(str(sequences[1].seq))

        dashLoc = list(np.argwhere(alignedRef == '-').flatten())
        dashLoc.reverse()

        for i in dashLoc:
            alignedSeq2.pop(i)
        
        tmpList = ''
        tmpList.join(alignedSeq2)
        alignedSeq2 = tmpList

        self.alignedSequences.append(alignedSeq2)

    def writeAligned(self, outputfile):
        f = open(outputfile, 'w')

        for i in range(0, len(self.AccessionIDs)):
            string1 = '>' + self.AccessionIDs[i] + '\n'
            string2 = str(self.alignedSequences[i]) + '\n'

            f.write(string1)
            f.write(string2)

        print("Done")
    
    def seperateMonth(self, metadataDF, metricDF):
        sequences = metricDF.index
        weeks = []
        months = []
        Accessions = np.array(metadataDF.index)
        metaDates = metadataDF['Collection date'].values
        for i in sequences:
            index = int(np.argwhere(Accessions == i).flatten()[0])
            date = metaDates[index]
            splitDate = date.split('-')
            # print(splitDate)
            if len(splitDate) == 3:
                weekNumber = int(datetime.date(int(splitDate[0]), int(splitDate[1]), int(splitDate[2])).isocalendar()[1])

                months.append(splitDate[1])
                weeks.append(weekNumber)
            
            else:
                months.append(-1)
                weeks.append(-1)

        return months, weeks
    
    def allPercentages(self, allData):
        APercent = stat.mean(allData['A Percent'].values)
        CPercent = stat.mean(allData['C Percent'].values)
        GPercent = stat.mean(allData['G Percent'].values)
        TPercent = stat.mean(allData['T Percent'].values)

        plt.figure()
        ax = sns.barplot(x = ['A', 'C', 'G', 'T'], y = [APercent, CPercent, GPercent, TPercent])
        ax.set_xticklabels(['A', 'C', 'G', 'T'])
        plt.ylabel("Average Percentages")
        plt.savefig("AllPercentages.png", bbox_inches="tight", pad_inches=1)
        plt.clf()
    
    def monthPercentages(self, allData, refPercent):
        monthValues = allData['Month'].values
        AValues = allData['A Percent'].values
        CValues = allData['C Percent'].values
        GValues = allData['G Percent'].values
        TValues = allData['T Percent'].values

        WholeAPercent = stat.mean(allData['A Percent'].values)
        WholeCPercent = stat.mean(allData['C Percent'].values)
        WholeGPercent = stat.mean(allData['G Percent'].values)
        WholeTPercent = stat.mean(allData['T Percent'].values)
        
        WholePercent = [[WholeAPercent], [WholeCPercent], [WholeGPercent], [WholeTPercent]]

        perMonth = []
        for i in range(0, 12):
            perMonth.append([[], [], [], []])

        Deleted = 0

        for i in range(0, len(monthValues)):
            if monthValues[i] != -1:
                perMonth[monthValues[i] - 1][0].append(AValues[i])
                perMonth[monthValues[i] - 1][1].append(CValues[i])
                perMonth[monthValues[i] - 1][2].append(GValues[i])
                perMonth[monthValues[i] - 1][3].append(TValues[i])
            
            else:
                Deleted+=1

        previousYear = []
        refPercent = [[refPercent[0]], [refPercent[1]], [refPercent[2]], [refPercent[3]]]
        # WholePercent = [[WholePercent[0]], [WholePercent[1]], [WholePercent[2]], [WholePercent[3]]]
        previousYear.append(refPercent)
        previousYear.append(WholePercent)
        previousYear.append(perMonth[-1])
        thisYear = perMonth[:-3]

        for i in thisYear:
            previousYear.append(i)

        perMonth = previousYear

        finalValues = []
        for i in perMonth:
            if len(i[0]) == 0:
                finalValues.append([0, 0, 0, 0])
            else:
                AverageA = stat.mean(i[0])
                AverageC = stat.mean(i[1])
                AverageG = stat.mean(i[2])
                AverageT = stat.mean(i[3])
                finalValues.append([AverageA, AverageC, AverageG, AverageT])

        ABar = []
        CBar = []
        GBar = []
        TBar = []

        width = 0.35
        ind = np.arange(12)

        for i in finalValues:
            ABar.append(i[0])
            CBar.append(i[1])
            GBar.append(i[2])
            TBar.append(i[3])

        bars = np.add(ABar, CBar).tolist()
        bars1 = np.add(bars, GBar).tolist()

        p1 = plt.bar(ind, ABar, width)
        p2 = plt.bar(ind, CBar, width, bottom=ABar)
        p3 = plt.bar(ind, GBar, width, bottom=bars)
        p4 = plt.bar(ind, TBar, width, bottom=bars1)

        xlabelticks = ['December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'Reference Sequence','All Sequences']

        xticklabels = ['Reference Sequence','All Sequences', 'December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September']
        xticklabels0 = ['Reference Sequence','All Sequences', 'December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September']
        xticklabels1 = ['Reference Sequence','All Sequences', 'December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September']
        xticklabels2 = ['Reference Sequence','All Sequences', 'December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September']

        print('Number of Sequences without a Month:', Deleted)

        plt.ylabel('Percentages')
        plt.title('Percentages by Month')
        plt.xticks(ind, ('Reference Sequence','All Sequences', 'December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September'), rotation=90)
        # plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A Percent', 'C Percent', 'G Percent', 'T Percent'))
        plt.savefig("PerMonthPercentages.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(ABar)):
            if ABar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del ABar[i]
                del xticklabels[i]

        ind = np.arange(len(xticklabels))

        ref1 = []
        ref2 = []
        ref3 = []
        ref4 = []
        for i in range(0, len(xticklabels)):
            ref1.append(refPercent[0][0])
            ref2.append(refPercent[1][0])
            ref3.append(refPercent[2][0])
            ref4.append(refPercent[3][0])

        plt.figure()
        sns.lineplot(x = xticklabels, y = ABar)
        sns.lineplot(x = xticklabels, y = ref1)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage A')
        plt.xlabel('Month Number')
        plt.title('Percentages By Month')
        plt.savefig('AMonthLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(CBar)):
            if CBar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del CBar[i]
                del xticklabels0[i]

        ind = np.arange(len(xticklabels0))

        plt.figure()
        sns.lineplot(x = xticklabels0, y = CBar)
        sns.lineplot(x = xticklabels0, y = ref2)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage C')
        plt.xlabel('Month Number')
        plt.title('Percentages By Month')
        plt.savefig('CMonthLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(GBar)):
            if GBar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del GBar[i]
                del xticklabels1[i]

        ind = np.arange(len(xticklabels1))

        plt.figure()
        sns.lineplot(x = xticklabels1, y = GBar)
        sns.lineplot(x = xticklabels1, y = ref3)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage G')
        plt.xlabel('Month Number')
        plt.title('Percentages By Month')
        plt.savefig('GMonthLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(TBar)):
            if TBar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del TBar[i]
                del xticklabels2[i]

        ind = np.arange(len(xticklabels2))

        plt.figure()
        sns.lineplot(x = xticklabels2, y = TBar)
        sns.lineplot(x = xticklabels2, y = ref4)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage T')
        plt.xlabel('Month Number')
        plt.title('Percentages By Month')
        plt.savefig('TMonthLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        Bars = [ABar, CBar, GBar, TBar]

        monthDF = pd.DataFrame(Bars, index = ['A', 'C', 'G', 'T'], columns = (xticklabels))
        monthDF.to_csv('MonthValues.csv')

    def WeekPercentage(self, allData, refPercent):
        weekNumbers = allData['Week Number'].values
        values = []
        Deleted = 0
        for i in range(0, 52):
            values.append([[], [], [], []])
        
        AValues = allData['A Percent'].values
        CValues = allData['C Percent'].values
        GValues = allData['G Percent'].values
        TValues = allData['T Percent'].values

        for i in range(0, len(weekNumbers)):
            if weekNumbers[i] != -1:
                values[weekNumbers[i] - 1][0].append(AValues[i])
                values[weekNumbers[i] - 1][1].append(CValues[i])
                values[weekNumbers[i] - 1][2].append(GValues[i])
                values[weekNumbers[i] - 1][3].append(TValues[i])
            
            else:
                Deleted += 1
        
        Averages = []
        for i in values:
            if len(i[0]) != 0:
                AAverage = stat.mean(i[0])
                CAverage = stat.mean(i[1])
                GAverage = stat.mean(i[2])
                TAverage = stat.mean(i[3])
            else:
                AAverage = 0
                CAverage = 0
                GAverage = 0
                TAverage = 0
            Averages.append([AAverage, CAverage, GAverage, TAverage])
        
        WholeAPercent = stat.mean(allData['A Percent'].values)
        WholeCPercent = stat.mean(allData['C Percent'].values)
        WholeGPercent = stat.mean(allData['G Percent'].values)
        WholeTPercent = stat.mean(allData['T Percent'].values)

        WholePercent = [WholeAPercent, WholeCPercent, WholeGPercent, WholeTPercent]

        newvals = [WholePercent, refPercent]
        for i in Averages:
            newvals.append(i)

        ABar = []
        CBar = []
        GBar = []
        TBar = []

        width = 0.35
        ind = np.arange(39)

        for i in newvals[:39]:
            ABar.append(i[0])
            CBar.append(i[1])
            GBar.append(i[2])
            TBar.append(i[3])

        bars = np.add(ABar, CBar).tolist()
        bars1 = np.add(bars, GBar).tolist()

        p1 = plt.bar(ind, ABar, width)
        p2 = plt.bar(ind, CBar, width, bottom=ABar)
        p3 = plt.bar(ind, GBar, width, bottom=bars)
        p4 = plt.bar(ind, TBar, width, bottom=bars1)

        xticklabels = ['Reference Sequence', 'All Sequences']
        xticklabels0 = ['Reference Sequence', 'All Sequences']
        xticklabels1 = ['Reference Sequence', 'All Sequences']
        xticklabels2 = ['Reference Sequence', 'All Sequences']
        for i in range(1, 38):
            xticklabels.append(str(i))
            xticklabels0.append(str(i))
            xticklabels1.append(str(i))
            xticklabels2.append(str(i))

        xlabelticks = list(np.arange(1,37))
        xlabelticks.append('Reference Sequence')
        xlabelticks.append('All Sequences')

        print("Number of sequences without a week:", Deleted)

        plt.ylabel('Percentages')
        plt.xlabel('Week Number')
        plt.title('Percentages by Week')
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        # plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A Percent', 'C Percent', 'G Percent', 'T Percent'))
        plt.savefig("PerWeekPercentages.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(ABar)):
            if ABar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del ABar[i]
                del xticklabels[i]

        ind = np.arange(len(xticklabels))

        ref1 = []
        ref2 = []
        ref3 = []
        ref4 = []
        for i in range(0, len(xticklabels)):
            ref1.append(refPercent[0])

        plt.figure()
        sns.lineplot(x = xticklabels, y = ABar)
        sns.lineplot(x = xticklabels, y = ref1)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage A')
        plt.xlabel('Week Number')
        plt.title('Percentages By Week')
        plt.savefig('AWeekLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(CBar)):
            if CBar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del CBar[i]
                del xticklabels0[i]

        ind = np.arange(len(xticklabels0))

        for i in range(0, len(xticklabels0)):
            ref2.append(refPercent[1])

        plt.figure()
        sns.lineplot(x = xticklabels0, y = CBar)
        sns.lineplot(x = xticklabels0, y = ref2)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage C')
        plt.xlabel('Week Number')
        plt.title('Percentages By Week')
        plt.savefig('CWeekLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(GBar)):
            if GBar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del GBar[i]
                del xticklabels1[i]

        ind = np.arange(len(xticklabels1))

        for i in range(0, len(xticklabels1)):
            ref3.append(refPercent[2])

        plt.figure()
        sns.lineplot(x = xticklabels1, y = GBar)
        sns.lineplot(x = xticklabels1, y = ref3)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage G')
        plt.xlabel('Week Number')
        plt.title('Percentages By Week')
        plt.savefig('GWeekLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        zeros = []
        for i in range(0, len(TBar)):
            if TBar[i] == 0:
                zeros.append(i)

        zeros.reverse()
        
        if len(zeros) != 0:
            for i in zeros:
                del TBar[i]
                del xticklabels2[i]

        ind = np.arange(len(xticklabels2))

        for i in range(0, len(xticklabels2)):
            ref4.append(refPercent[3])

        plt.figure()
        sns.lineplot(x = xticklabels2, y = TBar)
        sns.lineplot(x = xticklabels2, y = ref4)
        plt.xticks(ind, np.array(xlabelticks), rotation=90)
        plt.ylabel('Percentage T')
        plt.xlabel('Week Number')
        plt.title('Percentages By Week')
        plt.savefig('TWeekLine.png', bbox_inches="tight", pad_inches=1)
        plt.clf()

        Bars = [ABar, CBar, GBar, TBar]

        weekDF = pd.DataFrame(Bars, index = ['A', 'C', 'G', 'T'], columns = (xticklabels))
        weekDF.to_csv('WeekValues.csv')
    
    def countryselect(self, metaDataDF, country):
        location = metaDataDF['Location'].values

        countries = []
        for i in location:
            tmp = i.split(' / ')
            tmp[1] = tmp[1].replace(' ', '')
            tmp[0] = tmp[0].replace(' ', '')
            countries.append(tmp[1])

        countries = np.array(countries)
        
        meetselected = list(np.argwhere(countries == country).flatten())
        AccessionIDs = metaDataDF['Accession ID'].values

        goodIDs = []
        for i in meetselected:
            goodIDs.append(AccessionIDs[i])
        
        return goodIDs


        
    def countries(self, allData, refPercent, metaDataDF):
        location = metaDataDF['Location'].values
        countries = []
        continents = []
        for i in location:
            tmp = i.split(' / ')
            tmp[1] = tmp[1].replace(" ", "")
            tmp[0] = tmp[0].replace(" ", "")
            countries.append(tmp[1])
            continents.append(tmp[0])
        
        countries = np.array(countries)
        countries = np.intersect1d(countries, countries)
        continents = np.array(continents)
        continents = np.intersect1d(continents, continents)
        continents = np.array(continents)

        sequences = list(allData.index)

        seqcountries = []
        seqcontinents = []

        for i in sequences:
            tmp = metaDataDF['Location'][i].split('/')
            tmp[1] = tmp[1].replace(" ", "")
            tmp[0] = tmp[0].replace(" ", "")
            seqcountries.append(tmp[1])
            seqcontinents.append(tmp[0])

        contValues = []
        for i in continents:
            contValues.append([[], [], [], []])
        
        countryValues = []
        for i in countries:
            countryValues.append([[], [], [], []])
        
        for i in range(0, len(sequences)):
            APercent = allData['A Percent'][sequences[i]]
            CPercent = allData['C Percent'][sequences[i]]
            GPercent = allData['G Percent'][sequences[i]]
            TPercent = allData['T Percent'][sequences[i]]

            continentLoc = np.argwhere(continents == seqcontinents[i]).flatten()[0]
            countryLoc = np.argwhere(countries == seqcountries[i]).flatten()[0]

            contValues[continentLoc][0].append(APercent)
            contValues[continentLoc][1].append(CPercent)
            contValues[continentLoc][2].append(GPercent)
            contValues[continentLoc][3].append(TPercent)

            countryValues[countryLoc][0].append(APercent)
            countryValues[countryLoc][1].append(CPercent)
            countryValues[countryLoc][2].append(GPercent)
            countryValues[countryLoc][3].append(TPercent)

        Averages = []
        for i in contValues:
            if len(i[0]) != 0:
                AAverage = stat.mean(i[0])
                CAverage = stat.mean(i[1])
                GAverage = stat.mean(i[2])
                TAverage = stat.mean(i[3])
            else:
                AAverage = 0
                CAverage = 0
                GAverage = 0
                TAverage = 0
            Averages.append([AAverage, CAverage, GAverage, TAverage])

        Averages2 = []
        for i in countryValues:
            if len(i[0]) != 0:
                AAverage = stat.mean(i[0])
                CAverage = stat.mean(i[1])
                GAverage = stat.mean(i[2])
                TAverage = stat.mean(i[3])
            else:
                AAverage = 0
                CAverage = 0
                GAverage = 0
                TAverage = 0
            Averages2.append([AAverage, CAverage, GAverage, TAverage])

        WholeAPercent = stat.mean(allData['A Percent'].values)
        WholeCPercent = stat.mean(allData['C Percent'].values)
        WholeGPercent = stat.mean(allData['G Percent'].values)
        WholeTPercent = stat.mean(allData['T Percent'].values)

        WholePercent = [WholeAPercent, WholeCPercent, WholeGPercent, WholeTPercent]
        
        newvals = [WholePercent, refPercent]
        for i in Averages:
            newvals.append(i)

        newvals2 = [WholePercent, refPercent]
        for i in Averages2:
            newvals2.append(i)

        ABar = []
        CBar = []
        GBar = []
        TBar = []

        width = 0.35
        ind = np.arange(len(continents) + 2)

        for i in newvals:
            ABar.append(i[0])
            CBar.append(i[1])
            GBar.append(i[2])
            TBar.append(i[3])

        bars = np.add(ABar, CBar).tolist()
        bars1 = np.add(bars, GBar).tolist()

        p1 = plt.bar(ind, ABar, width)
        p2 = plt.bar(ind, CBar, width, bottom=ABar)
        p3 = plt.bar(ind, GBar, width, bottom=bars)
        p4 = plt.bar(ind, TBar, width, bottom=bars1)

        xticklabels = ['Reference Sequence', 'All Sequences']
        for i in continents:
            xticklabels.append(str(i))

        plt.ylabel('Percentages')
        plt.xlabel('Continent')
        plt.title('Percentages by Continent')
        plt.xticks(ind, np.array(xticklabels), rotation=90)
        # plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A Percent', 'C Percent', 'G Percent', 'T Percent'))
        plt.savefig("PerContinentPercentages.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        Bars = [ABar, CBar, GBar, TBar]

        continentDF = pd.DataFrame(Bars, index = ['A', 'C', 'G', 'T'], columns = (xticklabels))
        continentDF.to_csv('ContinentValues.csv')

        ABar = []
        CBar = []
        GBar = []
        TBar = []

        width = 0.35

        toDrop = []
        for i in range(0, len(newvals2)):
            if newvals2[i][0] == 0:
                toDrop.append(i)
            else:
                ABar.append(newvals2[i][0])
                CBar.append(newvals2[i][1])
                GBar.append(newvals2[i][2])
                TBar.append(newvals2[i][3])

        xticklabels = ['Reference Sequence', 'All Sequences']
        for i in countries:
            xticklabels.append(str(i))

        toDrop.reverse()

        for i in toDrop:
            del xticklabels[i]

        ind = np.arange(len(xticklabels))

        bars = np.add(ABar, CBar).tolist()
        bars1 = np.add(bars, GBar).tolist()

        p1 = plt.bar(ind, ABar, width)
        p2 = plt.bar(ind, CBar, width, bottom=ABar)
        p3 = plt.bar(ind, GBar, width, bottom=bars)
        p4 = plt.bar(ind, TBar, width, bottom=bars1)

        plt.ylabel('Percentages')
        plt.xlabel('Country')
        plt.title('Percentages by Country')
        plt.xticks(ind, np.array(xticklabels), rotation=90, fontsize = 8)
        # plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A Percent', 'C Percent', 'G Percent', 'T Percent'))
        plt.savefig("PerCountryPercentages.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        Bars = [ABar, CBar, GBar, TBar]

        countryDF = pd.DataFrame(Bars, index = ['A', 'C', 'G', 'T'], columns = (xticklabels))
        countryDF.to_csv('CountryValues.csv')

    def monthData(self, allData, refseq):
        monthValues = allData['Month'].values
        moreA = allData['more A'].values
        moreC = allData['more C'].values
        moreG = allData['more G'].values
        moreT = allData['more T'].values

        perMonth = []
        for i in range(0, 12):
            perMonth.append([[], [], [], []])

        for i in range(0, len(monthValues)):
            if monthValues[i] != -1:
                AValues = moreA[i]
                CValues = moreC[i]
                GValues = moreG[i]
                TValues = moreT[i]

                perMonth[monthValues[i] - 1][0].append(AValues)
                perMonth[monthValues[i] - 1][1].append(CValues)
                perMonth[monthValues[i] - 1][2].append(GValues)
                perMonth[monthValues[i] - 1][3].append(TValues)

        AValues = []
        CValues = []
        GValues = []
        TValues = []
        for i in perMonth:
            if len(i[0]) > 1:
                AValues.append(stat.mean(i[0]))
                CValues.append(stat.mean(i[1]))
                GValues.append(stat.mean(i[2]))
                TValues.append(stat.mean(i[3]))
            else:
                AValues.append(0)
                CValues.append(0)
                GValues.append(0)
                TValues.append(0)
        
        final = AValues[-1]
        tmpAValues = AValues
        AValues = [final]
        for i in tmpAValues[:-1]:
            AValues.append(i)

        final = CValues[-1]
        tmpCValues = CValues
        CValues = [final]
        for i in tmpCValues[:-1]:
            CValues.append(i)

        final = GValues[-1]
        tmpGValues = GValues
        GValues = [final]
        for i in tmpGValues[:-1]:
            GValues.append(i)

        final = TValues[-1]
        tmpTValues = TValues
        TValues = [final]
        for i in tmpTValues[:-1]:
            TValues.append(i)

        xticklabels = ['December', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September']
        AValues = AValues[:-2]
        CValues = CValues[:-2]
        GValues = GValues[:-2]
        TValues = TValues[:-2]

        ind = np.arange(0, 10)

        plt.figure()
        sns.lineplot(x = np.arange(0, 10), y = AValues)
        plt.ylabel('Sequence A - Reference A')
        plt.xlabel('Months')
        plt.xticks(ind, np.array(xticklabels), rotation=90)
        plt.savefig("PerMonthADifference.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        plt.figure()
        sns.lineplot(x = np.arange(0, 10), y = CValues)
        plt.ylabel('Sequence C - Reference C')
        plt.xlabel('Months')
        plt.xticks(ind, np.array(xticklabels), rotation=90)
        plt.savefig("PerMonthCDifference.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        plt.figure()
        sns.lineplot(x = np.arange(0, 10), y = GValues)
        plt.ylabel('Sequence G - Reference G')
        plt.xlabel('Months')
        plt.xticks(ind, np.array(xticklabels), rotation=90)
        plt.savefig("PerMonthGDifference.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        plt.figure()
        sns.lineplot(x = np.arange(0, 10), y = TValues)
        plt.ylabel('Sequence T - Reference T')
        plt.xlabel('Months')
        plt.xticks(ind, np.array(xticklabels), rotation=90)
        plt.savefig("PerMonthTDifference.png", bbox_inches="tight", pad_inches=1)
        plt.clf()

        values = [AValues, CValues, GValues, TValues]

        DataFrame = pd.DataFrame(values, index = ['A', 'C', 'G', 'T'], columns = xticklabels)
        DataFrame.to_csv('PerMonthDifference.csv')

    def USASeg(self, MetaDF, outdir):
        # AccessionIDs = MetaDF['Accession ID'].values
        # AccessionIDs = list(np.array(MetaDF.index).flatten())
        location = MetaDF['Location'].values

        states = []
        USPositive = []
        for i in range(0, len(location)):
            tmp = location[i].split(' / ')
            if len(tmp) > 2:
                if tmp[1] == 'USA':
                    tmp[2] = tmp[2].replace(" ", "")
                    states.append(tmp[2])
                    USPositive.append(i)

        goodMetaDF = MetaDF.index[USPositive]

        MetaDF = MetaDF.T[goodMetaDF]
        tmpMeta = MetaDF.T
        AccessionIDs = list(tmpMeta.index)

        days = []

        metaDates = tmpMeta['Collection date'].values
        goodIndices = []
        for i in range(0, len(metaDates)):
            tmp = metaDates[i].split('-')
            for j in tmp:
                j = j.replace(' ', '')
            if len(tmp) > 2:
                goodIndices.append(i)
                time = datetime.datetime(int(tmp[0]), int(tmp[1]), int(tmp[2]))
                day = int(time.strftime('%j'))
                days.append(day)

        AccessionIDs = list(np.array(AccessionIDs)[goodIndices])
        states = list(np.array(states)[goodIndices])

        states = np.array(states)

        liststates = np.intersect1d(states, states)

        PerState = []
        for i in liststates:
            PerState.append([])

        print(len(states), len(AccessionIDs))

        for i in range(0, len(states)):
            PerState[np.argwhere(np.array(liststates) == states[i]).flatten()[0]].append([AccessionIDs[i], days[i]])

        for i in range(0, len(PerState)):
            PerState[i] = sorted(PerState[i], key=itemgetter(1))

        StateDF = pd.DataFrame(PerState, index = liststates)
        StateDF.to_csv(outdir + 'PerState.csv')

        StateDF = StateDF.T
        for i in liststates:
            print(i)
            singleState = StateDF[i].values
            singleStateTimeline = []
            singleStateLengths = []
            for j in singleState:
                if j != None:
                    singleStateTimeline.append(j[1])
                    singleStateLengths.append(1)

            StateValues = []
            for j in range(0, max(singleStateTimeline)+1):
                StateValues.append([])
                
            for j in range(0, len(singleStateTimeline)):
                StateValues[singleStateTimeline[j]].append(singleStateLengths[j])

            forLine = []
            for j in StateValues:
                if len(j) > 0:
                    forLine.append(len(j))
                else:
                    forLine.append(0)
            
            plt.figure()
            sns.lineplot(x = range(0, len(forLine)), y = forLine)
            plt.title(i + ' Cases per Day')
            plt.xlabel('Day')
            plt.ylabel('Case Count')
            plt.savefig(outdir + i + "Cases.png", bbox_inches="tight", pad_inches=1)
            plt.clf()