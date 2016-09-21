#!/usr/bin/env python
import sys

def getOptions(configFile):
    module = __import__(configFile, globals(), locals(), [], -1)
    lastalDict={"-D":module.alignExpect,
    "-E":module.maxEG2,
    "-r":module.matchScore,
    "-q":module.misMatch,
    "-p":module.scoreMatrix,
    "-a":module.gapExistCost,
    "-b":module.gapExtendCost,
    "-A":module.insertExistCost,
    "-B":module.insertExtendCost,
    "-c":module.afflineGap,
    "-x":module.maxScoreDropGap,
    "-y":module.maxScoreDropGapless,
    "-z":module.maxScoreDropGapFinal,
    "-d":module.minScoreGap,
    "-e":module.minAlignScore,
    "-m":module.multiplicity,
    "-l":module.minInitMatchLength,
    "-L":module.maxInitMatchLength,
    "-k":module.searchStep,
    "-W":module.minQuerySize,
    "-s":module.queryStrand,
    "-S":module.dnaStrand,
    "-K":module.numOverlapQuery,
    "-C":module.numOverlapExtend,
    "-P":module.numThreads,
    "-i":module.batchSize,
    "-M":module.minDiffAligns,
    "-T":module.alignType,
    "-n":module.alignPerQueryPosition,
    "-R":module.repeatHandling,
    "-u":module.treatLower,
    "-w":module.timeKludge,
    "-G":module.alternateCode,
    "-t":module.scoreToLikelihood,
    "-g":module.gammaVal,
    "-j":module.output,
    "-Q":module.qualityTreatment}
    
    optionString=''
    for k in lastalDict.keys():
        if lastalDict[k]:
            #print k,lastalDict[k],
            optionString+=('''%s %s ''' % (k,lastalDict[k]))
    return optionString


