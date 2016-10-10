#!/usr/bin/env python
import sys

def getOptions(configFile):
    module = __import__(configFile, globals(), locals(), [], -1)
    lastdbDict={"-R":module.repeatHandling,
    "-c":module.softMask,
    "-u":module.seedScheme,
    "-w":module.stepSize,
    "-W":module.window,
    "-s":module.volumeSize,
    "-Q":module.inputFormat,
    "-P":module.numThreads,
    "-m":module.seedPattern,
    "-a":module.alphabet,
    "-i":module.matches,
    "-b":module.bucketDepth,
    "-C":module.childTable,
    "-x":module.letterCounting,
    "-v":module.verbose,
    "-V":module.version
    }

    optionString=''
    for k in lastdbDict.keys():
        if lastdbDict[k]:
            #print k,lastalDict[k],
            optionString+=('''%s %s ''' % (k,lastdbDict[k]))
    return optionString
