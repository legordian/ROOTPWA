import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def calcIntegrals(ampFileList, dataFileName, maxNmbEvents=0, weightFileName="", otfBin={}):
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	ampFiles = []
	ampMetas = []
	for waveName in ampFileList:
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		ampFiles.append(ampFile)
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
			return None
		ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		ampMetas.append(ampMeta)
		dataFile = ROOT.TFile.Open(dataFileName, "READ")
		if not dataFile:
			pyRootPwa.utils.printErr("could not open data file '" + dataFileName + "'.")
			return None
		dataMeta = pyRootPwa.core.eventMetadata.readEventFile(dataFile)
		if not ampMeta:
			pyRootPwa.utils.printErr("could not read amplitude file '" + ampFileName + "'.")
	if not integralMatrix.integrate(ampMetas, maxNmbEvents, weightFileName, dataMeta, otfBin):
		pyRootPwa.utils.printErr("could not run integration")
		del ampFiles
		return None
	del ampFiles
	return integralMatrix
