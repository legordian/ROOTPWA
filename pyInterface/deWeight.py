#!/usr/bin/env python

import argparse
import sys

import numpy

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Deweights ROOTPWA .root file."
	                                 )
	parser.add_argument("inputFileName", type=str, metavar="<inputFile>", help="input RootPwa file")
	parser.add_argument("outputFileName", type=str, metavar="<outputFile>", help="deweighted output RootPwa file")
	parser.add_argument("-f", "--weightFactor", type=float, default=1., metavar="#", help="weight factor (default = 1)")
	parser.add_argument("-s", "--seed", type=int, metavar="#", dest="seed", default=2451083, help="random number seed (default = 2451083)")

	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo

	inputFile = pyRootPwa.ROOT.TFile.Open(args.inputFileName, "READ")
	if not inputFile:
		printErr("error opening input file '" + args.inputFileName + "'. Aborting...")
		sys.exit(1)

	metaData = pyRootPwa.core.eventMetadata.readEventFile(inputFile)
	if metaData == 0:
		printErr("error reading metaData. Input file is not a RootPWA root file.")
	inputTree = metaData.eventTree()
	maxWeight = 0.

	additionalVariableLabels = metaData.additionalSavedVariableLables()
	additionalVariables = [numpy.zeros(1, dtype=float)] * len(additionalVariableLabels)
	weightIndex = additionalVariableLabels.index("weight")
	for i, additionalVariableLabel in enumerate(additionalVariableLabels):
		inputTree.SetBranchAddress(additionalVariableLabel, additionalVariables[i])

	for i in xrange(inputTree.GetEntries()):
		inputTree.GetEntry(i)
		if additionalVariables[weightIndex] > maxWeight:
			maxWeight = additionalVariables[weightIndex]

	printInfo("maxWeight: " + str(maxWeight))
	maxWeight *= args.weightFactor
	printInfo("adjusted maxWeight to " + str(maxWeight))

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if not outputFile:
		printErr("output file '" + args.outputFileName + "' already exists. Aborting...")
		sys.exit(1)

	fileWriter = pyRootPwa.core.eventFileWriter()
	if not fileWriter.initialize(outputFile,
	                             metaData.userString(),
	                             pyRootPwa.core.eventMetadata.REAL,
	                             metaData.productionKinematicsParticleNames(),
	                             metaData.decayKinematicsParticleNames(),
	                             metaData.binningMap(),
	                             additionalVariableLabels):
		printErr("could not initialize fileWriter. Aborting...")
		inputFile.Close()
		sys.exit(1)
	pyRootPwa.ROOT.gRandom.SetSeed(args.seed)
	acceptedEntries = 0
	overallEntries = 0

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	inputTree.SetBranchAddress(pyRootPwa.core.eventMetadata.productionKinematicsMomentaBranchName, prodKinMomenta)
	inputTree.SetBranchAddress(pyRootPwa.core.eventMetadata.decayKinematicsMomentaBranchName, decayKinMomenta)

	for i in xrange(inputTree.GetEntries()):
		inputTree.GetEntry(i)
		normWeight = additionalVariables[weightIndex] / maxWeight
		cut = pyRootPwa.ROOT.gRandom.Rndm()
		if normWeight > cut:
			fileWriter.addEvent(prodKinMomenta, decayKinMomenta, [float(variable) for variable in additionalVariables])
			acceptedEntries += 1
		overallEntries += 1
	fileWriter.finalize()
	inputFile.Close()

	printSucc("efficiency = %.2f" % (acceptedEntries*100./overallEntries) + "%")
	printSucc("# accepted events = " + str(acceptedEntries))
	printSucc("successfully deweighted input file '" + args.inputFileName + "' and saved as '" + args.outputFileName + "'.")
