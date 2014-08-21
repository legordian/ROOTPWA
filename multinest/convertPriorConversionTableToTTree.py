#!/usr/bin/env python

import argparse
import os
import sys

import ROOT

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="convert ASCII integral table to ROOT")
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="input file to convert")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="output root file")
	parser.add_argument("--mathematica-list", type=str, metavar="mathematicaList", dest="mathematicaList",
	                    help="produce a file with a list of mathematica coordinates for ListPlot3D[]")
	args = parser.parse_args()
	sys.argv = [sys.argv[0], '-b'] # set root to batch mode

	inputFileName = os.path.abspath(args.inputFile)
	if not os.path.isfile(inputFileName):
		print("ERROR: input file not found. Aborting...")
		sys.exit(1)

	outputFileName = os.path.abspath(args.outputFile)
	if os.path.exists(outputFileName):
		print("ERROR: output file already exists. Aborting...")
		sys.exit(1)

	inputTable = ""
	with open(inputFileName, 'r') as inputFile:
		inputTable = inputFile.read()

	inputLines = inputTable.split('\n')
	nEvents = int(inputLines[0])
	u = int(inputLines[1])
	name = str(nEvents) + "_" + str(u)

	newInputLines = []
	for line in inputLines[3:]:
		if line == "":
			continue
		(r, v, x) = line.split(" ")
		r = float(r)
		v = float(v)
		x = float(x)
		if not newInputLines or r != newInputLines[-1][-1][0]:
			newInputLines.append([])
		newInputLines[-1].append((r, v, x))
	inputLines = newInputLines

	newInputLines = []
	for rLine in inputLines:
		newRLine = []
		for i in range(len(rLine)-1):
			line = rLine[i]
			v = line[1]
			nextV = rLine[i+1][1]
			if v == 0. and nextV == 0.:
				continue
			newRLine.append(line)
		newRLine.append(rLine[-1])
		newInputLines.append(newRLine)
	inputLines = newInputLines

	newInputLines = []
	for rLine in inputLines:
		newRLine = rLine
		for line in reversed(rLine[:-1]):
			newRLine.append((line[0], 1.-line[1], -line[2]))
		newInputLines.append(newRLine)
	inputLines = newInputLines

	cleanedTable = ""
	for rLine in inputLines:
		for line in rLine:
			cleanedTable += str(line[0]) + " " + str(line[1]) + " " + str(line[2]) + "\n"

	rootStream = ROOT.istringstream(cleanedTable)

	outputFile = ROOT.TFile(outputFileName, "NEW")

	tree = ROOT.TTree(name, name)
	tree.ReadStream(rootStream, "r/D:u:x")

	if args.mathematicaList:
		import numpy
		r = numpy.zeros(1, dtype=float)
		v = numpy.zeros(1, dtype=float)
		x = numpy.zeros(1, dtype=float)
		tree.SetBranchAddress("r", r)
		tree.SetBranchAddress("u", v)
		tree.SetBranchAddress("x", x)
		with open(args.mathematicaList, 'w') as mathematicaFile:
			mathematicaFile.write("{")
			for i in range(tree.GetEntries()):
				tree.GetEntry(i)
				line = "{" + str(r[0]) + ", " + str(x[0]) + ", " + str(v[0]) + "}"
				line = line.replace("e", "*^")
				mathematicaFile.write(line)
				if i < tree.GetEntries() - 1:
					mathematicaFile.write(",")
				mathematicaFile.write("\n")
			mathematicaFile.write("}")

	tree.Write()
	outputFile.Close()
