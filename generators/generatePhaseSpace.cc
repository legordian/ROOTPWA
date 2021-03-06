///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

/** @brief Simple partial wave event generator (homogeneous in m)
 */


#include<cstdlib>
#include<iostream>
#include<fstream>
#include<getopt.h>
#include<unistd.h>

#include <boost/progress.hpp>

#include<fileUtils.hpp>
#include<generator.h>
#include<generatorManager.h>
#include<particleDataTable.h>
#include<randomNumberGenerator.h>
#include<reportingUtils.hpp>
#include<reportingUtilsEnvironment.h>


using namespace rpwa;
using namespace std;


void printUsage(char* prog, int errCode = 0)
{
	cerr << "usage:" << endl
	     << prog
	     << " -n # [-a # -m # -M # -B # -s #] -o <file> -p <file> -w <file> -k <path> -i <file> -r <file>" << endl
	     << "    where:" << endl
	     << "        -n #       (max) number of events to generate (default: 100)" << endl
	     << "        -a #       (max) number of attempts to do (default: infinity)" << endl
	     << "        -o <file>  ASCII output file (if not specified, generated automatically)" << endl
	     << "        -p         path to particle data table file (default: ./particleDataTable.txt)" << endl
	     << "        -c         if 1 a comgeant eventfile (.fort.26) is written with same naming as the root file (default 0)" << endl
	     << "        -r <file>  reaction config file" << endl
	     << "        -s #       set seed " << endl
	     << "        -M #       lower boundary of mass range in MeV (overwrites values from config file)" << endl
	     << "        -B #       width of mass bin in MeV" << endl
	     << "        --beamfile <file> path to beam file (overrides values from config file)" << endl
	     << "        --noRandomBeam     read the events from the beamfile sequentially" << endl
	     << "        --randomBlockBeam  like --noRandomBeam but with random starting position" << endl
	     << endl;
	exit(errCode);
}


int main(int argc, char** argv)
{
	printCompilerInfo();
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

	unsigned int nEvents = 100;
	unsigned int maxAttempts = 0;
	string outputEvtFileName = "";
	string outputWhtFileName = "";
	string outputComgeantFileName = "";
	string pdgFileName = "./particleDataTable.txt";
	string reactionFile;
	int seed = 123456;
	bool seedSet = false;
	int massLower = 0;
	int massBinWidth = 0;
	bool overrideMass = false;
	bool writeComgeantOut = false;
	string beamfileNameOverride = "";
	int readBeamfileSequentially = 0;
	int readBeamfileRandomBlock = 0;

	static struct option longOptions[] =
	    {
	         { "beamfile", required_argument, 0, 10000 },
	         { "noRandomBeam", no_argument, &readBeamfileSequentially, 1 },
	         { "randomBlockBeam", no_argument, &readBeamfileRandomBlock, 1 },
	         { 0, 0, 0, 0 }
	    };


	int c;
	int optionIndex;
	while ((c = getopt_long(argc, argv, "n:a:o:p:w:k:i:r:m:s:M:B:hc", longOptions, &optionIndex)) != -1) {
		switch (c) {
			case 0:
				if(longOptions[optionIndex].flag != 0) {
					break;
				}
			case 'n':
				nEvents = atoi(optarg);
				break;
			case 'a':
				maxAttempts = atoi(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				seedSet = true;
				break;
			case 'o':
				outputEvtFileName = optarg;
				break;
			case 'p':
				pdgFileName = optarg;
				break;
			case 'r':
				reactionFile = optarg;
				break;
			case 'c':
				writeComgeantOut = true;
				break;
			case 'M':
				massLower = atoi(optarg);
				overrideMass = true;
				break;
			case 'B':
				massBinWidth = atoi(optarg);
				overrideMass = true;
				break;
			case 10000:
				beamfileNameOverride = optarg;
				break;

			case 'h':
				printUsage(argv[0]);
				break;
			default:
				printUsage(argv[0], 5);
				break;
		}
	}

	if(maxAttempts && (maxAttempts < nEvents)) {
		printWarn << "Maximum attempts is smaller than the number of events. Setting it to infinity." << endl;
		maxAttempts = 0;
	}

	if(overrideMass && not (massLower && massBinWidth)) {
		printErr << "'-M' and '-B' can only be set together. Aborting..." << endl;
		exit(2);
	}

	if(not seedSet) {
		printInfo << "Setting random seed to " << seed << endl;
	}
	randomNumberGenerator::instance()->setSeed(seed);

	particleDataTable::readFile(pdgFileName);
	generatorManager generatorMgr;
	if(beamfileNameOverride != "") {
		generatorMgr.overrideBeamFile(beamfileNameOverride);
	}

	if(not generatorMgr.readReactionFile(reactionFile)) {
		printErr << "could not read reaction file. Aborting..." << endl;
		exit(1);
	}

	if(overrideMass) {
		generatorMgr.overrideMassRange(massLower / 1000., (massLower + massBinWidth) / 1000.);
	}
	if(readBeamfileSequentially == 1) {
		generatorMgr.readBeamfileSequentially();
	}
	if(readBeamfileRandomBlock == 1) {
		generatorMgr.readBeamfileSequentially();
		generatorMgr.randomizeBeamfileStartingPosition();
	}

	if(not generatorMgr.initializeGenerator()) {
		printErr << "could not initialize generator. Aborting..." << endl;
		exit(1);
	}

	if(outputEvtFileName == "") {
		stringstream fileName;
		fileName << massLower << "." << massLower + massBinWidth << ".genbod.evt";
		outputEvtFileName = fileName.str();
	}
	ofstream outputEvtFile(outputEvtFileName.c_str());
	printInfo << "output event file: " << outputEvtFileName << endl;

	ofstream outputComgeantFile;
	if(writeComgeantOut) {
		outputComgeantFileName = changeFileExtension(outputEvtFileName, ".fort.26");
		printInfo << "output comgeant file: " << outputComgeantFileName << endl;
		outputComgeantFile.open(outputComgeantFileName.c_str());
	}

	generatorMgr.print(printInfo);

	boost::progress_display* progressIndicator = new boost::progress_display(nEvents, cout, "");

	unsigned int attempts = 0;
	unsigned int eventsGenerated = 0;
	for(; eventsGenerated < nEvents; ++eventsGenerated) {

		attempts += generatorMgr.event();
		if(maxAttempts && (attempts > maxAttempts)) {
			printWarn << "reached maximum attempts. Aborting..." << endl;
			break;
		}
		const generator& gen = generatorMgr.getGenerator();
		rpwa::particle beam = gen.getGeneratedBeam();
		std::vector<rpwa::particle> finalState = gen.getGeneratedFinalState();
		generator::convertEventToAscii(outputEvtFile, beam, finalState);
		if(writeComgeantOut) {
			rpwa::particle recoil = gen.getGeneratedRecoil();
			TVector3 vertex = gen.getGeneratedVertex();
			generator::convertEventToComgeant(outputComgeantFile, beam, recoil, vertex, finalState, false);
		}
		++(*progressIndicator);

	}

	outputEvtFile.close();
	if(writeComgeantOut) {
		outputComgeantFile.close();
	}

	printSucc << "generated " << eventsGenerated << " events." << endl;
	printInfo << "attempts: " << attempts << endl;
	printInfo << "efficiency: " << setprecision(3) << 100 * ((double)eventsGenerated / attempts) << "%" << endl;

}
