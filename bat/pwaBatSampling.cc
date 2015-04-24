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
//-----------------------------------------------------------
//
// Description:
//      fitting program for rootpwa
//      minimizes pwaLikelihood function
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <cassert>
#include <time.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <BAT/BCModel.h>
#include <BAT/BCParameter.h>
#include <BAT/BCModelOutput.h>

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"
#include "amplitudeTreeLeaf.h"
#include "partialWaveFitHelper.h"


using namespace std;
using namespace ROOT::Math;
using namespace rpwa;


template<typename complexT>
class batLikelihood : public BCModel {

  public:

	batLikelihood(const pwaLikelihood<complexT>& likeli)
		: _likeli(likeli) { }

	double LogAPrioriProbability(const vector<double>& parameters) { return 0.; }

	double LogLikelihood(const vector<double>& parameters) {
		return -_likeli.DoEval(parameters.data());
	}

  private:

	const pwaLikelihood<complexT>& _likeli;

};


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "performs PWA fit for given mass bin and list of waves" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l # -u # -f fitResut [-d amplitude directory -R -o outfile -s seed -x [startvalue] -N -n normfile"
	     << " -a normfile -A # normalisation events -i # -r rank -t # -m # -C -q -z -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -f file    path to fit result file with start values and wave list" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	     << "        -R         use .root amplitude files (default: false)" << endl
#else
	     << "        -R         use .root amplitude files [not supported; ROOT version too low]" << endl
#endif
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -s #       seed for random start values (default: 1234567)" << endl
	     << "        -x #       use fixed instead of random start values (default: 0.01)" << endl

	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to (default: use number of events from acceptance integral file)" << endl
	     << "        -i #       number of iterations (default: 100)" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -t #       relative parameter tolerance (default: 0.0001)" << endl
	     << "        -m #       absolute likelihood tolerance (default: 0.000001)" << endl
	     << "        -C         use half-Cauchy priors (default: false)" << endl
	     << "        -q         run quietly (default: false)" << endl
	     << "        -z         save space by not saving integral and covariance matrices (default: false)" << endl
	     << "        -h         print help" << endl
	     << endl;
	exit(errCode);
}


int
main(int    argc,
     char** argv)
{
	printCompilerInfo();
	printLibraryInfo ();
	printGitHash     ();
	cout << endl;

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");

	// ---------------------------------------------------------------------------
	// internal parameters
	const string       valTreeName           = "pwa";
	const string       valBranchName         = "fitResult_v2";
	double             defaultStartValue     = 0.01;
	bool               useFixedStartValues   = false;
//	const unsigned int maxNmbOfIterations    = 50000;
	int                startValSeed          = 1234567;

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName            = argv[0];
	double       massBinMin          = 0;                      // [MeV/c^2]
	double       massBinMax          = 0;                      // [MeV/c^2]
	string       inputFitResultFileName    = "";                     // wavelist filename
	string       ampDirName          = ".";                    // decay amplitude directory name
	bool         useRootAmps         = false;                  // if true .root amplitude files are read
	//bool         useRootAmps         = true;                   // if true .root amplitude files are read
	string       outFileName         = "fitresult.root";       // output filename
	bool         useNormalizedAmps   = false;                  // if true normalized amplitudes are used
	string       normIntFileName     = "";                     // file with normalization integrals
	string       accIntFileName      = "";                     // file with acceptance integrals
	unsigned int numbAccEvents       = 0;                      // number of events used for acceptance integrals
	unsigned int rank                = 1;                      // rank of fit
	double       minimizerTolerance  = 1e-4;                   // minimizer tolerance
	double       likelihoodTolerance = 1e-6;                   // tolerance of likelihood function
	bool         cauchy              = false;
	bool         quiet               = false;
	bool         saveSpace           = false;
	unsigned int iterations          = 100;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:f:d:Ro:s:x::Nn:a:A:i:r:t:m:Cqzh")) != -1)
		switch (c) {
		case 'l':
			massBinMin = atof(optarg);
			break;
		case 'u':
			massBinMax = atof(optarg);
			break;
		case 'f':
			inputFitResultFileName = optarg;
			break;
		case 'd':
			ampDirName = optarg;
			break;
		case 'R':
#ifdef USE_STD_COMPLEX_TREE_LEAFS
			useRootAmps = true;
#endif
			break;
		case 'o':
			outFileName = optarg;
			break;
		case 's':
			startValSeed = atoi(optarg);
			break;
		case 'x':
			if (optarg)
				defaultStartValue = atof(optarg);
			useFixedStartValues = true;
			break;
		case 'N':
			useNormalizedAmps = true;
			break;
		case 'n':
			normIntFileName = optarg;
			break;
		case 'a':
			accIntFileName = optarg;
			break;
		case 'A':
			numbAccEvents = atoi(optarg);
			break;
		case 'i':
			iterations = atoi(optarg);
			break;
		case 'r':
			rank = atoi(optarg);
			break;
		case 't':
			minimizerTolerance = atof(optarg);
			break;
		case 'm':
			likelihoodTolerance = atof(optarg);
			break;
		case 'C':
			cauchy = true;
			break;
		case 'q':
			quiet = true;
			break;
		case 'z':
			saveSpace = true;
			break;
		case 'h':
			usage(progName);
			break;
		}
	if (normIntFileName.length() <= 1) {
		normIntFileName = "norm.int";
		printWarn << "using default normalization integral file '" << normIntFileName << "'" << endl;
	}
	if (accIntFileName.length() <= 1) {
		accIntFileName = "norm.int";
		printWarn << "using default acceptance normalization integral file "
		          << "'" << accIntFileName << "'" << endl;
	}
	if (inputFitResultFileName.length() <= 1) {
		printErr << "no fit result file specified. Aborting..." << endl;
		usage(progName, 1);
	}
	TFile* inFile = TFile::Open(inputFitResultFileName.c_str(), "READ");
	if(not inFile || inFile->IsZombie()) {
		printErr << "could not open input file '" << inputFitResultFileName << "'. Aborting..." << endl;
		return 1;
	}
	TTree* inTree = 0;
	inFile->GetObject(valTreeName.c_str(), inTree);
	if(not inTree) {
		printErr << "could not find result tree '" << valTreeName << "' in input file '" << inputFitResultFileName << "'. Aborting..." << endl;
		return 1;
	}
	fitResult* result = 0;
	inTree->SetBranchAddress(valBranchName.c_str(), &result);

	if(inTree->GetEntries() != 1) {
		printErr << "result tree '" << valTreeName << "' has more than one entry, NOT IMPLEMENTED." << endl;
		return 1;
	}
	inTree->GetEntry(0);
	if(not result->hasHessian()) {
		printErr << "fit result in input file '" << inputFitResultFileName << "' does not seem to contain errors. Aborting..." << endl;
		return 1;
	}
	// create temporary file to store the wavelist
	char tempFileName[] = "XXXXXX";
	close(mkstemp(tempFileName));
	const string waveListFileName(tempFileName);
	ofstream waveListFile(waveListFileName.c_str());
	rpwa::partialWaveFitHelper::extractWaveList(*result, waveListFile);
	waveListFile.close();

	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to wave list file ......................... '" << waveListFileName << "'" << endl
	     << "    path to amplitude directory .................... '" << ampDirName       << "'" << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)      << endl
	     << "    path to output file ............................ '" << outFileName      << "'" << endl
	     << "    seed for random start values ................... "  << startValSeed            << endl;
	if (useFixedStartValues)
		cout << "    using fixed instead of random start values ..... " << defaultStartValue << endl;
	cout << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName  << "'" << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName   << "'" << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents    << endl
	     << "    rank of spin density matrix .................... "  << rank                    << endl
	     << "    relative parameter tolerance.................... "  << minimizerTolerance << endl
	     << "    absolute likelihood tolerance................... "  << likelihoodTolerance << endl
	     << "    using half-Cauchy priors........................ "  << yesNo(cauchy) << endl
	     << "    saving integral and covariance matrices......... "  << yesNo(not saveSpace) << endl
	     << "    quiet .......................................... "  << yesNo(quiet) << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	const double massBinCenter  = (massBinMin + massBinMax) / 2;
	printInfo << "creating and setting up likelihood function" << endl;
	pwaLikelihood<complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
	L.init(rank, massBinCenter, waveListFileName, normIntFileName, accIntFileName,
	       ampDirName, numbAccEvents, useRootAmps);
	if (not quiet)
		cout << L << endl;
	const unsigned int nmbPar  = L.NDim();
	const unsigned int nmbEvts = L.nmbEvents();
	const double sqrtNmbEvts = sqrt((double)nmbEvts);

	if (cauchy)
		L.setPriorType(L.HALF_CAUCHY);

	printInfo << "using prior: ";
	switch(L.priorType())
	{
		case pwaLikelihood<complex<double> >::FLAT:
			cout << "flat" << endl;
			break;
		case pwaLikelihood<complex<double> >::HALF_CAUCHY:
			cout << "half-cauchy" << endl;
			break;
	}
	unsigned int maxParNameLength = 0;
	for (unsigned int i = 0; i < nmbPar; ++i) {
				const string parName = L.parName(i);
				if (parName.length() > maxParNameLength) {
					maxParNameLength = parName.length();
				}
	}

	batLikelihood<complex<double> >* model = new batLikelihood<complex<double> >(L);
	const unsigned int nChains = model->MCMCGetNChains();
	vector<double> startValues(nChains*nmbPar, 0.);
	printInfo << "setting parameters: " << endl;
	for(unsigned int i = 0; i < nmbPar; ++i) {
		const string& parName = L.parName(i);
		for(unsigned int chainIndex = 0; chainIndex < nChains; ++chainIndex) {
			startValues[chainIndex*nmbPar + i] = result->fitParameter(parName);
		}
		double err = result->fitParameterErr(parName);
		cout << "    parameter [" << setw(3) << i << "] "
		     << setw(maxParNameLength) << parName << " = "
		     << setw(12) << maxPrecisionAlign(startValues[i]) << " +- "
		     << setw(12) << maxPrecisionAlign(err);
		double lowerLimit = startValues[i] - 20.*err;
		double upperLimit = startValues[i] + 20.*err;
		if(parName == "V_flat") {
			lowerLimit = -3.*sqrtNmbEvts;
			upperLimit = 3.*sqrtNmbEvts;
		} else {
			if(lowerLimit < -3.*sqrtNmbEvts) {
				lowerLimit = -3.*sqrtNmbEvts;
			}
			if(upperLimit > 3.*sqrtNmbEvts) {
				upperLimit = 3.*sqrtNmbEvts;
			}
		}
		cout << " (limits: [" << lowerLimit << ", " << upperLimit << "])." << endl;
		model->AddParameter(new BCParameter(L.parName(i).c_str(), lowerLimit, upperLimit));
	}
	model->MCMCSetInitialPositions(startValues);
	BCModelOutput* mout = new BCModelOutput(model, outFileName.c_str());
	mout->WriteMarkovChain(true);
	model->MCMCSetFlagPreRun(false);
	model->MCMCSetNIterationsRun(iterations);
	model->MarginalizeAll();
	mout->WriteMarginalizedDistributions();
	mout->Close();
	delete model;
	delete mout;

	return 0;
}
