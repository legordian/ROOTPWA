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
#include <time.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TDecompChol.h"

#include "reportingUtilsEnvironment.h"
#include "conversionUtils.hpp"
#include "pwaLikelihood.h"
#include "fitResult.h"

#include "amplitudeTreeLeaf.h"

#include <multinest.h>


using namespace std;
using namespace rpwa;


namespace {

	class transformationTableContainer
	{

	  public:

		transformationTableContainer(const string& filename, const unsigned int& nEvents, const unsigned int& e)
		{
			TFile* integralFile = TFile::Open(filename.c_str(), "READ");
			if(not integralFile) {
				printErr << "could not open integral file '" << filename << "'. Aborting..." << endl;
				throw;
			}
			stringstream treeName;
			treeName << nEvents << "_" << e;
			TTree* integralTree = (TTree*)integralFile->Get(treeName.str().c_str());
			if(not integralTree) {
				printErr << "could not find integral TTree '" << treeName.str()
				         << "' in integral file '" << filename << ".' Aborting..." << endl;
				throw;
			}
			double r;
			double u;
			double x;
			integralTree->SetBranchAddress("r", &r);
			integralTree->SetBranchAddress("u", &u);
			integralTree->SetBranchAddress("x", &x);
			map<double, vector<pair<double, double> > > rMap;
			for(int i = 0; i < integralTree->GetEntries(); ++i) {
				integralTree->GetEntry(i);
				map<double, vector<pair<double, double> > >::const_iterator rMapIt = rMap.find(r);
				if(rMapIt == rMap.end()) {
					_rValues.push_back(r);
				}
				rMap[r].push_back(pair<double, double>(u, x));
			}
			for(unsigned int i = 0; i < _rValues.size(); ++i) {
				_valueVectors.push_back(rMap[_rValues[i]]);
			}
			assert(_rValues.size() == _valueVectors.size());
			integralFile->Close();
		}

		double transformCoordinate(const double& r, const double& u) const
		{
			pair<unsigned int, unsigned int> rBracket = getClosestRIndices(r, _rValues);
			const std::vector<std::pair<double, double> >& firstValueVector = _valueVectors[rBracket.first];
			pair<unsigned int, unsigned int> firstRValueIndices = getClosestUIndices(u, firstValueVector);
			pair<double, double> firstPoint;
			if(firstRValueIndices.first == firstRValueIndices.second) {
				firstPoint = firstValueVector[firstRValueIndices.first];
			} else {
				firstPoint = interpolate(u, firstValueVector[firstRValueIndices.first], firstValueVector[firstRValueIndices.second]);
			}
			if(rBracket.first == rBracket.second) {
				return firstPoint.second;
			}
			const std::vector<std::pair<double, double> >& secondValueVector = _valueVectors[rBracket.second];
			pair<unsigned int, unsigned int> secondRValueIndices = getClosestUIndices(u, secondValueVector);
			pair<double, double> secondPoint;
			if(secondRValueIndices.first == secondRValueIndices.second) {
				secondPoint = secondValueVector[secondRValueIndices.first];
			} else {
				secondPoint = interpolate(u, secondValueVector[secondRValueIndices.first], secondValueVector[secondRValueIndices.second]);
			}
			const double returnValue = interpolate(r,
			                                       pair<double,double>(rBracket.first, firstPoint.second),
			                                       pair<double,double>(rBracket.second, secondPoint.second)
			                                      ).second;
			return returnValue;
		}

	  private:

		static pair<unsigned int, unsigned int> getClosestRIndices(const double& p, const vector<double>& values) {
			unsigned int highIndex = 0;
			unsigned int lowIndex = 0;
			if(p <= values[0]) {
				// do nothing here
			} else if(p > values[values.size()-1]) {
				lowIndex = highIndex = values.size() - 1;
			} else {
				for(; p > values[highIndex]; ++highIndex);
				lowIndex = highIndex - 1;
			}
			return pair<unsigned int, unsigned int>(lowIndex, highIndex);
		}

		static pair<unsigned int, unsigned int> getClosestUIndices(const double& p,
		                                                          const vector<pair<double, double> >& values)
		{
			unsigned int highIndex = 0;
			unsigned int lowIndex = 0;
			if(p <= values[0].first) {
				// do nothing here
			} else if(p > values[values.size()-1].first) {
				lowIndex = highIndex = values.size() - 1;
			} else {
				for(; p > values[highIndex].first; ++highIndex);
				lowIndex = highIndex - 1;
			}
			return pair<unsigned int, unsigned int>(lowIndex, highIndex);
		}

		static pair<double, double> interpolate(const double& x,
		                          const pair<double, double>& firstPoint,
		                          const pair<double, double>& secondPoint) {
			return pair<double, double>(x, firstPoint.second + (x-firstPoint.first) * ((secondPoint.second - firstPoint.second)/(secondPoint.first - firstPoint.first)));
		}

		vector<vector<pair<double, double> > > _valueVectors;
		vector<double> _rValues;

	};

	class multiNestLogLike : public pwaLikelihood<complex<double> >
	{

	  public:

		multiNestLogLike()
		  : pwaLikelihood<complex<double> >(),
		  _transformationTables(),
		  _transformationMatrix(),
		  _uniformPriors(false) { }

		~multiNestLogLike() {
			for(unsigned int i = 0; i < _transformationTables.size(); ++i) {
				delete _transformationTables[i];
			}
			_transformationTables.clear();
		}

		// operations
		void init(const unsigned int rank,
		          const double& massBinCenter,
		          const std::string& waveListFileName,
		          const std::string& normIntFileName,
		          const std::string& accIntFileName,
		          const std::string& ampDirName    = ".",
		          const unsigned int numbAccEvents = 0,
		          const bool         useRootAmps   = false)  ///< prepares all internal data structures
		{
			pwaLikelihood<complex<double> >::init(rank,
			                                      massBinCenter,
			                                      waveListFileName,
			                                      normIntFileName,
			                                      accIntFileName,
			                                      ampDirName,
			                                      numbAccEvents,
			                                      useRootAmps);
			complexMatrix accIntegral (nmbWaves() + 1, nmbWaves() + 1);  // normalization integral over full phase space with acceptance
			{
				complexMatrix normIntegral(nmbWaves() + 1, nmbWaves() + 1);  // normalization integral over full phase space without acceptance
				vector<double> phaseSpaceIntegral;
				pwaLikelihood<complex<double> >::getIntegralMatrices(normIntegral, accIntegral, phaseSpaceIntegral);
			}
			const unsigned int tempDims = accIntegral.nRows();
			TMatrixT<double> realAccIntegral(2*tempDims, 2*tempDims);
			for(unsigned int i = 0; i < accIntegral.nRows(); ++i) {
				for(unsigned int j = 0; j < accIntegral.nCols(); ++j) {
					const complex<double>& elem = accIntegral.get(i, j);
					const double& real = elem.real();
					const double& imag = elem.imag();
					realAccIntegral[2*i  ][2*j  ] =  real;
					realAccIntegral[2*i+1][2*j+1] =  real;
					realAccIntegral[2*i+1][2*j  ] =  imag;
					realAccIntegral[2*i  ][2*j+1] = -imag;
				}
			}
			TDecompChol decomposition(realAccIntegral);
			if(not decomposition.Decompose()) {
				printErr << "cholesky decomposition failed. Aborting..." << endl;
				throw;
			}
			_transformationMatrix.ResizeTo(2*tempDims, 2*tempDims);
			_transformationMatrix = decomposition.GetU();
			_transformationMatrix = _transformationMatrix.Invert();
		}

		void readTransformationTables(string tableDirectory) {
			if(NDim() == 0) {
				printErr << "likelihood not initialized. Aborting..." << endl;
				throw;
			}
			_transformationTables.resize(NDim(), 0);
			for(unsigned int i = 0; i < NDim(); ++i) {
				stringstream strStr;
				strStr << tableDirectory << "/" << i << ".root";
				_transformationTables[i] = new transformationTableContainer(strStr.str(), nmbEvents(), i);
			}
		}

		void convertCoordinates(const unsigned int& nDim, double* coordinates) const
		{
			if(_uniformPriors) {
				static const double stretchFactor = nmbEvents() * 1.4;
				assert(nDim == NDim());
				for(unsigned int i = 0; i < nDim; ++i) {
					coordinates[i] = 2. * (coordinates[i] - 0.5) * stretchFactor;
				}
				return;
			}
			double r = 0.;
			static std::vector<double> newCoordinates(nDim, 0.);
			for(unsigned int i = 0; i < nDim; ++i) {
				const unsigned int e = nDim - i - 1;
				const double newCoordinate = _transformationTables[e]->transformCoordinate(sqrt(r), coordinates[i]);
				r += newCoordinate * newCoordinate;
				newCoordinates[i] = newCoordinate;
			}
			unsigned int parameterIndexI = 0;
			for(int i = 0; i < _transformationMatrix.GetNrows()-1; ++i) {
				if(i == 1) {
					continue;
				}
				unsigned int parameterIndexJ = 0;
				for(int j = 0; j < _transformationMatrix.GetNcols()-1; ++j) {
					if(j == 1) {
						continue;
					}
					coordinates[parameterIndexI] += _transformationMatrix[i][j] * newCoordinates[parameterIndexJ];
					parameterIndexJ++;
				}
				parameterIndexI++;
			}
		}

		void setUniformPriors(bool uniformPriors = true)
		{
			_uniformPriors = uniformPriors;
			if(_uniformPriors) {
				printInfo << "using uniform priors" << endl;
			} else {
				printInfo << "using Poisson priors" << endl;
			}
		}

		double DoEval(const double* par) const
		{
			// Mind the minus!
			return -pwaLikelihood<complex<double> >::DoEval(par);
		}

	  private:

		std::vector<transformationTableContainer*> _transformationTables;
		TMatrixT<double> _transformationMatrix;
		bool _uniformPriors;

	};

}


void LogLike(double* Cube, int& ndim, int& npars, double& lnew, void* context)
{
	assert(ndim == npars);
	const unsigned int nDim = (unsigned int)ndim;
	multiNestLogLike* L = (multiNestLogLike*)context;
	assert(nDim == L->NDim());
	L->convertCoordinates(nDim, Cube);
	lnew = L->DoEval(Cube);
}

void dumper(int &nSamples,
            int &nlive,
            int &nPar,
            double** physLive,
            double** posterior,
            double** paramConstr,
            double &maxLogLike,
            double &logZ,
            double &INSlogZ,
            double &logZerr,
            void* context)
{

}


void
usage(const string& progName,
      const int     errCode = 0)
{

	cerr << "performs PWA fit for given mass bin and list of waves" << endl
	     << endl
	     << "usage:" << endl
	     << progName
	     << " -l # -u # -w wavelist [-d amplitude directory -R -o outfile -S start value file -N -n normfile"
	     << " [-a normfile] -r rank -M minimizer [-m algorithm -g strategy -t #] -q -h]" << endl
	     << "    where:" << endl
	     << "        -l #       lower edge of mass bin [MeV/c^2]" << endl
	     << "        -u #       upper edge of mass bin [MeV/c^2]" << endl
	     << "        -w file    path to wavelist file" << endl
	     << "        -d dir     path to directory with decay amplitude files (default: '.')" << endl
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	     << "        -R         use .root amplitude files (default: false)" << endl
#else
	     << "        -R         use .root amplitude files [not supported; ROOT version too low]" << endl
#endif
	     << "        -o file    path to output file (default: 'fitresult.root')" << endl
	     << "        -N         use normalization of decay amplitudes (default: false)" << endl
	     << "        -n file    path to normalization integral file (default: 'norm.int')" << endl
	     << "        -a file    path to acceptance integral file (default: 'norm.int')" << endl
	     << "        -A #       number of input events to normalize acceptance to" << endl
	     << "        -r #       rank of spin density matrix (default: 1)" << endl
	     << "        -t         transformation table directory" << endl
	     << "        -U         use uniform instead of Poisson priors" << endl
	     << "        -q         run quietly (default: false)" << endl
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

	// ---------------------------------------------------------------------------
	// parse command line options
	const string progName           = argv[0];
	double       massBinMin         = 0;                      // [MeV/c^2]
	double       massBinMax         = 0;                      // [MeV/c^2]
	string       waveListFileName   = "";                     // wavelist filename
	string       ampDirName         = ".";                    // decay amplitude directory name
	bool         useRootAmps        = false;                  // if true .root amplitude files are read
	string       outFileName        = "fitresult.root";       // output filename
	bool         useNormalizedAmps  = false;                  // if true normalized amplitudes are used
	string       normIntFileName    = "";                     // file with normalization integrals
	string       accIntFileName     = "";                     // file with acceptance integrals
	unsigned int numbAccEvents      = 0;                      // number of events used for acceptance integrals
	unsigned int rank               = 1;                      // rank of fit
	string       tableDirectory     = "";                     // directory with transformation tables
	bool         uniformPriors      = false;                  // use uniform instead of poisson priors
	bool         quiet              = false;
	extern char* optarg;
	// extern int optind;
	int c;
	while ((c = getopt(argc, argv, "l:u:w:d:Ro:S:s:x::Nn:a:A:r:M:m:g:t:Ucqh")) != -1)
		switch (c) {
		case 'l':
			massBinMin = atof(optarg);
			break;
		case 'u':
			massBinMax = atof(optarg);
			break;
		case 'w':
			waveListFileName = optarg;
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
		case 'r':
			rank = atoi(optarg);
			break;
		case 'c':
			break;
		case 't':
			tableDirectory = optarg;
			break;
		case 'U':
			uniformPriors = true;
			break;
		case 'q':
			quiet = true;
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
	if (waveListFileName.length() <= 1) {
		printErr << "no wavelist file specified. aborting." << endl;
		usage(progName, 1);
	}
	// report parameters
	printInfo << "running " << progName << " with the following parameters:" << endl;
	cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << endl
	     << "    path to wave list file ......................... '" << waveListFileName << "'"  << endl
	     << "    path to amplitude directory .................... '" << ampDirName       << "'"  << endl
	     << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)       << endl
	     << "    path to output file ............................ '" << outFileName      << "'"  << endl
	     << "    use normalization .............................. "  << yesNo(useNormalizedAmps) << endl
	     << "        path to file with normalization integral ... '" << normIntFileName  << "'"  << endl
	     << "        path to file with acceptance integral ...... '" << accIntFileName   << "'"  << endl
	     << "        number of acceptance norm. events .......... "  << numbAccEvents            << endl
	     << "    rank of spin density matrix .................... "  << rank                     << endl
	     << "    quiet .......................................... "  << yesNo(quiet)             << endl;

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << endl;
	multiNestLogLike L;
	if(quiet) {
		L.setQuiet();
	}
	L.useNormalizedAmps(useNormalizedAmps);
	const double massBinCenter = (massBinMax + massBinMin) / 2.;
	L.init(rank, massBinCenter, waveListFileName, normIntFileName,
	       accIntFileName, ampDirName, numbAccEvents, useRootAmps);
	L.readTransformationTables(tableDirectory);
	if(not uniformPriors) {
		L.readTransformationTables(tableDirectory);
	}
	L.setUniformPriors(uniformPriors);
	if(not quiet) {
		cout << L << endl;
	}
	const unsigned int nmbPar  = L.NDim();
	const unsigned int nLivePoints = 1000 + 50 * nmbPar;

	if(outFileName.size() > 99) {
		printErr << "output file string is too long (remember, we're interacting with FORTRAN, *sigh*). Aborting..." << endl;
		return 1;
	}
	const char* outFileNameRaw = outFileName.c_str();

	printInfo << "Using " << nLivePoints << " live points" << endl;

	int     IS       = 1;            // do Nested Importance Sampling?
	int     mmodal   = 0;            // do mode separation?
	int     ceff     = 0;            // run in constant efficiency mode?
	int     nlive    = nLivePoints;  // number of live points
	double  efr      = 0.8;          // set the required efficiency
	double  tol      = 0.5;          // tol, defines the stopping criteria
	int     ndims    = nmbPar;       // dimensionality (no. of free parameters)
	int     nPar     = nmbPar;       // total no. of parameters including free & derived parameters
	int     nClsPar  = nmbPar;       // no. of parameters to do mode separation on
	int     updInt   = 1000;         // after how many iterations feedback is required & the output files should be updated
	                                 // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	double  Ztol     = -1E90;        // all the modes with logZ < Ztol are ignored
	int     maxModes = 100;          // expected max no. of modes (used only for memory allocation)

	int     pWrap[ndims];            // which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) {
		pWrap[i] = 0;
	}

//	char    root[100] = outFileNameRaw;     // root for output files
	int     seed      = -1;                 // random no. generator seed, if < 0 then take the seed from system clock
	int     fb        = 1;                  // need feedback on standard output?
	int     resume    = 1;                  // resume from a previous job?
	int     outfile   = 0;                  // write output files?
	int     initMPI   = 1;                  // initialize MPI routines?, relevant only if compiling with MPI
	                                        // set it to F if you want your main program to handle MPI initialization
	double  logZero   = -1E90;              // points with loglike < logZero will be ignored by MultiNest
	int     maxiter   = 0;                  // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
	                                        // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	void*   context   = (void*)&L;          // not required by MultiNest, any additional information user wants to pass

	// calling MultiNest
	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, outFileNameRaw, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLike, dumper, context);

	return 0;
}
