
#include "calcAmplitude.h"

#include <boost/progress.hpp>

#include <TClonesArray.h>
#include <TTree.h>
#include <TTreePerfStats.h>

#include <reportingUtils.hpp>

using namespace std;
using namespace rpwa;


vector<complex<double> >
rpwa::hli::calcAmplitude(const eventMetadata&      eventMeta,
                         const isobarAmplitudePtr& amplitude,
                         const long int            maxNmbEvents,
                         const bool                printProgress,
                         const string&             treePerfStatOutFileName,         // root file name for tree performance result
                         const long int            treeCacheSize,
                         const map<string, pair<double, double> >& otfBin)
{
	vector<complex<double> > retval;

	if(not amplitude) {
		printWarn << "null pointer to isobar decay amplitude. cannot process tree." << endl;
		return retval;
	}
	// initialize amplitude
	amplitude->init();
	const isobarDecayTopologyPtr& decayTopo = amplitude->decayTopology();

	TTree* tree = eventMeta.eventTree();
	if(not tree) {
		printErr << "event tree not found." << endl;
		return retval;
	}

	// create branch pointers and leaf variables
	TBranch*      prodKinMomentaBr  = 0;
	TBranch*      decayKinMomentaBr = 0;
	TClonesArray* prodKinMomenta    = 0;
	TClonesArray* decayKinMomenta   = 0;

	// connect leaf variables to tree branches
	tree->SetBranchAddress(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
	tree->SetBranchAddress(eventMetadata::decayKinematicsMomentaBranchName.c_str(), &decayKinMomenta, &decayKinMomentaBr);
	tree->SetCacheSize(treeCacheSize);
	tree->AddBranchToCache(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  true);
	tree->AddBranchToCache(eventMetadata::decayKinematicsMomentaBranchName.c_str(), true);
	tree->StopCacheLearningPhase();
	TTreePerfStats* treePerfStats = 0;
	if(treePerfStatOutFileName != "") {
		treePerfStats = new TTreePerfStats("ioPerf", tree);
	}

	bool otfBinning = not otfBin.empty();
	vector<double> binningVariables(otfBin.size());
	vector<pair<double, double> > bounds(otfBin.size());
	if(otfBinning) {
		unsigned int otfBinIndex = 0;
		for(map<string, pair<double, double> >::const_iterator elem = otfBin.begin(); elem != otfBin.end(); ++elem) {
			printInfo << "using on-the-fly bin '"
			          << elem->first << ": ["
			          << elem->second.first << ", "
			          << elem->second.second << "]'." << endl;
			int err = tree->SetBranchAddress(elem->first.c_str(), &binningVariables[otfBinIndex]);
			bounds[otfBinIndex] = elem->second;
			++otfBinIndex;
			if(err < 0) {
				printErr << "could not set branch address for branch '" << elem->first << "' (error code " << err << ")." << endl;
				return vector<complex<double> >();
			}
		}
		if(binningVariables.size() != bounds.size()) {
			printErr << "size mismatch between binning variables and bounds ("
			         << binningVariables.size() << "!= " << bounds.size() << ")." << endl;
			return vector<complex<double> >();
		}
	}

	// loop over events
	if(not decayTopo->initKinematicsData(eventMeta.productionKinematicsParticleNames(), eventMeta.decayKinematicsParticleNames())) {
		printWarn << "problems initializing input data. cannot read input data." << endl;
		return retval;
	}
	unsigned long     eventCounter      = 0;
	const long int    nmbEventsTree     = tree->GetEntries();
	const long int    nmbEvents         = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
	                                       : nmbEventsTree);
	boost::progress_display* progressIndicator = (printProgress) ? new boost::progress_display(nmbEvents, cout, "") : 0;
	for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		if(progressIndicator) {
			++(*progressIndicator);
		}

		tree->GetEntry(eventIndex);

		if(otfBinning) {
			bool veto = false;
			for(unsigned int iBinVar = 0; iBinVar < binningVariables.size(); ++iBinVar) {
				if(binningVariables[iBinVar] < bounds[iBinVar].first or binningVariables[iBinVar] >= bounds[iBinVar].second) {
					veto = true;
					break;
				}
			}
			if(veto) {
				retval.push_back(std::complex<double>(0., 0.));
				continue;
			}
		}
		++eventCounter;

		if(not prodKinMomenta or not decayKinMomenta) {
			printWarn << "at least one of the input data arrays is a null pointer: "
			          << "        production kinematics: " << "momenta = " << prodKinMomenta  << endl
			          << "        decay kinematics:      " << "momenta = " << decayKinMomenta << endl
			          << "skipping event." << endl;
			return vector<complex<double> >();
		}

		if(decayTopo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {
			retval.push_back((*amplitude)());
		} else {
			printWarn << "problems reading event[" << eventIndex << "]" << endl;
			return vector<complex<double> >();
		}
	}

	if(printProgress) {
		tree->PrintCacheStats();
	}
	if(treePerfStats) {
		treePerfStats->SaveAs(treePerfStatOutFileName.c_str());
		delete treePerfStats;
	}
	if(otfBinning) {
		printInfo << eventCounter << " events found in on-the-fly bin." << endl;
	}
	return retval;
}
