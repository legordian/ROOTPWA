#include <lorentzFactorFunctionGenerator.h>

#include <TF2.h>
#include <TFile.h>
#include <TTree.h>

#include <reportingUtils.hpp>
#include <primeNumbers.h>
#include <TLSContrib.h>
#include <TJSS.h>
#include <TFhh.h>

using namespace std;
namespace lzf = rpwa::lorentzfactors;


lzf::lorentzFactors* lzf::lorentzFactors::_instance = 0;
const string lzf::lorentzFactors::_lorentzFactorFunctionDirectory = "/home/kbicker/analysis/lorentzFactorCache";
const string lzf::lorentzFactors::_primeNumberFileName = "/opt/rootpwa/primeNumberCache_big.root";


map<lzf::lorentzFactorKey, std::vector<TF2> > lzf::lorentzFactors::getLorentzFactorFunctionsFromRelampl(const lorentzFactorKey& key)
{
	TJSS jss(key.J, key.P, key.s1, key.p1, key.s2, key.p2);
	jss.CalcAmpl();
	const vector<TFhh*>& fhhs = jss.fhh();
	map<lzf::lorentzFactorKey, std::vector<TF2> > retval;
	for(unsigned int i = 0; i < fhhs.size(); ++i) {
		const vector<TLSContrib*>& contribs = fhhs[i]->GetLSt();
		for(unsigned j = 0; j < contribs.size(); ++j) {
			lorentzFactorKey functionKey = key;

/*			//TODO: correct key here
			functionKey.M       = key.M;                // <- where should that come from?
			functionKey.lambda1 = fhhs[i]->GetLambda(); // <- is that right??
			functionKey.lambda2 = fhhs[i]->GetNu();     // <- is that right??
			functionKey.L       = contribs[i]->GetL();
			functionKey.S       = contribs[i]->GetS();
*/
			const string keyName = functionKey.name();
			stringstream sstr;
			sstr << keyName << "_" << retval[functionKey].size();
			TF2 function = convertContribToTF(sstr.str(), *contribs[i]);
			retval[functionKey].push_back(function);
		}
	}
	return retval;
}


TF2 lzf::lorentzFactors::convertContribToTF(const string& functionName, const TLSContrib& contrib)
{
	//TODO: function limits?
	return TF2(functionName.c_str(), "1*x/x*y/y");
}


const vector<TF2>& lzf::lorentzFactors::lorentzFactorFunction(const lorentzFactorKey& key)
{
	if(_lorentzFactorStorage.find(key) == _lorentzFactorStorage.end()) {
		vector<TF2> functions;
		if(not readLorentzFactorFunctionFromFile(key, functions)) {
			std::map<lorentzFactorKey, std::vector<TF2> > newElements = getLorentzFactorFunctionsFromRelampl(key);
			//TODO: handle empty vector<TF2>s
			if(not writeLorentzFactorFunctionsToFiles(newElements)) {
				printErr << "could not write lorentz factor functions to file. Aborting..." << endl;
				throw;
			}
			if(not readLorentzFactorFunctionFromFile(key, functions)) {
				printErr << "could not read lorentz factor functions from file even after writing them. Aborting..." << endl;
				throw;
			}
		}
		_lorentzFactorStorage[key] = functions;
	}
	return _lorentzFactorStorage[key];
}


lzf::lorentzFactors* lzf::lorentzFactors::instance()
{
	if(not _instance) {
		_instance = new lzf::lorentzFactors();
		if(not rpwa::primeNumbers::instance().readCacheFile(_primeNumberFileName)) {
			printErr << "could not initialize prime numbers chach with file '"
			         << _primeNumberFileName << "'. Aborting..." << endl;
			throw;
		}
	}
	return _instance;
}


bool lzf::lorentzFactors::readLorentzFactorFunctionFromFile(const lorentzFactorKey& key,
                                                                  std::vector<TF2>& functions)
{
	if(not functions.empty()) {
		printErr << "got a non-empty vector. Aborting..." << endl;
		throw;
	}
	const string keyName = key.name();
	const string fileName = getFileNameFromKey(key);
	TFile* file = TFile::Open(fileName.c_str(), "READ");
	if(not file) {
		return false;
	}
	TTree* inTree = 0;
	file->GetObject(keyName.c_str(), inTree);
	if(not inTree) {
		printErr << "could not read tree '" << keyName << "' in file '" << fileName << "'. Aborting..." << endl;
		throw;
	}
	if(inTree->GetEntries() < 1) {
		printErr << "tree '" << keyName << "' in file '" << fileName << "' has less than one entry. Aborting..." << endl;
		throw;
	}
	TF2* func = 0;
	if(inTree->SetBranchAddress(keyName.c_str(), &func) != 0) {
		printErr << "could not set address for branch '" << keyName << "' in file '" << fileName << "'. Aborting..." << endl;
		throw;
	}
	functions.resize(inTree->GetEntries());
	for(long i = 0; i < inTree->GetEntries(); ++i) {
		inTree->GetEntry(i);
		functions[i] = TF2(*func);
	}
	file->Close();
	return true;
}


bool lzf::lorentzFactors::writeLorentzFactorFunctionsToFiles(const map<lorentzFactorKey, vector<TF2> >& functions)
{
	for(map<lorentzFactorKey, vector<TF2> >::const_iterator it = functions.begin(); it != functions.end(); ++it)
	{
		const lorentzFactorKey& key  = it->first;
		const vector<TF2>& functions = it->second;
		const string keyName = key.name();
		const string fileName = getFileNameFromKey(key);
		TFile* file = TFile::Open(fileName.c_str(), "NEW");
		if(not file) {
			printWarn << "could not open file '" << fileName << "' for writing. Aborting..." << endl;
			return false;
		}
		TTree* tree = new TTree(keyName.c_str(), keyName.c_str());
		TF2* function = 0;
		tree->Branch(keyName.c_str(), &function);
		for(unsigned int i = 0; i < functions.size(); ++i) {
			function = new TF2(functions[i]);
			tree->Fill();
		}
		file->cd();
		tree->Write();
		file->Close();
	}
	return true;
}


std::string lzf::lorentzFactors::getFileNameFromKey(const lorentzFactorKey& key)
{
	stringstream sstr;
	sstr << _lorentzFactorFunctionDirectory << "/" << key.name() << ".root";
	return sstr.str();
}
