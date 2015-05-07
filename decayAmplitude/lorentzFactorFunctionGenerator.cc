#include <lorentzFactorFunctionGenerator.h>

#include <iostream>

#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>

#include <reportingUtils.hpp>

#include <primeNumbers.h>
#include <TLSContrib.h>
#include <TJSS.h>
#include <TFhh.h>

using namespace std;
namespace lzf = rpwa::lorentzfactors;


lzf::lorentzFactors* lzf::lorentzFactors::_instance = 0;
const vector<TF2> lzf::lorentzFactors::_zeroFunction = vector<TF2>(1, TF2("zero_func", "0.*x*y"));
bool lzf::lorentzFactors::_debug = true;
const string lzf::lorentzFactors::_lorentzFactorFunctionDirectory = "/home/kbicker/analysis/lorentzFactorCache";
const string lzf::lorentzFactors::_primeNumberFileName = "/opt/rootpwa/primeNumberCache_big.root";


map<lzf::lorentzFactorKey, std::vector<TF2> > lzf::lorentzFactors::getLorentzFactorFunctionsFromRelampl(const lorentzFactorKey& key)
{
	TJSS jss(key.J, key.P, key.s1, key.p1, key.s2, key.p2);
	jss.CalcAmpl();
	const vector<TFhh*>& fhhs = jss.fhh();
	map<lzf::lorentzFactorKey, std::vector<TF2> > retval;
	if(_debug) {
		printDebug << "called with key " << key.name() << endl;
	}
	for(unsigned int i = 0; i < fhhs.size(); ++i) {
		const vector<TLSContrib*>& contribs = fhhs[i]->GetLSt();
		for(unsigned j = 0; j < contribs.size(); ++j) {
			lorentzFactorKey functionKey = key;
			functionKey.lambda1 = fhhs[i]->GetLambda();
			functionKey.lambda2 = fhhs[i]->GetNu();
			functionKey.L       = contribs[j]->GetL();
			functionKey.S       = contribs[j]->GetS();
			TF2 function = convertContribToTF(functionKey, *contribs[j], retval[functionKey].size());
			retval[functionKey].push_back(function);
			if(_debug) {
				printDebug << "added          function for " << functionKey.name() << endl;
			}
			lorentzFactorKey parityMirroredKey = functionKey;
			parityMirroredKey.lambda1 *= -1.;
			parityMirroredKey.lambda2 *= -1.;
			if(functionKey != parityMirroredKey) {
				TF2 parityMirroredFunction = convertContribToTF(parityMirroredKey,
				                                                *contribs[j],
				                                                retval[functionKey].size());
				retval[parityMirroredKey].push_back(parityMirroredFunction);
				if(_debug) {
					printDebug << "added mirrored function for " << parityMirroredKey.name() << endl;
				}
			}
		}
	}
	return retval;
}


TF2 lzf::lorentzFactors::convertContribToTF(const lorentzFactorKey& key,
                                            const TLSContrib& contrib,
                                            const unsigned int& index)
{
	//TODO: function limits?
	return TF2(functionName.c_str(), "1*x/x*y/y");
}


const vector<TF2>& lzf::lorentzFactors::lorentzFactorFunction(const lorentzFactorKey& key)
{
	if(_lorentzFactorStorage.find(key) == _lorentzFactorStorage.end()) {
		map<lorentzFactorKey, vector<TF2> > elementsFromFile;
		if(not readLorentzFactorFunctionFromFile(key, elementsFromFile)) {
			map<lorentzFactorKey, vector<TF2> > newElements = getLorentzFactorFunctionsFromRelampl(key);
			//TODO: handle empty vector<TF2>s
			if(not writeLorentzFactorFunctionsToFiles(newElements)) {
				printErr << "could not write lorentz factor functions to file. Aborting..." << endl;
				throw;
			}
			if(not readLorentzFactorFunctionFromFile(key, elementsFromFile)) {
				printErr << "could not read lorentz factor functions from file even after writing them. Aborting..." << endl;
				throw;
			}
		}
		bool oneFound = false;
		bool oneNotFound = false;
		for(map<lorentzFactorKey, vector<TF2> >::const_iterator it = elementsFromFile.begin(); it != elementsFromFile.end(); ++it)
		{
			if(_lorentzFactorStorage.find(it->first) != _lorentzFactorStorage.end()) {
				oneFound = true;
			} else {
				oneNotFound = true;
			}
		}
		if(oneFound and oneNotFound) {
			printErr << "lorentz factor cache seems to be in an inconsistent state. Aborting..." << endl;
			throw;
		}
		_lorentzFactorStorage.insert(elementsFromFile.begin(), elementsFromFile.end());
		if(_lorentzFactorStorage.find(key) == _lorentzFactorStorage.end()) {
			_lorentzFactorStorage[key] = _zeroFunction;
		}
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
                                                            map<lorentzFactorKey, vector<TF2> >& elementsFromFile)
{
	if(not elementsFromFile.empty()) {
		printErr << "got a non-empty map. Aborting..." << endl;
		throw;
	}
	const string fileName = getFileNameFromKey(key);
	TFile* file = TFile::Open(fileName.c_str(), "READ");
	if(not file) {
		return false;
	}
	const TList* namesInFile = file->GetListOfKeys();
	for(int i = 0; i < namesInFile->GetEntries(); ++i) {
		string nameInFile(((TKey*)namesInFile->At(i))->GetName());
		if(nameInFile.length() <= 4) {
			printErr << "got an invalid key '" << nameInFile << "' from file '" << fileName << "'. Aborting..." << endl;
			throw;
		}
		if(nameInFile.substr( nameInFile.length() - 4 ) == "_key") {
			continue;
		}
		if(_debug) {
			printDebug << "found key '" << nameInFile << "' in file '" << fileName << "'." << endl;
		}
		TTree* inTree = 0;
		file->GetObject(nameInFile.c_str(), inTree);
		if(not inTree) {
			printErr << "could not read tree '" << nameInFile << "' in file '" << fileName << "'. Aborting..." << endl;
			throw;
		}
		if(inTree->GetEntries() < 1) {
			printErr << "tree '" << nameInFile << "' in file '" << fileName << "' has less than one entry. Aborting..." << endl;
			throw;
		}
		string nameOfKeyInFile = nameInFile + "_key";
		lorentzFactorKey* newKey = 0;
		file->GetObject(nameOfKeyInFile.c_str(), newKey);
		if(not newKey) {
			printErr << "could not read key '" << nameOfKeyInFile << "' in file '" << fileName << "'. Aborting..." << endl;
			throw;
		}
		TF2* func = 0;
		if(inTree->SetBranchAddress(nameInFile.c_str(), &func) != 0) {
			printErr << "could not set address for branch '" << nameInFile << "' in file '" << fileName << "'. Aborting..." << endl;
			throw;
		}
		elementsFromFile[*newKey].resize(inTree->GetEntries());
		for(long i = 0; i < inTree->GetEntries(); ++i) {
			inTree->GetEntry(i);
			elementsFromFile[lorentzFactorKey(*newKey)][i] = TF2(*func);
		}
	}
	file->Close();
	return true;
}


bool lzf::lorentzFactors::writeLorentzFactorFunctionsToFiles(const map<lorentzFactorKey, vector<TF2> >& functions)
{
	const lorentzFactorKey& referenceKey = functions.begin()->first;
	const string& fileName = getFileNameFromKey(referenceKey);
	TFile* file = TFile::Open(fileName.c_str(), "NEW");
	if(not file) {
		printWarn << "could not open file '" << fileName << "' for writing. Aborting..." << endl;
		return false;
	}
	for(map<lorentzFactorKey, vector<TF2> >::const_iterator it = functions.begin(); it != functions.end(); ++it)
	{
		const lorentzFactorKey& key  = it->first;
		if(not key.sameCalcAmplCall(referenceKey)) {
			printErr << "trying to write key '" << key << "' which does not belong "
			         << "into the same file as key '" << referenceKey << "'. Aborting..." << endl;
			throw;
		}
		const vector<TF2>& functions = it->second;
		const string keyName = key.name();
		TTree* tree = new TTree(keyName.c_str(), keyName.c_str());
		TF2* function = 0;
		tree->Branch(keyName.c_str(), &function);
		for(unsigned int i = 0; i < functions.size(); ++i) {
			function = new TF2(functions[i]);
			tree->Fill();
		}
		file->cd();
		{
			lorentzFactorKey writeKey = key;
			stringstream sstr;
			sstr << keyName << "_key";
			writeKey.SetName(sstr.str().c_str());
			writeKey.Write();
		}
		tree->Write();
	}
	file->Close();
	return true;
}


std::string lzf::lorentzFactors::getFileNameFromKey(const lorentzFactorKey& key)
{
	stringstream sstr;
	sstr << _lorentzFactorFunctionDirectory << "/"
	     << "lorentzFactorKey_"
	     << "J" << key.J << "_"
	     << "P" << key.P << "_"
	     << "sA" << key.s1 << "_"
	     << "sB" << key.s2 << "_"
	     << "pA" << key.p1 << "_"
	     << "pB" << key.p2 << ".root";
	return sstr.str();
}
