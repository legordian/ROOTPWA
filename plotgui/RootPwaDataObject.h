///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev:: 836                         $: revision of last commit
// $Author:: schmeing                 $: author of last commit
// $Date:: 2011-12-21 12:31:38 +0100 #$: date of last commit
//
// Description:
//      Header file for the RootPwaDataObject class that provides
//		functionality to read in fit results and integrals from
//		root files generated by rootpwa
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#ifndef RootPwaDataObject_H
#define RootPwaDataObject_H

#include <fstream>
#include <string>
#include <map>
#include <list>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>

#include "fitResult.h"

namespace rpwa{

	class RootPwaDataObject{
	private:
		// Variables
		TFile *_DataFile; // Pointer to the root file containing the data
		TTree *_DataTree; // Pointer to the root tree containing the data
		TBranch *_DataBranch; // Pointer to the root branch containing the data

		std::string _DataFileName;
		std::string _DataTreeName;
		std::string _DataBranchName; ///< Gives the name of the branch which contains the fitResults in the root tree

		std::map<double, fitResult *> _TreeMap;

		static bool _Debug; ///< If set to true, debug messages are printed

		// Functions
		bool FileNotLoaded(); ///< Returns if a File has not been yet loaded and prints an error message if that is the case
		void EmptyTreeMap(); ///< Deletes all elements in TreeMap and clears the map
		void XAxisParameters( double& from, double& to ) const; ///< Determines the start and end value of the x-axis for the histograms

	public:
		// Constructors + Destructors
		RootPwaDataObject(); ///< Constructor
		~RootPwaDataObject(); ///< Destructor

		// Get && Set
		const std::string& DataFileName() const; ///< Returns the filename of the loaded file
		const std::string& DataTreeName() const; ///< Returns the treename of the selected tree
		const std::string& DataBranchName() const; ///< Returns the branchname of the selected branch

		static bool Debug() { return _Debug; } ///< returns debug flag
		static void SetDebug(const bool Debug = true) { _Debug = Debug; } ///< sets debug flag

		// Functions
		TFile *LoadFile( std::string FileName ); ///< Loading the given root tree file
		unsigned int TreesInFile( std::list<std::string>& TreeList ); ///< Reads the names of all trees in the file and returns it in a std::list& given as parameter, returns the number of found trees
		TTree *SelectTree( std::string TreeName ); ///< Selecting root tree in _DataFile by name: Returns Null pointer if tree with given name does not exist
		unsigned int BranchesInTree( std::list<std::string>& BranchList ); ///< Reads the names of all branches in the tree and returns it in a std::list& given as parameter, returns the number of found branches
		TBranch *SelectBranch( std::string BranchName ); ///< Selecting branch in _DataTree by name and returns if selection was successful (if branch exists)

		const std::vector<std::string> *WavesInTree() const;
		unsigned int MapTreeByMassWithHighestLikelihood(); ///< Creates _TreeMap as map of fitResults out of _DataTree sorted by massBinCenter only inserting the fitResults with the Highest logLikelihood for each massBinCenter, returns number of entries in the map
		TH1F *IntensityHist( const std::string& WaveName, const std::string& TitleXAxis ) const; ///< Creates a root histogram of the intensity out of the fitResults in map over the sorting parameter and the spacing between elements has to be constant (The deletion of the histogram is responsibility of calling function)
		TH1F *CoherentSumHist( const std::list<std::string>& WaveNames, const std::string& TitleXAxis ) const; ///< Creates a root histogram of the coherent sum out of the fitResults in map over the sorting parameter and the spacing between elements has to be constant (The deletion of the histogram is responsibility of calling function)
		TH1F *PhaseShiftHist( const std::string& WaveNameA, const std::string& WaveNameB, const std::string& TitleXAxis ) const; ///< Creates a root histogram of the phase shift between wave A and B out of the fitResults in map over the sorting parameter and the spacing between elements has to be constant (The deletion of the histogram is responsibility of calling function)
		void Clear(); ///< Clears the data object without calling destructor of the _DataObject
		void Empty(); ///< Clears the data object and calls destructor of the _DataObject
		std::ostream& Print( std::ostream& Out ) const; ///< Prints all important variables of class
	};

	inline std::ostream& operator <<( std::ostream& Out, const RootPwaDataObject& DataObject ){
		return DataObject.Print(Out);
	}

} // namespace rpwa

#endif /* RootPwaDataObject_H */
