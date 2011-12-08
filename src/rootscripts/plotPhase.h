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
// File and Version Information:
// $Id$
//
// Description:
//      Draws phase angle of (wave A - wave B) from tree.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef PLOTPHASE_HH
#define PLOTPHASE_HH


#include <string>
#include <vector>

#include "TTree.h"
#include "TMultiGraph.h"

#include "fitResult.h"


// .............................................................................
// signatures with wave indices
TMultiGraph*
plotPhase(const unsigned int nmbTrees,             // number of fitResult trees
          TTree**            trees,                // array of fitResult trees
          const int          waveIndexA,           // index of first wave
          const int          waveIndexB,           // index of second wave
          const bool         saveEps     = false,  // if set, EPS file with name wave ID is created
          const int*         graphColors = NULL,   // array of colors for graph line and marker
          const bool         drawLegend  = true,   // if set legend is drawn
          const std::string& graphTitle  = "",     // name and title of graph (default is wave IDs)
          const char*        drawOption  = "AP",   // draw option for graph
          const std::string& selectExpr  = "",     // TTree::Draw() selection expression
          const std::string& branchName  = "fitResult_v2");  // fitResult branch name


inline
TMultiGraph*
plotPhase(std::vector<TTree*>&    trees,                // array of fitResult trees
          const int               waveIndexA,           // index of first wave
          const int               waveIndexB,           // index of second wave
          const bool              saveEps     = false,  // if set, EPS file with name wave ID is created
          const std::vector<int>& graphColors = std::vector<int>(),  // array of colors for graph line and marker
          const bool              drawLegend  = true,   // if set legend is drawn
          const std::string&      graphTitle  = "",     // name and title of graph (default is wave IDs)
          const char*             drawOption  = "AP",   // draw option for graph
          const std::string&      selectExpr  = "",     // TTree::Draw() selection expression
          const std::string&      branchName  = "fitResult_v2")  // fitResult branch name
{
	return plotPhase(trees.size(), &(*(trees.begin())), waveIndexA, waveIndexB, saveEps,
	                 &(*(graphColors.begin())), drawLegend, graphTitle, drawOption,
	                 selectExpr, branchName);
}


inline
TMultiGraph*
plotPhase(TTree*             tree,                 // fitResult tree
          const int          waveIndexA,           // index of first wave
          const int          waveIndexB,           // index of second wave
          const bool         saveEps    = false,   // if set, EPS file with name wave ID is created
          const int          graphColor = kBlack,  // array of colors for graph line and marker
          const bool         drawLegend = false,   // if set legend is drawn
          const std::string& graphTitle = "",      // name and title of graph (default is wave IDs)
          const char*        drawOption = "AP",    // draw option for graph
          const std::string& selectExpr = "",      // TTree::Draw() selection expression
          const std::string& branchName = "fitResult_v2")  // fitResult branch name
{
	return plotPhase(1, &tree, waveIndexA, waveIndexB, saveEps, &graphColor,
	                 drawLegend, graphTitle, drawOption, selectExpr, branchName);
}


// .............................................................................
// signatures with wave names
TMultiGraph*
plotPhase(const unsigned int nmbTrees,             // number of fitResult trees
          TTree**            trees,                // array of fitResult trees
          const std::string& waveNameA,            // name of first wave
          const std::string& waveNameB,            // name of second wave
          const bool         saveEps     = false,  // if set, EPS file with name wave ID is created
          const int*         graphColors = NULL,   // array of colors for graph line and marker
          const bool         drawLegend  = true,   // if set legend is drawn
          const std::string& graphTitle  = "",     // name and title of graph (default is wave IDs)
          const char*        drawOption  = "AP",   // draw option for graph
          const std::string& selectExpr  = "",     // TTree::Draw() selection expression
          const std::string& branchName  = "fitResult_v2");  // fitResult branch name


inline
TMultiGraph*
plotPhase(std::vector<TTree*>&    trees,                // array of fitResult trees
          const std::string&      waveNameA,            // name of first wave
          const std::string&      waveNameB,            // name of second wave
          const bool              saveEps     = false,  // if set, EPS file with name wave ID is created
          const std::vector<int>& graphColors = std::vector<int>(),  // array of colors for graph line and marker
          const bool              drawLegend  = true,   // if set legend is drawn
          const std::string&      graphTitle  = "",     // name and title of graph (default is wave IDs)
          const char*             drawOption  = "AP",   // draw option for graph
          const std::string&      selectExpr  = "",     // TTree::Draw() selection expression
          const std::string&      branchName  = "fitResult_v2")  // fitResult branch name
{
	return plotPhase(trees.size(), &(*(trees.begin())), waveNameA, waveNameB, saveEps,
	                 &(*(graphColors.begin())), drawLegend, graphTitle, drawOption,
	                 selectExpr, branchName);
}


inline
TMultiGraph*
plotPhase(TTree*             tree,                 // fitResult tree
          const std::string& waveNameA,            // name of first wave
          const std::string& waveNameB,            // name of second wave
          const bool         saveEps    = false,   // if set, EPS file with name wave ID is created
          const int          graphColor = kBlack,  // array of colors for graph line and marker
          const bool         drawLegend = false,   // if set legend is drawn
          const std::string& graphTitle = "",      // name and title of graph (default is wave IDs)
          const char*        drawOption = "AP",    // draw option for graph
          const std::string& selectExpr = "",      // TTree::Draw() selection expression
          const std::string& branchName = "fitResult_v2")  // fitResult branch name
{
	return plotPhase(1, &tree, waveNameA, waveNameB, saveEps, &graphColor,
	                 drawLegend, graphTitle, drawOption, selectExpr, branchName);
}


#endif  // PLOTPHASE_HH
