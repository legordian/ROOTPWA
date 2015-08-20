#include "calcAmplitude_py.h"

#include <boost/python.hpp>

#include "calcAmplitude.h"

namespace bp = boost::python;


namespace {

	bp::list calcAmplitude(rpwa::eventMetadata&            eventMeta,
	                       const rpwa::isobarAmplitudePtr& amplitude,
	                       const long int                  maxNmbEvents,
	                       const bool                      printProgress,
	                       const std::string&              treePerfStatOutFileName,
	                       const long int                  treeCacheSize,
	                       const bp::dict&                 pyOtfBin)
	{
		std::map<std::string, std::pair<double, double> > otfBin;
		if(bp::len(pyOtfBin.keys()) > 0) {
			bp::list keys = pyOtfBin.keys();
			for(unsigned int i = 0; i < len(keys); ++i) {
				otfBin[bp::extract<std::string>(keys[i])] =
					std::pair<double, double>(bp::extract<double>(pyOtfBin[pyOtfBin.keys()[i]][0]),
					                          bp::extract<double>(pyOtfBin[pyOtfBin.keys()[i]][1]));
			}
		}
		return bp::list(rpwa::hli::calcAmplitude(eventMeta,
		                                         amplitude,
		                                         maxNmbEvents,
		                                         printProgress,
		                                         treePerfStatOutFileName,
		                                         treeCacheSize,
		                                         otfBin));
	}

}


void rpwa::py::exportCalcAmplitude()
{

	bp::def(
		"calcAmplitude"
		, &::calcAmplitude
		, (bp::arg("eventMeta"),
		   bp::arg("amplitude"),
		   bp::arg("maxNmbEvents") = -1,
		   bp::arg("printProgress") = true,
		   bp::arg("treePerfStatOutFileName") = "",
		   bp::arg("treeCacheSize") = 25000000,
		   bp::arg("otfBin") = bp::dict())
	);

}
