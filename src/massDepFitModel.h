//-----------------------------------------------------------
//
// Description:
//    (BW) Component of mass dependent fit
//      
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITMODEL_HH
#define MASSDEPFITMODEL_HH

#include <complex>
#include <vector>

class TF1;

namespace rpwa {


	class massDepFitComponent;


	class massDepFitModel {

	public:

		massDepFitModel() : _numpar(0), _fsmdFunction(NULL), _fsmdFixed(false) {}
		~massDepFitModel() {}

		void add(massDepFitComponent* comp);

		TF1* getFsmdFunction() const { return _fsmdFunction; }
		void setFsmdFunction(TF1* fsmdFunction);

		bool init(const std::vector<std::string>& waveNames,
		          const std::vector<double>& massBinCenters);





    size_t n() const {return _comp.size();}
    unsigned int numPar() const {return _numpar;}
    
    void setPar(const double* par); // set parameters
    void getPar(double* par);       // return parameters 
    unsigned int nFreeFsmdPar() const {return _fsmdFreeParameters.size();}
    double getFreeFsmdPar(unsigned int i) const;
    void getFreeFsmdLimits(unsigned int i, double& lower, double& upper) const;


    const massDepFitComponent* operator[](size_t i) const {return _comp[i];}
    const std::vector<std::pair<size_t, size_t> >&
      getCompChannel(size_t idx) const { return _compChannel[idx]; }


		std::complex<double> productionAmplitude(const size_t idxWave,
		                                         const double mass,
		                                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;
		double intensity(const size_t idxWave,
		                 const double mass,
		                 const size_t idxMass = std::numeric_limits<size_t>::max()) const;
		double phaseAbsolute(const size_t idxWave,
		                     const double mass,
		                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;
		std::complex<double> spinDensityMatrix(const size_t idxWave,
		                                       const size_t jdxWave,
		                                       const double mass,
		                                       const size_t idxMass = std::numeric_limits<size_t>::max()) const;
		double phase(const size_t idxWave,
		             const size_t jdxWave,
		             const double mass,
		             const size_t idxMass = std::numeric_limits<size_t>::max()) const;

		double calcFsmd(const double mass,
		                const size_t idxMass = std::numeric_limits<size_t>::max()) const;

		std::ostream& print(std::ostream& out) const;

	private:

		bool initMapping();
		bool initFsmd(const std::vector<double>& massBinCenters);

		std::vector<std::string> _waveNames;

    std::vector<massDepFitComponent*> _comp;
    unsigned int _numpar;

		TF1* _fsmdFunction;
		bool _fsmdFixed;
		std::vector<unsigned int> _fsmdFreeParameters;
		std::vector<double> _fsmdValues;

    // mapping for wave -> which components with which channel
    // wavelist in same order as given by wavelist
    std::vector<std::vector<std::pair<size_t, size_t> > > _compChannel;    

	};


} // end namespace rpwa

inline std::ostream& operator<< (std::ostream& out, const rpwa::massDepFitModel& fitModel) {
	return fitModel.print(out);
}

#endif // MASSDEPFITMODEL_HH
