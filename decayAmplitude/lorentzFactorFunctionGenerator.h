#ifndef LORENTZFACTORFUNCTIONGENERATOR_HH
#define LORENTZFACTORFUNCTIONGENERATOR_HH

#include <map>
#include <vector>

#include <lorentzFactorKey.hpp>

class TF2;
class TLSContrib;


namespace rpwa {

	namespace lorentzfactors {

		class lorentzFactors {

		  public:

			const std::vector<TF2>& lorentzFactorFunction(const lorentzFactorKey& key);

			static lorentzFactors* instance();

			static void setDebug(const bool& debug = true) { _debug = debug; }
			static bool debug   ()                         { return _debug; }

		  private:

			lorentzFactors()
				: _lorentzFactorStorage() { }

			static TF2  convertContribToTF(const lorentzFactorKey& key,
			                               const TLSContrib& contrib,
			                               const unsigned int& index);
			static bool readLorentzFactorFunctionFromFile(const lorentzFactorKey& key,
			                                              std::map<lorentzFactorKey, std::vector<TF2> >& functions);
			static bool writeLorentzFactorFunctionsToFiles(const std::map<lorentzFactorKey, std::vector<TF2> >& functions);

			static std::map<lorentzFactorKey, std::vector<TF2> > getLorentzFactorFunctionsFromRelampl(const lorentzFactorKey& key);
			static std::string getFileNameFromKey(const lorentzFactorKey& key);

			std::map<lorentzFactorKey, std::vector<TF2> > _lorentzFactorStorage;
			const static std::vector<TF2> _zeroFunction;

			const static std::string _lorentzFactorFunctionDirectory;
			const static std::string _primeNumberFileName;

			static bool _debug;

			static lorentzFactors* _instance;


		};

	}

}

#endif
