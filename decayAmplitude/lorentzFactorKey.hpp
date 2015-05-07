#ifndef LORENTZFACTORKEY_HH
#define LORENTZFACTORKEY_HH


#include<sstream>
#include<string>


namespace rpwa {
	namespace lorentzfactors {

		struct lorentzFactorKey {

			int J;
			int P;
			int L;
			int S;
			int s1;
			int s2;
			int p1;
			int p2;
			int lambda1;
			int lambda2;

			lorentzFactorKey(const int& J_,
			                 const int& P_,
			                 const int& L_,
			                 const int& S_,
			                 const int& s1_,
			                 const int& s2_,
			                 const int& p1_,
			                 const int& p2_,
			                 const int& lambda1_,
			                 const int& lambda2_)
				: J(J_),
				  P(P_),
				  L(L_),
				  S(S_),
				  s1(s1_),
				  s2(s2_),
				  p1(p1_),
				  p2(p2_),
				  lambda1(lambda1_),
				  lambda2(lambda2_) { }

			std::string name() const;

		};

		inline
		bool operator<(const lorentzFactorKey& lhs, const lorentzFactorKey& rhs)
		{
			const unsigned int nVars = 10;
			const int* lhsVars[nVars] = { &lhs.J, &lhs.P, &lhs.L, &lhs.S, &lhs.s1,
			                              &lhs.s2, &lhs.p1, &lhs.p2, &lhs.lambda1, &lhs.lambda2 };
			const int* rhsVars[nVars] = { &rhs.J, &rhs.P, &rhs.L, &rhs.S, &rhs.s1,
			                              &rhs.s2, &rhs.p1, &rhs.p2, &rhs.lambda1, &rhs.lambda2 };
			for(unsigned int i = 0; i < nVars; ++i) {
				if(*(lhsVars[i]) != *(rhsVars[i])) {
					return *(lhsVars[i]) < *(rhsVars[i]);
				}
			}
			return false;
		}

		inline
		std::string lorentzFactorKey::name() const
		{
			std::stringstream sstr;
			sstr << "lorentzFactorKey_"
			     << "J" << J << "_"
			     << "P" << P << "_"
			     << "L" << L << "_"
			     << "S" << S << "_"
			     << "sA" << s1 << "_"
			     << "sB" << s2 << "_"
			     << "pA" << p1 << "_"
			     << "pB" << p2 << "_"
			     << "lbA" << lambda1 << "_"
			     << "lbB" << lambda2;
			return sstr.str();
		}

	}
}


#endif
