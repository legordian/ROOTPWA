#ifndef LORENTZFACTORKEY_HH
#define LORENTZFACTORKEY_HH


#include<sstream>
#include<string>


namespace rpwa {
	namespace lorentzfactors {

		struct lorentzFactorKey {

			int J;
			int P;
			int M;
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
			                 const int& M_,
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
				  M(M_),
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
			const unsigned int nVars = 11;
			const int* lhsVars[nVars] = { &lhs.J, &lhs.P, &lhs.M, &lhs.L, &lhs.S, &lhs.s1,
			                              &lhs.s2, &lhs.p1, &lhs.p2, &lhs.lambda1, &lhs.lambda2 };
			const int* rhsVars[nVars] = { &rhs.J, &rhs.P, &rhs.M, &rhs.L, &rhs.S, &rhs.s1,
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
			sstr << "lorentzFactorKey_";
			sstr << J << "_";
			sstr << P << "_";
			sstr << M << "_";
			sstr << L << "_";
			sstr << S << "_";
			sstr << s1 << "_";
			sstr << s2 << "_";
			sstr << p1 << "_";
			sstr << p2 << "_";
			sstr << lambda1 << "_";
			sstr << lambda2;
			return sstr.str();
		}

	}
}


#endif
