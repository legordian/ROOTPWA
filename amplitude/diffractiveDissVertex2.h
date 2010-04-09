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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      class that describes production vertex in diffractive dissociation
//      beam-Reggeon-(X-system) vertex has exactly one incoming beam and
//      one outgoing X particle, which unambiguously defines the Reggeon
//      kinematics
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DIFFRACTIVEDISSVERTEX2_H
#define DIFFRACTIVEDISSVERTEX2_H


#include <boost/shared_ptr.hpp>

#include "interactionVertex2.h"


namespace rpwa {

  class diffractiveDissVertex2 : public interactionVertex2 {

  public:
  
    diffractiveDissVertex2(const particlePtr& beam,
			   const particlePtr& XSystem);  ///< force vertex to have exactly one incoming (beam) and one outgoing particle (X system)
    diffractiveDissVertex2(const diffractiveDissVertex2& vert);
    virtual ~diffractiveDissVertex2();
		
    virtual diffractiveDissVertex2* clone(const bool cloneInParticles  = false,
					  const bool cloneOutParticles = false) const;

    virtual bool addInParticle (const particlePtr&) { return false; }  ///< disabled; only 1 incoming particle (beam) is allowed
    virtual bool addOutParticle(const particlePtr&) { return false; }  ///< disabled; only 1 outgoing particle (X-system) is allowed

    // diffractive dissociation specific accessors
    // inline particlePtr&       beam   ()       { return inParticles ()[0]; }  ///< returns beam particle
    // inline particlePtr&       XSystem()       { return outParticles()[0]; }  ///< returns X particle
    inline const particlePtr& beam   () const { return inParticles ()[0]; }  ///< returns beam particle
    inline const particlePtr& XSystem() const { return outParticles()[0]; }  ///< returns X particle

    virtual std::ostream& print        (std::ostream& out) const;  ///< prints vertex parameters in human-readable form
    virtual std::ostream& dump         (std::ostream& out) const;  ///< prints all vertex data in human-readable form
    virtual std::ostream& printPointers(std::ostream& out) const;  ///< prints particle pointers strored in vertex

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    static bool _debug;  ///< if set to true, debug messages are printed

  };


  typedef boost::shared_ptr<diffractiveDissVertex2> diffractiveDissVertexPtr;


  inline
  diffractiveDissVertexPtr
  createDiffractiveDissVertex(const particlePtr& beam,
			      const particlePtr& XSystem)
  {
    diffractiveDissVertexPtr v(new diffractiveDissVertex2(beam, XSystem));
    return v;
  }


  inline
  std::ostream&
  operator <<(std::ostream&                 out,
	      const diffractiveDissVertex2& vert)
  {
    return vert.print(out);
  }


}  // namespace rpwa


#endif  // DIFFRACTIVEDISSVERTEX2_H
