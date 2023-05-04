/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "stab.hpp"

namespace libp {

void stab_t::setTypes(const Stab::SolverType sType, 
                      const Stab::DetectorType dType, 
                      const Stab::StabType stType) {

// Set Solver Dependent Parameters
  if (sType==Stab::HJS) {
    solverType = Stab::HJS;
    // Local Hamilton Jacobi Solver with 2 fields i.e. qp and qm;  
    dNfields = 2;        // number of fields to be detected
    sNfields = 2;        // number of fields to be stabilized
  } else if (sType==Stab::CNS) {
    solverType = Stab::CNS;
    dNfields   = 1; 
    sNfields   = mesh.dim==2 ? 4:5; // non-isothermal flow
  } else if (sType==Stab::INS) {
     solverType = Stab::INS;
      // LIBP_FORCE_ABORT("Stabilization solver type: " << sType <<" is not implemented");
  } else {
    // LIBP_FORCE_ABORT("Unknown solver type: " << sType);
  }

// Set Detector Type
if (dType==Stab::KLOCKNER) {
    detectorType = Stab::KLOCKNER;
 }else if (dType==Stab::PERSSON) {
    detectorType = Stab::PERSSON;
 }else {
    LIBP_FORCE_ABORT("Unknown solver type: " << dType);
  }


// Set Stabilization Type
if (stType==Stab::FILTER) {
    stabType = Stab::FILTER;
 }else if (stType==Stab::LIMITER) {
     stabType = Stab::LIMITER;
 }else if (stType==Stab::ARTDIFF) {
     stabType = Stab::ARTDIFF;
 }else if (stType==Stab::SUBCELL) {
     stabType = Stab::SUBCELL;
 }else {
    LIBP_FORCE_ABORT("Unknown solver type: " << stType);
  }


std::cout<<"Solver Type: "<< solverType<<std::endl; 
std::cout<<"Detector Type: "<< detectorType<<std::endl; 
std::cout<<"Stabilization Type: "<< stabType<<std::endl; 

props["defines/" "p_sNfields"]= sNfields;
props["defines/" "p_dNfields"]= dNfields;

}


} //namespace libp
