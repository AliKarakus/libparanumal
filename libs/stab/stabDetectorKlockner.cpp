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

void stab_t::DetectorSetupHJSKlockner(){

  // Initialize Required Memeory
  eList.malloc(mesh.Nelements*dNfields); 
  o_eList = platform.malloc<dlong>(eList); 

  qd.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*dNfields); 
  o_qd = platform.malloc<dfloat>(qd); 


  // Compute Art. Diff. Activation function as well!!!!
  // if(stabType==Stab::ARTDIFF){
    viscRamp.malloc(mesh.Nelements*dNfields, 0.0); 
    o_viscRamp = platform.malloc<dfloat>(viscRamp); 
  // }

  memory<dfloat> invV, invVT;
  // Compute 2D to 1D mode map 
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      ModeInfoKlocknerTri2D(mesh.N, modeMap); 
      // Compute invV 
      mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, invV);
      break; 
    case Mesh::QUADRILATERALS:
      ModeInfoKlocknerQuad2D(mesh.N, modeMap); 
      mesh.VandermondeQuad2D(mesh.N, mesh.r, mesh.s, invV);
      break; 
    case Mesh::TETRAHEDRA:
      ModeInfoKlocknerTet3D(mesh.N, modeMap); 
      mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, invV);
      break; 
    case Mesh::HEXAHEDRA:
      ModeInfoKlocknerHex3D(mesh.N, modeMap); 
      mesh.VandermondeHex3D(mesh.N, mesh.r, mesh.s, mesh.t, invV);
      break; 
  }

  invVT.malloc(mesh.Np*mesh.Np); 
  linAlg_t::matrixInverse(mesh.Np, invV); 
  linAlg_t::matrixTranspose(mesh.Np, mesh.Np, invV, mesh.Np, invVT, mesh.Np);

   o_modeMap = platform.malloc<int>(modeMap); 
   o_invV    = platform.malloc<dfloat>(invVT); 

   // Klockner 1D mode operations 
   LeastSquaresFitKlockner(mesh.N, LSF); 
   BaseLineDecayKlockner(mesh.N, BLD); 
   o_LSF = platform.malloc<dfloat>(LSF); 
   o_BLD = platform.malloc<dfloat>(BLD); 

   props["defines/" "p_dNfields"]= dNfields;
   props["defines/" "p_sNfields"]= sNfields;
   props["defines/" "p_Nq"]= mesh.N+1;


   // set kernel name suffix
    std::string suffix;
    if(mesh.elementType==Mesh::TRIANGLES)
      suffix = "Tri2D";
    if(mesh.elementType==Mesh::QUADRILATERALS)
      suffix = "Quad2D";
    if(mesh.elementType==Mesh::TETRAHEDRA)
      suffix = "Tet3D";
    if(mesh.elementType==Mesh::HEXAHEDRA)
      suffix = "Hex3D";

  std::string oklFilePrefix = STAB_DIR "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName      = oklFilePrefix + "detect" + suffix + oklFileSuffix;
  kernelName    = "detectKlocknerHJS" + suffix;
  detectKernel  = platform.buildKernel(fileName, kernelName, props);

  fileName        = oklFilePrefix + "utilities" + oklFileSuffix; 
  kernelName      = "copyFloat";
  copyFloatKernel = platform.buildKernel(fileName, kernelName, props);  
}


void stab_t::DetectorApplyHJSKlockner(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

// Directly copy from field for HJS
o_qd.copyFrom(o_Q); 

// Detect elements for each fields i.e. 2
detectKernel(mesh.Nelements, 
             mesh.o_vgeo, 
             o_modeMap, 
             mesh.o_MM, 
             o_invV, 
             o_LSF, 
             o_BLD, 
             o_qd, 
             o_viscRamp, 
             o_eList); 

}




void stab_t::ModeInfoKlocknerTri2D(int _N, memory<int>& _modeMap){
  const int _Np = (_N+1)*(_N+2)/2; 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
  if(i+j == id)
    _modeMap[n++] = sk; 
  sk++;
      }
    }
    sk=0;
  }

}


void stab_t::ModeInfoKlocknerQuad2D(int _N, memory<int>&_modeMap){
  const int _Np = (_N+1)*(_N+1); 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int j=0; j<_Nmodes1D; j++){
      for(int i=0; i<_Nmodes1D; i++){
         if(std::max(i,j) == id){
            _modeMap[n++] = sk; 
          }
        sk++;
         }
       }
    sk=0;
  }
}

void stab_t::ModeInfoKlocknerTet3D(int _N, memory<int>& _modeMap){
  const int _Np = (_N+1)*(_N+2)*(_N+3)/6; 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        for(int k=0; k<_Nmodes1D-i-j; k++){
          if(i+j+k == id){
            _modeMap[n++] = sk;
          } 
          sk++;
        }
      }
    }
    sk=0;
  }
}

void stab_t::ModeInfoKlocknerHex3D(int _N, memory<int>& _modeMap){
 const int _Np = (_N+1)*(_N+1)*(_N+1); 
 const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int j=0; j<_Nmodes1D; j++){
      for(int i=0; i<_Nmodes1D; i++){
        for(int k=0; k<_Nmodes1D; k++){
         if(std::max(std::max(i,j),k) == id){
            _modeMap[n++] = sk; 
          }
        sk++;
         }
        }
      }
    sk=0;
  }
}



// Least squares fit for 1D modes
void stab_t::LeastSquaresFitKlockner(int _N, memory<dfloat>& _LSF){
  memory<dfloat>  tmp(2*_N); 
  _LSF.malloc( _N); 

  for(int n=0; n<_N; n++){
    const dfloat logmode = log10(n+1); 
    tmp[2*n + 0] = logmode; 
    tmp[2*n + 1] = 1.0; 
  }

  linAlg_t::matrixPseudoInverse(_N, 2, tmp);

  for(int n=0; n<_N; n++){
    _LSF[n] = tmp[n];
  }
}


//  baseline decay (squared)
void stab_t::BaseLineDecayKlockner(int _N, memory<dfloat>& _BLD){
  _BLD.malloc(_N+1);   

  dfloat bsum = 0.0; 
  for(int j=1; j<_N+1; j++)
    bsum +=1.0/pow(j, 2*_N); 

  bsum = 1.0/sqrt(bsum); 

  BLD[0] = 0.0; 
  // baseline decay (squared) 
  for(int n=1; n<_N+1; n++){
    const dfloat bdecay = bsum*1.0/(pow(n,_N));
    _BLD[n] = bdecay*bdecay;
  }
}





} //namespace libp
