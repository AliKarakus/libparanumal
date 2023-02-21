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

void stab_t::detectSetupKlockner(){

  // Initialize Required Memeory
  eList.malloc((mesh.Nelements+mesh.totalHaloPairs)*dNfields); 
  o_eList = platform.malloc<dlong>(eList); 

  // Initialize Required Memeory: Remove LaterAK!
  efList.malloc((mesh.Nelements+mesh.totalHaloPairs)*dNfields); 
  o_efList = platform.malloc<dfloat>(efList); 

  qd.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*dNfields); 
  o_qd = platform.malloc<dfloat>(qd); 

  // Compute Art. Diff. Activation function as well
  if(stabType==Stab::ARTDIFF){ // Vertex based to make it continuous
    // viscActivation.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Nverts*dNfields, 0.0); 
    viscActivation.malloc((mesh.Nelements+mesh.totalHaloPairs)*dNfields, 0.0); 
    o_viscActivation = platform.malloc<dfloat>(viscActivation); 
  }

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

  invVT.malloc(mesh.Np*mesh.Np, 0.0); 
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
   props["defines/" "p_sNverts"] = mesh.Nverts;  
   props["defines/" "p_Nq"]= mesh.N+1;

  if(stabType==Stab::SUBCELL){
     props["defines/" "s_DGDG_TYPE"] = int(0); 
     props["defines/" "s_FVFV_TYPE"] = int(1); 
     props["defines/" "s_DGFV_TYPE"] = int(2); 
   }


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
  if(stabType==Stab::ARTDIFF){
    kernelName    = "detectKlocknerDiffusion" + suffix;
  }else{
    kernelName    = "detectKlockner" + suffix;
  }
  detectKernel  = platform.buildKernel(fileName, kernelName, props);

  if(stabType==Stab::SUBCELL){
    fileName      = oklFilePrefix + "subcell" + oklFileSuffix;
    kernelName    = "detectFindNeigh" + suffix;
    findNeighKernel =  platform.buildKernel(fileName, kernelName, props);
  }

  fileName        = oklFilePrefix + "utilities" + oklFileSuffix; 
  kernelName      = "copyFloat";
  copyFloatKernel = platform.buildKernel(fileName, kernelName, props);  

  kernelName      = "copyInt";
  copyIntKernel   = platform.buildKernel(fileName, kernelName, props); 

  kernelName         = "extractField";
  extractFieldKernel = platform.buildKernel(fileName, kernelName, props); 


}


void stab_t::detectApplyKlockner(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

// Directly copy from field for HJS
if(solverType==Stab::HJS){
  o_qd.copyFrom(o_Q); 
}else if(solverType==Stab::CNS){
  // Use Density field for now AK:
  const int field_id = 0; 
  extractFieldKernel(mesh.Nelements,
                     field_id, 
                     o_Q,
                     o_qd); 
}

if(stabType==Stab::ARTDIFF){
  // Detect elements for each fields i.e. 2
  detectKernel(mesh.Nelements, 
               mesh.o_vgeo, 
               o_modeMap, 
               mesh.o_MM, 
               o_invV, 
               o_LSF, 
               o_BLD, 
               o_qd, 
               o_viscActivation,
               o_eList); 
  
  }else if(stabType==Stab::SUBCELL){
  // // Detect elements for each fields i.e. 2
  // detectKernel(mesh.Nelements, 
  //              mesh.o_vgeo, 
  //              o_modeMap, 
  //              mesh.o_MM, 
  //              o_invV, 
  //              o_LSF, 
  //              o_BLD, 
  //              o_qd, 
  //              o_eList); 

  // findNeighKernel(mesh.Nelements, 
  //                 mesh.o_vmapP, 
  //                o_eList); 



  // Detect elements for each fields i.e. 2
  detectKernel(mesh.Nelements, 
               mesh.o_vgeo, 
               o_modeMap, 
               mesh.o_MM, 
               o_invV, 
               o_LSF, 
               o_BLD, 
               o_qd, 
               o_efList); 


  mesh.halo.Exchange(o_efList, dNfields); 


  findNeighKernel(mesh.Nelements, 
                  mesh.o_vmapP, 
                  o_efList); 

 

  copyIntKernel((mesh.Nelements+mesh.totalHaloPairs)*dNfields, 
              o_efList, 
              o_eList); 

}else{
  // // Detect elements for each fields i.e. 2
  // detectKernel(mesh.Nelements, 
  //              mesh.o_vgeo, 
  //              o_modeMap, 
  //              mesh.o_MM, 
  //              o_invV, 
  //              o_LSF, 
  //              o_BLD, 
  //              o_qd, 
  //              o_eList); 
}


// mesh.halo.Exchange(o_efList, dNfields); 




}


void stab_t::ModeInfoKlocknerTri2D(int _N, memory<int>& _modeMap){
  const int _Np = (_N+1)*(_N+2)/2; 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        if(i+j == id){
          _modeMap[n++] = sk; 
        }
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
  _BLD.malloc(_N+1, 0.0);   

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
