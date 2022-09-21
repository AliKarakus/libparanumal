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

void stab_t::DetectorSetupHJSPersson(){

  // Initialize Required Memeory
  eList.malloc(mesh.Nelements*dNfields); 
  o_eList = platform.malloc<dlong>(eList); 

  // Project to modal space spans N-1
  projectNm1.malloc(mesh.Np*mesh.Np);

  // 
  qd.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*dNfields); 
  o_qd = platform.malloc<dfloat>(qd); 

  // Compute Art. Diff. Activation function as well!!!!
  // if(stabType==Stab::ARTDIFF){
  viscRamp.malloc(mesh.Nelements*dNfields, 0.0); 
  o_viscRamp = platform.malloc<dfloat>(viscRamp); 
  // }


  memory<dfloat> V, invV, tmp, truncModes;
  // Compute 2D to 1D mode map 
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      ModeInfoPerssonTri2D(mesh.N, truncModes); 
      mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, V);
      break; 
    case Mesh::QUADRILATERALS:
      ModeInfoPerssonQuad2D(mesh.N, truncModes); 
      mesh.VandermondeQuad2D(mesh.N, mesh.r, mesh.s, V);
      break; 
    case Mesh::TETRAHEDRA:
      ModeInfoPerssonTet3D(mesh.N, truncModes); 
      mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
      break; 
    case Mesh::HEXAHEDRA:
      ModeInfoPerssonHex3D(mesh.N, truncModes); 
      mesh.VandermondeHex3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
      break; 
  }

  // Copy Vandermonde and invert
  invV.malloc(mesh.Np*mesh.Np); invV.copyFrom(V); // = V.clone(); 
  linAlg_t::matrixInverse(mesh.Np, invV); 

  // Cut-off Highest Modes
  tmp.malloc(mesh.Np*mesh.Np); 

  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Np; j++){
      tmp[i*mesh.Np + j] = truncModes[i]*invV[i*mesh.Np + j]; 
    }
  }

  // Transponse of Projection Operator
  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Np; j++){
      dfloat sum = 0; 
      for(int m=0; m<mesh.Np; m++){
          sum += V[i*mesh.Np +m]*tmp[m*mesh.Np + j]; 
      }
      projectNm1[j*mesh.Np + i] = sum; 
    }
  }

   o_projectNm1    = platform.malloc<dfloat>(projectNm1); 

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
  kernelName    = "detectPerssonHJS" + suffix;
  detectKernel  = platform.buildKernel(fileName, kernelName, props);

  fileName        = oklFilePrefix + "utilities" + oklFileSuffix; 
  kernelName      = "copyFloat";
  copyFloatKernel = platform.buildKernel(fileName, kernelName, props);  
}


void stab_t::DetectorApplyHJSPersson(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

// Directly copy from field for HJS
o_qd.copyFrom(o_Q); 

// Detect elements for each fields i.e. 2
detectKernel(mesh.Nelements, 
             mesh.o_vgeo, 
             mesh.o_MM, 
             o_projectNm1, 
             o_qd,
             o_viscRamp, 
             o_eList); 

}




void stab_t::ModeInfoPerssonTri2D(int _N, memory<dfloat>& _truncModes){
  const int _Np = (_N+1)*(_N+2)/2; 
  const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        if(i+j < (_Nmodes1D-1)){
          _truncModes[sk++] = 1.;
        }else{ 
          _truncModes[sk++] = 0.;
        }
      }
    }
}


void stab_t::ModeInfoPerssonQuad2D(int _N, memory<dfloat>& _truncModes){
  const int _Np = (_N+1)*(_N+1); 
  const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D; j++){
        if(std::max(i,j) < (_Nmodes1D-1)){
          _truncModes[sk++] = 1.;
        }else{ 
          _truncModes[sk++] = 0.;
        }
      }
    }
}



void stab_t::ModeInfoPerssonTet3D(int _N, memory<dfloat>& _truncModes){
  const int _Np = (_N+1)*(_N+2)*(_N+3)/6; 
  const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        for(int k=0; k<_Nmodes1D-i-j; k++){
          if((i+j+k)< (_Nmodes1D-1)){
           _truncModes[sk++] = 1.;
            }else{ 
          _truncModes[sk++] = 0.;
          }
        }
      }
    }
}

void stab_t::ModeInfoPerssonHex3D(int _N, memory<dfloat>& _truncModes){
 const int _Np = (_N+1)*(_N+1)*(_N+1); 
 const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int j=0; j<_Nmodes1D; j++){
      for(int i=0; i<_Nmodes1D; i++){
        for(int k=0; k<_Nmodes1D; k++){
         if(std::max(std::max(i,j),k) < (_Nmodes1D-1)){
            _truncModes[sk++] = 1.;
            }else{ 
            _truncModes[sk++] = 0.;
          }
         }
        }
      }
  }



} //namespace libp
