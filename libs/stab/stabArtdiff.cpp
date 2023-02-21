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

void stab_t::stabSetupArtdiff(){

  // Compute viscosity scaling i.e. alpha * mesh_length/ mesh.N
  // so that visc = vmax* (Viscosity Scaling * Viscosity Ramp) 
  viscScale.malloc(mesh.Nelements); 
  for(dlong e=0; e < mesh.Nelements; e++){
    switch (mesh.elementType) {
      case Mesh::TRIANGLES:
        viscScale[e] = ElementViscosityScaleTri2D(e);
        break;
      case Mesh::QUADRILATERALS:
        viscScale[e] = ElementViscosityScaleQuad2D(e);
        break;
      case Mesh::TETRAHEDRA:
        viscScale[e] = ElementViscosityScaleTet3D(e);
        break;
      case Mesh::HEXAHEDRA:
        viscScale[e] = ElementViscosityScaleHex3D(e);
        break;
    }
  } 
  o_viscScale = platform.malloc<dfloat>(viscScale); 



  // Smooth Viscosity Using Vertex Values
  vertexVisc.malloc(mesh.Nelements*mesh.Nverts*dNfields, 0.0); 
  o_vertexVisc = platform.malloc<dfloat>(vertexVisc); 

  // Allocate Memory for Artificial Viscosity
  visc.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*dNfields, 0.0); 
  o_visc = platform.malloc<dfloat>(visc); 

   memory<dfloat> V, invV1, r1, s1, t1;
  // Compute projection matrix 
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      mesh.NodesTri2D(1, r1, s1);
      mesh.VandermondeTri2D(1, r1, s1, invV1);
      mesh.VandermondeTri2D(1, mesh.r, mesh.s, V);
      break; 
    case Mesh::QUADRILATERALS:
      mesh.NodesQuad2D(1, r1, s1);
      mesh.VandermondeQuad2D(1, r1, s1, invV1);
      break; 
    case Mesh::TETRAHEDRA:
      mesh.NodesTet3D(1, r1, s1, t1);
      mesh.VandermondeTet3D(1, r1, s1, t1, invV1);
      break; 
    case Mesh::HEXAHEDRA:
      mesh.NodesHex3D(1,r1, s1, t1); 
      mesh.VandermondeHex3D(1, mesh.r, mesh.s, mesh.t, invV1);
      break; 
  }

  // invert N=1 Vandermonde Matrix
  linAlg_t::matrixInverse(mesh.Nverts, invV1);

  projectVisc.malloc(mesh.Nverts*mesh.Np, 0.0); 
  // Transponse of Projection Operator
  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Nverts; j++){
      dfloat sum = 0; 
      for(int m=0; m<mesh.Nverts; m++){
          sum += V[i*mesh.Nverts + m]*invV1[m*mesh.Nverts + j]; 
      }

      // printf("%.5f ", sum);
      projectVisc[j*mesh.Np + i] = sum; 
    }
    // printf("\n");
  } 
   o_projectVisc = platform.malloc<dfloat>(projectVisc);

     
  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int sNblockV = std::max(1, blockMax/mesh.Np);
  props["defines/" "p_sNblockV"]= sNblockV;

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
  fileName      = oklFilePrefix + "artDiff" + oklFileSuffix;
  kernelName    = "computeViscosity";
  computeViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

  kernelName    = "projectViscosity" + suffix;
  projectViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

}


// void stab_t::stabilizerApplyArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

// // Detect Element ->o_elist!!!!
// detectorApply(o_Q, o_RHS, T); 

// const dfloat alpha = 9.0; 
// computeViscosityKernel(mesh.Nelements, 
//                        alpha, 
//                        o_viscRamp,
//                        o_viscScale,
//                        o_visc); 

// // // Filter solution @ o_eList
// // filterKernel(mesh.Nelements, 
// //              o_eList, 
// //              o_filterM, 
// //              o_Q); 
// }


// Compute h/N for all element types
dfloat stab_t::ElementViscosityScaleTri2D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();
  for(int f=0;f<mesh.Nfaces;++f){
    dlong sid = mesh.Nsgeo*(mesh.Nfaces*e + f);
    dfloat sJ   = mesh.sgeo[sid + mesh.SJID];
    dfloat invJ = mesh.sgeo[sid + mesh.IJID];

    // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);

    h = std::min(h, hest);
  }
  return h/mesh.N;
}

dfloat stab_t::ElementViscosityScaleQuad2D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();

  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<mesh.Np;n++)
    J += mesh.vgeo[mesh.Nvgeo*mesh.Np*e + n + mesh.Np*mesh.JWID];

  for(int f=0;f<mesh.Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<mesh.Nfp;i++)
      sJ += mesh.sgeo[mesh.Nsgeo*(mesh.Nfaces*mesh.Nfp*e + mesh.Nfp*f + i) + mesh.WSJID];

    // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
    // h = 1/(sJ/J)
    dfloat hest = J/sJ;

    h = std::min(h, hest);
  }
  return h/mesh.N;
}

dfloat stab_t::ElementViscosityScaleTet3D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();
  for(int f=0;f<mesh.Nfaces;++f){
    dlong sid   = mesh.Nsgeo*(mesh.Nfaces*e + f);
    dfloat sJ   = mesh.sgeo[sid + mesh.SJID];
    dfloat invJ = mesh.sgeo[sid + mesh.IJID];
    // sJ = A/2, J = 3*V/4 -> sJ/J = 2*A/3*V = 2*A/3*(A*h/3) = 2/h
    // h = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);

    h = std::min(h, hest);
  }
  return h/mesh.N;
}

dfloat stab_t::ElementViscosityScaleHex3D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();

  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<mesh.Np;n++)
    J += mesh.vgeo[mesh.Nvgeo*mesh.Np*e + n + mesh.Np*mesh.JWID];

  for(int f=0;f<mesh.Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<mesh.Nfp;i++)
      sJ += mesh.sgeo[mesh.Nsgeo*(mesh.Nfaces*mesh.Nfp*e + mesh.Nfp*f + i) + mesh.WSJID];

    // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
    // h = 1/(sJ/J)
    dfloat hest = J/sJ;

    h = std::min(h, hest);
  }
  return h/mesh.N;
}


} //namespace libp
