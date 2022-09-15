/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "lss.hpp"

void lss_t::Setup(platform_t& _platform, mesh_t& _mesh, stab_t _stab, 
                   lssSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  stab = _stab; 
  comm = mesh.comm;
  settings = _settings;


  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  advection    = settings.compareSetting("SOLVER TYPE", "ADVECTION")  ? 1:0; 
  redistance   = settings.compareSetting("SOLVER TYPE", "REDISTANCE") ? 1:0; 

  cubature     = settings.compareSetting("ADVECTION TYPE", "CUBATURE")? 1:0; 
  compressible = settings.compareSetting("FLOW TYPE", "COMPRESSIBLE") ? 1:0; 

  Nfields      = redistance ? 2:1; 

  //setup cubature
  if (cubature) {
    mesh.CubatureSetup();
    mesh.CubaturePhysicalNodes();
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "max"});

  /*setup trace halo exchange */
  qTraceHalo = mesh.HaloTraceSetup(Nfields); 
  vTraceHalo = mesh.HaloTraceSetup(mesh.dim); //velocity field

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  }

  // compute samples of q at interpolation nodes
  q.malloc((Nlocal+Nhalo)*Nfields);
  o_q = platform.malloc<dfloat>(q);

  if(redistance){
    // this holds the gradient of qM, qP
    gradq.malloc(Nlocal*Nfields*mesh.dim);
    o_gradq = platform.malloc<dfloat>(gradq);

    phi.malloc(Nlocal+Nhalo); 
    o_phi = platform.malloc<dfloat>(phi); 

    // this saves the history for qP, and qM for recontruction
    if(settings.compareSetting("TIME RECONSTRUCTION", "ENO2")){
      Nrecon = 4; // Second order reconstruction, 
    }else if(settings.compareSetting("TIME RECONSTRUCTION", "ENO3")){
      Nrecon = 6; // Third order reconstruction
    }

    
    reconstructTime.malloc(Nrecon); 
    o_reconstructTime = platform.malloc<dfloat>(reconstructTime); 

    phiH.malloc(Nlocal*Nfields*Nrecon); 
    o_phiH = platform.malloc<dfloat>(phiH); 

  }else if(advection){
    U.malloc((Nlocal+Nhalo)*mesh.dim);
    o_U = platform.malloc<dfloat>(U);
  }



  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(1); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; /*copy base occa properties*/

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;


  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  if(cubature){
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = std::max(1,blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = std::max(1,blockMax/cubMaxNodes);
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }
  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_Nrecon"] = Nrecon;


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

  std::string oklFilePrefix = DLSS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  if(advection){

    // kernels from volume file
    if(cubature){
      fileName   = oklFilePrefix + "lssAdvectionCubature" + suffix + oklFileSuffix;
      kernelName = "lssAdvectionCubatureVolume" + suffix;
    } else {
      fileName   = oklFilePrefix + "lssAdvection" + suffix + oklFileSuffix;
      kernelName = "lssAdvectionVolume" + suffix;
    }
    advectionVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
   
    // kernels from surface file
    if(cubature){
    fileName   = oklFilePrefix + "lssAdvectionCubature" + suffix + oklFileSuffix;
    kernelName = "lssAdvectionCubatureSurface" + suffix;
    } else {
    fileName   = oklFilePrefix + "lssAdvection" + suffix + oklFileSuffix;
    kernelName = "lssAdvectionSurface" + suffix;
   }
    advectionSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    if(compressible){
      fileName   = oklFilePrefix + "lssCompressible" + suffix + oklFileSuffix;

      kernelName = "lssCompressibleVolume" + suffix;
      compressibleVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "lssCompressibleSurface" + suffix;
      compressibleSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
    }

    if (mesh.dim==2) {
      fileName   = oklFilePrefix + "lssSetFlowField2D" + oklFileSuffix;
      kernelName = "lssSetFlowField2D";
    }else{
      fileName   = oklFilePrefix + "lssSetFlowField3D" + oklFileSuffix;
      kernelName = "lssSetFlowField3D";
    }

    setFlowFieldKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  }else if(redistance){

    fileName   = oklFilePrefix + "lssRedistance" + suffix + oklFileSuffix;
    kernelName = "lssRedistanceVolume" + suffix; 
    redistanceVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo); 

    fileName   = oklFilePrefix + "lssRedistance" + suffix + oklFileSuffix;
    kernelName = "lssRedistanceSurface" + suffix;
    redistanceSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo); 


    
    fileName   = oklFilePrefix + "lssRedistanceUtilites" + oklFileSuffix;
    kernelName = "lssRedistanceSetFields";
    redistanceSetFieldsKernel = platform.buildKernel(fileName, kernelName, kernelInfo); 

    if(settings.compareSetting("TIME RECONSTRUCTION", "ENO2")){ 
    kernelName = "redistanceReconstructENO2"; 
    timeReconstructKernel = platform.buildKernel(fileName, kernelName, kernelInfo); 


    kernelName = "redistanceInitialHistoryENO2"; 
    timeInitialHistoryKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

     }else{

    kernelName = "redistanceReconstructENO3"; 
    timeReconstructKernel = platform.buildKernel(fileName, kernelName, kernelInfo); 

    kernelName = "redistanceInitialHistoryENO3"; 
    timeInitialHistoryKernel = platform.buildKernel(fileName, kernelName, kernelInfo);
    }
  }

  // kernels from initialization files
  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "lssInitialCondition2D" + oklFileSuffix;
    kernelName = "lssInitialCondition2D";
  }else{
    fileName   = oklFilePrefix + "lssInitialCondition3D" + oklFileSuffix;
    kernelName = "lssInitialCondition3D";    
  }  
  initialConditionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

}
