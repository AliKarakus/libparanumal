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

#include "cns.hpp"

void cns_t::Setup(platform_t& _platform, mesh_t& _mesh, stab_t &_stab, 
                  cnsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  {
  stab = _stab; 
  int Nc = 1; 
  mesh_t meshC = mesh.SetupNewDegree(Nc);
  
  bool verbose = true,  unique  = true; 
  stab.ogs.Setup(mesh.Nelements*mesh.Nverts, meshC.globalIds, mesh.comm, 
                ogs::Signed, ogs::Auto, unique, verbose, platform); 


  stab.weight.malloc(mesh.Nelements*mesh.Nverts, 1.0);
  stab.ogs.GatherScatter(stab.weight, 1, ogs::Add, ogs::Sym); 
  for(int i=0; i < mesh.Nelements*mesh.Nverts; i++ ){ 
    stab.weight[i] = 1./stab.weight[i];
  }

    stab.o_weight = platform.malloc<dfloat>(stab.weight); 

  }


{
  //get physical paramters
  useNonDimensionalEqn = settings.compareSetting("NONDIMENSIONAL EQUATIONS","TRUE")? 1:0;
  
  Nph  = 6;  
  pCoeff.malloc(Nph,0.0);
  MUID = 0; 
  GMID = 1; 
  PRID = 2; 
  RRID = 3; 
  CPID = 4; 
  CVID = 5; 

  settings.getSetting("GAMMA", gamma);
  settings.getSetting("PRANDTL NUMBER", Pr); 

  if(useNonDimensionalEqn){
    settings.getSetting("REYNOLDS NUMBER", Re);
    settings.getSetting("MACH NUMBER", Ma);
    // Set viscosity to 1/Re 
    mu = 1.0/Re; 
    // Set gas constant to 1/(gamma*Ma^2) 
    R  = 1.0/(gamma*Ma*Ma); 
    // Reference Temperature // Tref = p_ref/r_ref * gamma
    // Assuming r_ref = 1, u_ref = 1, p_ref = r_ref * u_ref^2/(Ma^2*gamma)
    Tref = 1;
  }else{
    settings.getSetting("BULK VISCOSITY", mu);
    settings.getSetting("GAS CONSTANT", R); 
    settings.getSetting("REFERENCE TEMPERATURE", Tref); 
  }

  cp=R*gamma/(gamma-1.0); 
  cv=R/(gamma-1.0); 

    pCoeff[MUID] = mu; // Bulk Viscosity
    pCoeff[PRID] = Pr; // Prandtl Number
    pCoeff[RRID] = R;  // Specific Gas Constant
    pCoeff[GMID] = gamma; // 
    pCoeff[CPID] = cp;
    pCoeff[CVID] = cv;

  // // Set viscosity treatment
  // int Nph = 0; 
  if(settings.compareSetting("VISCOSITY TYPE","CONSTANT")){
    viscType = 1; // CONSTANT
  }else if(settings.compareSetting("VISCOSITY TYPE","SUTHERLAND")){
    viscType = 2; // SUTHERLAND temperature dependency
    Nph      += 4;  
    EXID = 6; 
    TRID = 7; 
    TSID = 8; 
    CSID = 9; 
    pCoeff.realloc(Nph);
    // // Tr = 273.15; Ts = 110.4; exp = 1.5;  
    pCoeff[EXID] = 1.5;            // exponent  
    pCoeff[TRID] = 1.0/Tref;       // inverse of reference temperature = 1/Tref   
    pCoeff[TSID] = 110.4/273.15;   // Ts/Tref approximately  
    pCoeff[CSID] = pow(pCoeff[TSID],pCoeff[EXID])*(1.0+pCoeff[TSID])/(2.0*pCoeff[TSID]*pCoeff[TSID]); // exponent  

  }else if(settings.compareSetting("VISCOSITY TYPE","POWER")){
    viscType = 3; // POWER LAW temperature dependency 
  }
}

  cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  isothermal = (settings.compareSetting("ISOTHERMAL", "TRUE")) ? 1:0;
  inviscid   =  mu < 1E-10 ? 1:0;

  //setup cubature
  if (cubature) {
    mesh.CubatureSetup();
    mesh.CubaturePhysicalNodes();
  }

  // Number of conserved fields
  Nfields   = (mesh.dim==3) ? 4:3;
  // include energy equation for non-isothermal model
  if (!isothermal) Nfields++; 

  if(stab.stabType==Stab::ARTDIFF){
    if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
      Ngrads = mesh.dim*Nfields; // for all conservative fields
    }else if(stab.settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
      Ngrads = mesh.dim*mesh.dim; // just for velocity field      
    }
  }else{
    Ngrads = mesh.dim*mesh.dim; // only for velocity
  }


  dlong NlocalFields = mesh.Nelements*mesh.Np*Nfields;
  dlong NhaloFields  = mesh.totalHaloPairs*mesh.Np*Nfields;
  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;

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

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "max", "amx", "set"});

  /*setup trace halo exchange */
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  // compute samples of q at interpolation nodes
  q.malloc(NlocalFields+NhaloFields);
  o_q = platform.malloc<dfloat>(q);

  gradq.malloc(NlocalGrads+NhaloGrads);
  o_gradq = platform.malloc<dfloat>(gradq);

  Vort.malloc(mesh.dim*mesh.Nelements*mesh.Np);
  o_Vort = platform.malloc<dfloat>(Vort);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // viscosity model
  o_pCoeff = platform.malloc<dfloat>(pCoeff); 


  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_Ngrads"]= Ngrads;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (cubature) {
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = std::max(1, blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = std::max(1, blockMax/cubMaxNodes);
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }

{
  kernelInfo["defines/" "p_viscType"]= viscType;

  kernelInfo["defines/" "p_MUID"]= MUID;
  kernelInfo["defines/" "p_GMID"]= GMID;
  kernelInfo["defines/" "p_PRID"]= PRID;
  kernelInfo["defines/" "p_RRID"]= RRID;
  kernelInfo["defines/" "p_CPID"]= CPID;
  kernelInfo["defines/" "p_CVID"]= CVID;
  // if(viscType==2){
  kernelInfo["defines/" "p_EXID"]= EXID;
  kernelInfo["defines/" "p_TRID"]= TRID;
  kernelInfo["defines/" "p_TSID"]= TSID;
  kernelInfo["defines/" "p_CSID"]= CSID;
  // }


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

  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  if (isothermal) {
    if (cubature) {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsIsothermalCubatureVolume" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalCubatureVolume" + suffix;

      cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
      // kernels from surface file
      fileName   = oklFilePrefix + "cnsIsothermalCubatureSurface" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalCubatureSurface" + suffix;

      cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
    } else {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsIsothermalVolume" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalVolume" + suffix;

      volumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
      // kernels from surface file
      fileName   = oklFilePrefix + "cnsIsothermalSurface" + suffix + oklFileSuffix;
      kernelName = "cnsIsothermalSurface" + suffix;

      surfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  } else {
    if (cubature) {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsCubatureVolume" + suffix + oklFileSuffix;

      if(stab.stabType==Stab::ARTDIFF){
        if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
          kernelName = "cnsCubatureVolumeADiffLaplace" + suffix;
        }else if(stab.settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
          kernelName = "cnsCubatureVolumeADiffPhysical" + suffix;
        }
      }else{
        kernelName = "cnsCubatureVolume" + suffix;
      }
      cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                               kernelInfo);


      fileName   = oklFilePrefix + "cnsCubatureSurface" + suffix + oklFileSuffix;
      if(stab.stabType==Stab::ARTDIFF){
        if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
          kernelName = "cnsCubatureSurfaceADiffLaplace" + suffix;

        }else if(stab.settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
          kernelName = "cnsCubatureSurfaceADiffPhysical" + suffix;

        }
      }else{
        kernelName = "cnsCubatureSurface" + suffix;
      }

      cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                               kernelInfo);
    } else {
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsVolume" + suffix + oklFileSuffix;

       if(stab.stabType==Stab::ARTDIFF){
        if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
          kernelName = "cnsVolumeADiffLaplace" + suffix;
        }else if(stab.settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
          kernelName = "cnsVolumeADiffPhysical" + suffix;
        }
        }else{
        kernelName = "cnsVolume" + suffix;
       }

      volumeKernel =  platform.buildKernel(fileName, kernelName,
                                             kernelInfo);


      // kernels from surface file
      fileName   = oklFilePrefix + "cnsSurface" + suffix + oklFileSuffix;
      if(stab.stabType==Stab::ARTDIFF){
        if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
          kernelName = "cnsSurfaceADiffLaplace" + suffix;

        }else if(stab.settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
          kernelName = "cnsSurfaceADiffPhysical" + suffix;

        }
      }else{
        kernelName = "cnsSurface" + suffix;
      }
      surfaceKernel = platform.buildKernel(fileName, kernelName,
                                             kernelInfo);
    }
  }


  // kernels from volume file Add isothermal version as well AK. 
  fileName   = oklFilePrefix + "cnsGradVolume" + suffix + oklFileSuffix;
 
  if(stab.stabType==Stab::ARTDIFF){
    if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
      kernelName = "cnsGradVolumeConservative" + suffix; // gradient of all conservative fields
    }else{
      kernelName = "cnsGradVolume" + suffix; // Compute gradient of velocity field
    }
  }else{
    kernelName = "cnsGradVolume" + suffix; 
  }

  gradVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  // kernels from surface file
  fileName   = oklFilePrefix + "cnsGradSurface" + suffix + oklFileSuffix;
  if(stab.stabType==Stab::ARTDIFF){
    if(stab.settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
      kernelName = "cnsGradSurfaceConservative" + suffix; // gradient of all conservative fields
    }else{
      kernelName = "cnsGradSurface" + suffix; // Compute gradient of velocity field
    }
  }else{
    kernelName = "cnsGradSurface" + suffix; 
  }

  gradSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    if (isothermal)
      kernelName = "cnsIsothermalInitialCondition2D";
    else
      kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    if (isothermal)
      kernelName = "cnsIsothermalInitialCondition3D";
    else
      kernelName = "cnsInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  if (isothermal) {
    kernelName = "cnsIsothermalMaxWaveSpeed" + suffix;
  } else {
    kernelName = "cnsMaxWaveSpeed" + suffix;
  }

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);
}
