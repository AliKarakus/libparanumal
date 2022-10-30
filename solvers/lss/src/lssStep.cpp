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

//evaluate ODE rhs = f(q,t)
void lss_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  // Solve Level-Set Advection  
  if(advection){
    Advection(o_Q, o_RHS, T);
  }else if(redistance){
    if(stab.stabType==Stab::FILTER){
      RedistanceFilter(o_Q, o_RHS, T);
    }else if(stab.stabType==Stab::ARTDIFF){
      RedistanceArtdiff(o_Q, o_RHS, T);
    }else if(stab.stabType==Stab::SUBCELL){
      RedistanceSubcell(o_Q, o_RHS, T);      
    }else{ // Default to no-stabilization
      Redistance(o_Q, o_RHS, T);
    }
  }

}

// //evaluate ODE rhs = f(q,t)

// void lss_t::reconstruct(deviceMemory<dfloat>& o_Q, const dfloat T){
  
//   bool outstep = timeStepper.isOutputStep(); 

//     if(outstep==false){
//       // Hold history of recontruction times
//       static int shiftIndex   = 0; 
//       static int historyIndex = 0; 
//       const dlong N = mesh.Nelements*mesh.Np*Nfields; 

//       reconstructTime[shiftIndex] = T; 
//       o_reconstructTime.copyFrom(reconstructTime); 

//       o_phiH.copyFrom(o_Q, N, shiftIndex*N, 0); 

//       if(historyIndex < (Nrecon/2) ){
//         std::cout<<"creating initial history: index: "<<historyIndex<<" at time "<<T<<std::endl;
//         // Create Initial History
//         timeInitialHistoryKernel(mesh.Nelements, 
//                                  shiftIndex, 
//                                  N, 
//                                  o_phiH);
//       }

//       if(historyIndex>=(Nrecon/2)){
//         timeReconstructKernel(mesh.Nelements, 
//                             shiftIndex, 
//                             N, 
//                             o_reconstructTime,
//                             o_phiH,
//                             o_phi);
//       }

//       // update indices
//       historyIndex +=1; 
//       shiftIndex= (shiftIndex+Nrecon-1)%Nrecon;
//   }

// }


//evaluate ODE rhs = f(q,t)

void lss_t::postStep(deviceMemory<dfloat>& o_Q, const dfloat T, const dfloat dt){
  
  bool outstep = timeStepper.isOutputStep(); 

    if(outstep==false){
      // Hold history of recontruction times
      static int shiftIndex   = 0; 
      static int historyIndex = 0; 
      const dlong N = mesh.Nelements*mesh.Np*Nfields; 

      // const dfloat dt = timeStepper.GetTimeStep();

      reconstructTime[shiftIndex] = T + dt; 
      o_reconstructTime.copyFrom(reconstructTime); 

      o_phiH.copyFrom(o_Q, N, shiftIndex*N, 0); 

      if(historyIndex < (Nrecon/2) ){
        std::cout<<"creating initial history: index: "<<historyIndex<<" at time "<<T<<std::endl;
        // Create Initial History
        timeInitialHistoryKernel(mesh.Nelements, 
                                 shiftIndex, 
                                 N, 
                                 o_phiH);
      }

      if(historyIndex>=(Nrecon/2)){
        timeReconstructKernel(mesh.Nelements, 
                            shiftIndex, 
                            N, 
                            o_reconstructTime,
                            o_phiH,
                            o_phi);
      }

      // update indices
      historyIndex +=1; 
      shiftIndex= (shiftIndex+Nrecon-1)%Nrecon;

      // stab.detectApply(o_Q, o_Q, T); 
  }

}


void lss_t::postStage(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat> & o_sQ, const dfloat T, const dfloat DT){


// if(stab.stabType==Stab::SUBCELL){
//  reconstructDGKernel(mesh.Nelements, 
//                      stab.o_eList,
//                      stab.o_RM,
//                      o_sQ, 
//                      o_Q);
// }

}



void lss_t::Advection(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
 
// extract q halo on DEVICE
  qTraceHalo.ExchangeStart(o_Q, 1);

  // This replaces flow solver
  setFlowFieldKernel(mesh.Nelements,  
                    T, 
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_U);

  // just in case velocity is not communicated
  vTraceHalo.ExchangeStart(o_U, 1);


  if(cubature){
     advectionVolumeKernel(mesh.Nelements,
                         T,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_cubx,
                         mesh.o_cuby,
                         mesh.o_cubz,
                         o_U,
                         o_Q,
                         o_RHS);
  }else{
  advectionVolumeKernel(mesh.Nelements,
                       T,
                       mesh.o_vgeo,
                       mesh.o_D,
                       mesh.o_x,
                       mesh.o_y,
                       mesh.o_z,
                       o_U,
                       o_Q,
                       o_RHS);
 }

 if(compressible){
// vTraceHalo.ExchangeStart(o_U, 1);
// Assumes velocity field is already communicated
compressibleVolumeKernel(mesh.Nelements,
                       T,
                       mesh.o_vgeo,
                       mesh.o_D,
                       o_U,
                       o_Q,
                       o_RHS);
}




  qTraceHalo.ExchangeFinish(o_Q, 1);
  vTraceHalo.ExchangeFinish(o_U, 1);
  
  if(cubature){
    advectionSurfaceKernel(mesh.Nelements,
                            T,
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            o_U,
                            o_Q,
                            o_RHS);

  }else{
  advectionSurfaceKernel(mesh.Nelements,
                mesh.o_sgeo,
                mesh.o_LIFT,
                mesh.o_vmapM,
                mesh.o_vmapP,
                mesh.o_EToB,
                T,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_U,
                o_Q,
                o_RHS);
 }

if(compressible){
compressibleSurfaceKernel(mesh.Nelements,
                mesh.o_sgeo,
                mesh.o_LIFT,
                mesh.o_vmapM,
                mesh.o_vmapP,
                mesh.o_EToB,
                T,
                mesh.o_x,
                mesh.o_y,
                mesh.o_z,
                o_U,
                o_Q,
                o_RHS);
}
}


// This uses strong form implementation
void lss_t::Redistance(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
 
  // extract q halo on DEVICE
  qTraceHalo.ExchangeStart(o_Q, 1);

 redistanceVolumeKernel(mesh.Nelements,
                       T,
                       mesh.o_vgeo,
                       mesh.o_DW,
                       o_Q,
                       o_gradq);

  qTraceHalo.ExchangeFinish(o_Q, 1);

  redistanceSurfaceKernel(mesh.Nelements,
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_Q,
                          o_gradq,
                          o_RHS);
}



// This form uses strong forms
void lss_t::RedistanceFilter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
 
// // Recontruct phi from the history
// reconstruct(o_Q, T); 

// Detect troubled elements
stab.detectApply(o_Q, o_RHS, T); 

  // Filter
  filterKernel(mesh.Nelements, 
               stab.o_eList, 
               stab.o_filterM, 
               o_Q); 
 // }

  // extract q halo on DEVICE
  qTraceHalo.ExchangeStart(o_Q, 1);

 redistanceVolumeKernel(mesh.Nelements,
                       T,
                       mesh.o_vgeo,
                       mesh.o_D,
                       o_Q,
                       o_gradq);

  qTraceHalo.ExchangeFinish(o_Q, 1);

  redistanceSurfaceKernel(mesh.Nelements,
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_Q,
                          o_gradq,
                          o_RHS);
}



// This solver uses strong from
void lss_t::RedistanceArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  

  // // Recontruct phi from the history
  // reconstruct(o_Q, T); 

  // Detect troubled elements
  stab.detectApply(o_Q, o_RHS, T); 


  const dfloat alpha = 3.0; 
  stab.computeViscosityKernel(mesh.Nelements, 
                              alpha, 
                              stab.o_viscRamp,
                              stab.o_viscScale,
                              stab.o_visc); 

  // extract q halo on DEVICE
  qTraceHalo.ExchangeStart(o_Q, 1);

  // compute gradient volume: should be in strong form
 redistanceVolumeKernel(mesh.Nelements,
                       T,
                       mesh.o_vgeo,
                       mesh.o_D,
                       o_Q,
                       o_gradq);

  // We can combine q and gradq : AK
  qTraceHalo.ExchangeFinish(o_Q, 1);

  // Surface Kernel
  redistanceSurfaceKernel(mesh.Nelements,
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_Q,
                          o_gradq,
                          o_RHS);

// // // *********************************************************/
//   // // qTraceHalo.ExchangeStart(o_Q, 1);

//   // // Compute Volume Contribution
//   // gradientVolumeKernel(mesh.Nelements,
//   //                     mesh.o_vgeo,
//   //                     mesh.o_D,
//   //                     o_Q,
//   //                     o_dGq);

  // qTraceHalo.ExchangeFinish(o_Q, 1);

  // Add Surface Compute Surface Conribution
  gradientSurfaceKernel(mesh.Nelements,
                       mesh.o_sgeo,
                       mesh.o_LIFT,
                       mesh.o_vmapM,
                       mesh.o_vmapP,
                       mesh.o_EToB,
                       T,
                       mesh.o_x,
                       mesh.o_y,
                       mesh.o_z,
                       o_Q,
                       o_gradq);

  gTraceHalo.ExchangeStart(o_gradq, 1); 

  divergenceVolumeKernel(mesh.Nelements,
                        mesh.o_vgeo,
                        mesh.o_D,
                        stab.o_visc,
                        o_gradq,
                        o_RHS); // ?

  gTraceHalo.ExchangeFinish(o_gradq, 1); 


  divergenceSurfaceKernel(mesh.Nelements,
                           mesh.o_sgeo,
                           mesh.o_LIFT,
                           mesh.o_vmapM,
                           mesh.o_vmapP,
                           mesh.o_EToB,
                           qTau,
                           T,
                           mesh.o_x,
                           mesh.o_y,
                           mesh.o_z,
                           stab.o_visc,
                           o_Q,
                           o_gradq,
                           o_RHS);
}





void lss_t::RedistanceSubcell(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

// Detect troubled elements
stab.detectApply(o_Q, o_RHS, T); 

// qTraceHalo.ExchangeStart(o_Q, 1);
mesh.halo.ExchangeStart(o_Q, Nfields*mesh.Np);
// mesh.halo.Exchange(o_Q, Nfields*mesh.Np);


 redistanceVolumeKernel(mesh.Nelements,
                       T,
                       stab.o_eList, 
                       mesh.o_vgeo,
                       mesh.o_DW,
                       o_Q,
                       o_gradq);


 mesh.halo.ExchangeFinish(o_Q, Nfields*mesh.Np);

  redistanceSurfaceKernel(mesh.Nelements,
                          stab.o_eList, 
                          mesh.o_sgeo,
                          mesh.o_LIFT,   
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_Q,
                          o_gradq,
                          o_RHS);

 // Get subcell face values 
  projectDGKernel(mesh.Nelements,
                   stab.o_eList, 
                   stab.o_mFaceNodes, // mesh.o_vmapM,    
                   stab.o_mFToE,
                   stab.o_mFToF,
                   stab.o_PFM,
                   o_Q,
                   o_sface);

  mesh.halo.Exchange(o_sface, stab.sNfields*stab.Nsubcells*stab.Nfaces);


  // Get cell averages from nodal solution
  projectFVKernel(mesh.Nelements + mesh.totalHaloPairs,
                          stab.o_eList, 
                          stab.o_PM,
                          o_Q,
                          o_sq);


// mesh.halo.Exchange(o_sq, stab.sNfields*stab.Nsubcells);

// Reconstruct face values for all subcells 
     reconstructFaceKernel(mesh.Nelements, 
                          stab.o_eList,
                          stab.o_vgeo,
                          stab.o_sgeo, 
                          stab.o_emapP,
                          stab.o_fmapP, 
                          o_Q, 
                          o_sq,
                          o_sface);  

  // mesh.halo.Exchange(o_sq, stab.sNfields*stab.Nsubcells);
  // mesh.halo.Exchange(o_sface, stab.sNfields*stab.Nsubcells*stab.Nfaces);

   // FV compute 
 subcellComputeKernel(mesh.Nelements, 
                          stab.o_eList,
                          stab.o_emapP,
                          stab.o_fmapP, 
                          stab.o_RM,
                          stab.o_vgeo,
                          stab.o_sgeo, 
                          o_Q,
                          o_sface, 
                          o_srhs);  


 reconstructDGKernel(mesh.Nelements, 
                     stab.o_eList,
                     stab.o_RM,
                     o_srhs, 
                     o_RHS);







// // Detect troubled elements
// stab.detectApply(o_Q, o_RHS, T); 

// // extract q halo on DEVICE
//   qTraceHalo.ExchangeStart(o_Q, 1);

//  redistanceVolumeKernel(mesh.Nelements,
//                        T,
//                        stab.o_eList, 
//                        mesh.o_vgeo,
//                        mesh.o_DW,
//                        o_Q,
//                        o_gradq);

//   qTraceHalo.ExchangeFinish(o_Q, 1);

//   // Get cell averages from nodal solution
//   projectDGKernel((mesh.Nelements + mesh.totalHaloPairs),
//                    stab.o_eList, 
//                    mesh.o_vmapM,
//                    stab.o_mFToE,
//                    stab.o_mFToF,
//                    stab.o_PFM,
//                    o_Q,
//                    o_sface);

//   redistanceSurfaceKernel(mesh.Nelements,
//                           stab.o_eList, 
//                           mesh.o_sgeo,
//                           mesh.o_LIFT,
//                           mesh.o_vmapM,
//                           mesh.o_vmapP,
//                           mesh.o_EToB,
//                           T,
//                           mesh.o_x,
//                           mesh.o_y,
//                           mesh.o_z,
//                           o_Q,
//                           o_gradq,
//                           o_RHS);

//     // Get cell averages from nodal solution
//     projectFVKernel((mesh.Nelements + mesh.totalHaloPairs),
//                   stab.o_eList, 
//                   stab.o_PM,
//                   o_Q,
//                   o_sq);


//      // Reconstruct face values for all subcells 
//      reconstructFaceKernel(mesh.Nelements, 
//                           stab.o_eList,
//                           stab.o_vgeo,
//                           stab.o_sgeo, 
//                           stab.o_emapP,
//                           stab.o_fmapP, 
//                           o_Q, 
//                           o_sq,
//                           o_sface);  


//        // FV compute 
//      subcellComputeKernel(mesh.Nelements, 
//                               stab.o_eList,
//                               stab.o_emapP,
//                               stab.o_fmapP, 
//                               stab.o_RM,
//                               stab.o_vgeo,
//                               stab.o_sgeo, 
//                               o_Q,
//                               o_sface, 
//                               o_srhs);  


//  reconstructDGKernel(mesh.Nelements, 
//                      stab.o_eList,
//                      stab.o_RM,
//                      o_srhs, 
//                      o_RHS);
}



void lss_t::rhsf_subcell(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_sQ,
                             deviceMemory<dfloat>& o_RHS, deviceMemory<dfloat>& o_sRHS, const dfloat T){

  // Detect troubled elements
  stab.detectApply(o_Q, o_RHS, T); 

  // extract q halo on DEVICE
  qTraceHalo.ExchangeStart(o_Q, 1);

 redistanceVolumeKernel(mesh.Nelements,
                       T,
                       stab.o_eList, 
                       mesh.o_vgeo,
                       mesh.o_DW,
                       o_Q,
                       o_gradq);

  qTraceHalo.ExchangeFinish(o_Q, 1);

  // Get cell averages from nodal solution
  projectDGKernel((mesh.Nelements + mesh.totalHaloPairs),
                   stab.o_eList, 
                   mesh.o_vmapM,
                   stab.o_mFToE,
                   stab.o_mFToF,
                   stab.o_PFM,
                   o_Q,
                   o_sface);

  redistanceSurfaceKernel(mesh.Nelements,
                          stab.o_eList, 
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          o_Q,
                          o_gradq,
                          o_RHS);

    // Get cell averages from nodal solution
    projectFVKernel((mesh.Nelements + mesh.totalHaloPairs),
                  stab.o_eList, 
                  stab.o_PM,
                  o_Q,
                  o_sQ);


     // Reconstruct face values for all subcells 
     reconstructFaceKernel(mesh.Nelements, 
                          stab.o_eList,
                          stab.o_vgeo,
                          stab.o_sgeo, 
                          stab.o_emapP,
                          stab.o_fmapP, 
                          o_Q, 
                          o_sQ,
                          o_sface);  


     // FV compute 
    subcellComputeKernel(mesh.Nelements, 
                            stab.o_eList,
                            stab.o_emapP,
                            stab.o_fmapP, 
                            stab.o_RM,
                            stab.o_vgeo,
                            stab.o_sgeo, 
                            o_Q,
                            o_sface, 
                            o_sRHS);  


     reconstructDGKernel(mesh.Nelements, 
                         stab.o_eList,
                         stab.o_RM,
                         o_sRHS, 
                         o_RHS);



}

