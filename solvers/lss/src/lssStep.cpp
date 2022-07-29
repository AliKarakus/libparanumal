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
    Redistance(o_Q, o_RHS, T);
  }

}

//evaluate ODE rhs = f(q,t)

void lss_t::postStep(deviceMemory<dfloat>& o_Q, const dfloat time, const dfloat dt){


if(redistance){
// if(timeStepper.outputStep==0){
   historyIndex +=1; 

  // const dfloat dt = timeStepper->GetTimeStep(); 
  reconstructTime[shiftIndex] = time + dt; 
  o_reconstructTime.copyFrom(reconstructTime); 


  const dlong N = mesh.Nelements*mesh.Np*Nfields; 
  o_phiH.copyFrom(o_Q, N, shiftIndex*N,0); 

      // std::cout<<"history index: "<<historyIndex<<" "<<time<<" "<<dt<<" "<< Nrecon<<std::endl;


   if(historyIndex==(Nrecon/2)){
    std::cout<<"history index: "<<historyIndex<<" "<<time<<" "<<dt<<" "<< Nrecon<<std::endl;
    // Create Initial History
    timeInitialHistoryKernel(mesh.Nelements, 
                             shiftIndex, 
                             N, 
                             o_phiH);
   }



  if(historyIndex>=(Nrecon/2))
  timeReconstructKernel(mesh.Nelements, 
                        shiftIndex, 
                        N, 
                        o_reconstructTime,
                        o_phiH,
                        o_phi);

  shiftIndex = (shiftIndex+Nrecon-1)%Nrecon;
  // std::cout<<shiftIndex<<std::endl;

 // DetectTroubledCells(o_Q, subcell->o_ElementList); 

// }


}


}

void lss_t::postStage(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat> & o_RHS, const dfloat T){

// std::cout<<"lssPostStage is called"<<std::endl; 

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

