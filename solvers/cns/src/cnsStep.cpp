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


//evaluate ODE rhs = f(q,t)
void cns_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
   if(stab.stabType==Stab::ARTDIFF){
      rhsArtDiff(o_Q, o_RHS, T); 
      // rhsNoStab(o_Q, o_RHS, T); 
    }else{
      rhsNoStab(o_Q, o_RHS, T); 
   }
}

dfloat cns_t::MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T){

  //Note: if this is on the critical path in the future, we should pre-allocate this
  deviceMemory<dfloat> o_maxSpeed = platform.malloc<dfloat>(mesh.Nelements);

  maxWaveSpeedKernel(mesh.Nelements,
                     mesh.o_vgeo,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_EToB,
                     o_pCoeff,
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_Q,
                     o_maxSpeed);

  const dfloat vmax = platform.linAlg().max(mesh.Nelements, o_maxSpeed, mesh.comm);

  return vmax;
}


void cns_t::postStage(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat> & o_sQ, const dfloat T, const dfloat DT){


}



void cns_t::postStep(deviceMemory<dfloat>& o_Q, const dfloat T, const dfloat dt){
  

}


//evaluate ODE rhs = f(q,t)
void cns_t::rhsArtDiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

{

  if(T>0.1){

  stab.detectApply(o_Q, o_RHS, T); 

  // const dfloat alpha = 0.025; 
  const dfloat alpha = 1.0; 
  stab.computeViscosityKernel(mesh.Nelements*mesh.Nverts, 
                             alpha, 
                             stab.o_viscActivation,
                             stab.o_viscScale,
                             stab.o_vertexVisc); 
  
  // stab.Report(time,tstep);
  stab.ogs.GatherScatter(stab.o_vertexVisc, 1, ogs::Add, ogs::Sym); 
  platform.linAlg().amx(mesh.Nelements*mesh.Nverts, 1, stab.o_weight, stab.o_vertexVisc); 

  stab.projectViscosityKernel(mesh.Nelements, 
                             stab.o_projectVisc,
                             stab.o_vertexVisc,
                             stab.o_visc); 

  // const dfloat zero = 0.0000; 
  // platform.linAlg().set(mesh.Nelements*mesh.Np, zero, stab.o_visc);

  // //Add physical viscosity to artificial one
  // platform.linAlg().add(mesh.Nelements*mesh.Np, mu, stab.o_visc); 

  // Change it later AK
   // mesh.halo.Exchange(stab.o_visc, mesh.Np);

   // stab.Report(T,T);
// 
   // Report(T, T);
   }else{
        // set time step
      dfloat hmin = mesh.MinCharacteristicLength();
      const dfloat visc = 1.0* hmin/mesh.N; 
      platform.linAlg().set(mesh.Nelements*mesh.Np, visc, stab.o_visc);


   }

  }

  // Now compute the viscous flux : Fv = mu* grad{q} 
  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

  // compute volume contributions to gradients
  gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_Q,
                   o_gradq);

  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    o_pCoeff, 
                    o_Q,
                    o_gradq);

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);

  // compute volume contribution to cns RHS
  if (cubature) {
    cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         T,
                         o_pCoeff, 
                         stab.o_visc,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 T,
                 mu,
                 stab.o_visc,
                 gamma,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
    // dfloat hmin = mesh.MinCharacteristicLength();
    const dfloat tau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      cubatureSurfaceKernel(mesh.Nelements,
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
                            T,
                            o_pCoeff, 
                            tau,
                            stab.o_visc,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      // dfloat hmin = mesh.MinCharacteristicLength();
    const dfloat tau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      surfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    tau,
                    stab.o_visc,
                    gamma,
                    o_Q,
                    o_gradq,
                    o_RHS);
    }


  // Report(T, T);
}


//evaluate ODE rhs = f(q,t)
void cns_t::rhsNoStab(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

  // compute volume contributions to gradients
  gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_Q,
                   o_gradq);

  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    gamma,
                    o_Q,
                    o_gradq);

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);

  // compute volume contribution to cns RHS
  if (cubature) {
    cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         T,
                         mu,
                         gamma,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 T,
                 mu,
                 gamma,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
      cubatureSurfaceKernel(mesh.Nelements,
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
                            T,
                            mu,
                            gamma,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      surfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    gamma,
                    o_Q,
                    o_gradq,
                    o_RHS);
    }
}

