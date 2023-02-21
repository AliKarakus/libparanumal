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

void cns_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         mu,
                         gamma,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed(o_q, startTime);

  dfloat dtAdv  = cfl/(vmax*(mesh.N+1.)*(mesh.N+1.));
  
  // For Euler tests optout viscous time resitriction 
  dfloat dtVisc = mu > 1E-12 ? cfl*pow(hmin, 2)/(pow(mesh.N+1,4)*mu): 1E12;
  // dfloat dtVisc = mu > 1E-12 ? cfl*pow(hmin, 2)/(pow(mesh.N+1,4)*hmin/mesh.N): 1E12;

  dfloat dt = std::min(dtAdv, dtVisc);

  printf("time step size: %.4e", dt);
  timeStepper.SetTimeStep(dt);


 //  stab.detectApply(o_q, o_q, 0.0); 

 //  const dfloat alpha = 1.0; 
 //  stab.computeViscosityKernel(mesh.Nelements*mesh.Nverts, 
 //                             alpha, 
 //                             stab.o_viscActivation,
 //                             stab.o_viscScale,
 //                             stab.o_vertexVisc); 

 // // stab.Report(0,0);

 //  stab.ogs.GatherScatter(stab.o_vertexVisc, 1, ogs::Add, ogs::Sym); 
 //  platform.linAlg().amx(mesh.Nelements*mesh.Nverts, 1, stab.o_weight, stab.o_vertexVisc); 

 //  stab.projectViscosityKernel(mesh.Nelements, 
 //                             stab.o_projectVisc,
 //                             stab.o_vertexVisc,
 //                             stab.o_visc); 
    


 //  stab.Report(0,0);

 //  Report(0,0);

  timeStepper.Run(*this, o_q, startTime, finalTime);

  // output norm of final solution
  {
    //compute q.M*q
    mesh.MassMatrixApply(o_q, o_Mq);

    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }

}
