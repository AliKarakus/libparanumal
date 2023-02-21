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

void cns_t::Report(dfloat time, int tstep){


  {

  stab.detectApply(o_q, o_q, time); 

  // const dfloat alpha = 0.0005; 
  const dfloat alpha = 1.0; 
  stab.computeViscosityKernel(mesh.Nelements*mesh.Nverts, 
                             alpha, 
                             stab.o_viscActivation,
                             stab.o_viscScale,
                             stab.o_vertexVisc); 
  
  // // stab.Report(time,tstep);
  stab.ogs.GatherScatter(stab.o_vertexVisc, 1, ogs::Add, ogs::Sym); 
  platform.linAlg().amx(mesh.Nelements*mesh.Nverts, 1, stab.o_weight, stab.o_vertexVisc); 

  stab.projectViscosityKernel(mesh.Nelements, 
                             stab.o_projectVisc,
                             stab.o_vertexVisc,
                             stab.o_visc); 
  stab.Report(time,tstep);
  }

  // stab.Report(time,time);

  

  static int frame=0;

  //compute vorticity
  vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, o_Vort);

  //compute q.M*q
  mesh.MassMatrixApply(o_q, o_Mq);

  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // copy data back to host
    o_q.copyTo(q);

    
   //  {
   //  char fname[BUFSIZ];
   //  sprintf(fname, "pressure_%d.dat", frame);
   //  std::string filename = std::string(fname); 
   //  FILE *fp = fopen(filename.c_str(), "w");

   //   for(int e=0; e< mesh.Nelements; e++){
   //    for(int f = 0; f< mesh.Nfaces; f++){

   //      int bc = mesh.EToB[e*mesh.Nfaces + f]; 

   //      if(bc == 1){
   //      for(int n=0; n<mesh.Nfp; n++){

   //        const dlong id = e*mesh.Nfaces*mesh.Nfp + f*mesh.Nfp + n; 
   //        const dlong vid = mesh.vmapM[id]; 

   //        // load traces
   //        const dlong eM = e;
   //        const int vidM = vid%mesh.Np;
   //        const dlong qbase = eM*mesh.Np*Nfields + vidM;

   //        const dfloat r  = q[qbase+0*mesh.Np];
   //        const dfloat ru = q[qbase+1*mesh.Np];
   //        const dfloat rv = q[qbase+2*mesh.Np];
   //        const dfloat E  = q[qbase+3*mesh.Np];

   //        // primitive variables
   //        const dfloat u = ru/r, v = rv/r;
   //        const dfloat pn = (gamma-1)*(E-0.5*r*(u*u+v*v));
   //        const dfloat xn = mesh.x[eM*mesh.Np + vidM]; 
   //        const dfloat yn = mesh.y[eM*mesh.Np + vidM]; 
   //        fprintf(fp, "%.8e %.8e %.8e\n", xn, yn, pn);

   //      }
   //      }
   //    }
   //   }
   //     fclose(fp);
   // }





    o_Vort.copyTo(Vort);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    PlotFields(q, Vort, std::string(fname));
  }
}
