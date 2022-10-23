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

void lss_t::Report(dfloat time, int tstep){

  static int frame=0;
  dfloat norm2 = 0; 
  if(advection){
    //compute q.M*q
    mesh.MassMatrixApply(o_q, o_Mq);\

    dlong Nentries = mesh.Nelements*mesh.Np;
    norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, comm));
  }else{
    //compute q.M*q
    mesh.MassMatrixApply(o_phi, o_Mq);

    dlong Nentries = mesh.Nelements*mesh.Np;
    norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_phi, o_Mq, comm));
  }

  if(mesh.rank==0)
    printf("%5.2f (%d), %.8e (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    // copy data back to host
    if(advection){
      o_q.copyTo(q);
      o_U.copyTo(U);
      PlotFields(q, std::string(fname));
    }else{
      // o_phi.copyTo(phi);
      dlong Nentries = mesh.Nelements*mesh.Np;
      // deviceMemory<dfloat> o_test; o_test = platform.malloc<dfloat>(); 
      o_gsphi.copyFrom(o_phi, Nentries);
      lssogs.GatherScatter(o_gsphi, 1, ogs::Add, ogs::Sym); 

      platform.linAlg().amx(Nentries, 1, o_weight, o_gsphi); 
      o_gsphi.copyTo(phi); 
      
      o_q.copyTo(q);
      PlotFields(phi, std::string(fname));
      stab.Report(time,tstep); 
    }
    
    

  }
}
