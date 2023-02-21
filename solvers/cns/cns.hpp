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

#ifndef CNS_HPP
#define CNS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "stab.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DCNS LIBP_DIR"/solvers/cns/"

using namespace libp;

class cnsSettings_t: public settings_t {
public:
  cnsSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     stabSettings_t& stabSettings,
                     const std::string filename);
};


// namespace  {
//   /*Viscosity Model*/
//   enum ViscosityType {
//     CONSTANT      =1,
//     SUTHERLAND    =2,
//     POWER         =3,
//   };
// } //namespace Stab

class cns_t: public solver_t {
public:
  mesh_t mesh;
  stab_t stab;

  int Nfields;
  int Ngrads;

  dfloat mu;
  dfloat gamma;

  dfloat Ma;  // Mach Number
  dfloat Re;  // Reynolds Number
  dfloat R;   // Gas Constant
  dfloat cp;  // Heat capacity for constant pressure
  dfloat cv;  // Heat capacity for constant volume
  dfloat Pr;  // Prandtl Number
  dfloat Tref;  // Reference Temperature

  int Nph; // number of physical parameters 
  int MUID, GMID, RRID, PRID, CPID, CVID; 
  int EXID, TRID, TSID, CSID; 


  int cubature;
  int isothermal;
  int inviscid; 

  int useNonDimensionalEqn; 
  int viscType; 

  memory<dfloat> pCoeff; //  Physical coefficients 
  deviceMemory<dfloat> o_pCoeff; 

  timeStepper_t timeStepper;

  ogs::halo_t fieldTraceHalo;
  ogs::halo_t gradTraceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  memory<dfloat> gradq;
  deviceMemory<dfloat> o_gradq;

  memory<dfloat> Vort;
  deviceMemory<dfloat> o_Vort;

  deviceMemory<dfloat> o_Mq;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t cubatureVolumeKernel;
  kernel_t cubatureSurfaceKernel;

  kernel_t gradVolumeKernel;
  kernel_t gradSurfaceKernel;

  kernel_t vorticityKernel;

  kernel_t constrainKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  cns_t() = default;
  cns_t(platform_t &_platform, mesh_t &_mesh, stab_t &_stab, 
              cnsSettings_t& _settings) {
    Setup(_platform, _mesh, _stab, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh, stab_t &_stab,
             cnsSettings_t& _settings);

  void Run();

  void Report(dfloat time, int tstep) override;

  void PlotFields(memory<dfloat> Q, memory<dfloat> V, std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  void rhsNoStab(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void rhsArtDiff(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);

  void postStep(deviceMemory<dfloat>& o_q,  const dfloat time, const dfloat dt);
  void postStage(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time, const dfloat dt);
};

#endif
