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

#ifndef LSS_HPP
#define LSS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "stab.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DLSS LIBP_DIR"/solvers/lss/"

using namespace libp;

class lssSettings_t: public settings_t {
public:
  lssSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     stabSettings_t& stabSettings,
                     const std::string filename);
};

class lss_t: public solver_t {
public:
  int cubature; 
  int advection;
  int redistance; 
  int compressible; // Type of flow field
  int Nfields; // Number of fields for level set 
  int Nrecon; // Number of history points for reconstruction  
  
  // dlong offset; 
  mesh_t mesh;
  stab_t stab; 
  timeStepper_t timeStepper;

  ogs::halo_t qTraceHalo;
  ogs::halo_t vTraceHalo;
  ogs::halo_t gTraceHalo;

  ogs::halo_t qHalo; 

  memory<dfloat> q;
  memory<dfloat> U; 
  memory<dfloat> gradq;

  deviceMemory<dfloat> o_q;
  deviceMemory<dfloat> o_U;
  deviceMemory<dfloat> o_gradq;

  memory<dfloat> phi; 
  deviceMemory<dfloat> o_phi; 

  memory<dfloat> phiH; 
  deviceMemory<dfloat> o_phiH; 

  memory<dfloat> reconstructTime; 
  deviceMemory<dfloat> o_reconstructTime; 

  deviceMemory<dfloat> o_Mq;

  kernel_t advectionVolumeKernel;
  kernel_t advectionSurfaceKernel;

  kernel_t advectionCubVolumeKernel;
  kernel_t advectionCubSurfaceKernel;

  // add compressible terms 
  kernel_t compressibleVolumeKernel;
  kernel_t compressibleSurfaceKernel;

  // Redistance
  kernel_t redistanceVolumeKernel;
  kernel_t redistanceSurfaceKernel;
  kernel_t timeReconstructKernel;
  kernel_t timeInitialHistoryKernel;
  kernel_t redistanceSetFieldsKernel; 

  kernel_t MassMatrixKernel;

  kernel_t initialConditionKernel;
  kernel_t setFlowFieldKernel;

  // Stabilization Related
  kernel_t filterKernel; 

  // ***************************ARTDIFF*********************//
  dfloat qTau; // Stabilization Parameter
  kernel_t gradientVolumeKernel;
  kernel_t gradientSurfaceKernel;
  kernel_t divergenceVolumeKernel; 
  kernel_t divergenceSurfaceKernel; 
  // ***************************SUBCELL*********************//
  ogs::ogs_t lssogs;
  deviceMemory<dfloat> o_gsphi; 

  memory<dfloat> sface; 
  deviceMemory<dfloat> o_sface;

  memory<dfloat> srhs; 
  deviceMemory<dfloat> o_srhs;

  memory<dfloat> sq; 
  deviceMemory<dfloat> o_sq;


  memory<dfloat> weight; 
  deviceMemory<dfloat> o_weight; 

  kernel_t projectFVKernel; // project DG -> Cell Centers
  kernel_t projectDGKernel; // project DG -> Cell Faces at DG-FV Interface

  kernel_t partialProjectFVKernel; // project DG -> Cell Centers
  kernel_t partialProjectDGKernel; // project DG -> Cell Faces at DG-FV Interface

  kernel_t reconstructFaceKernel; // compute face values
  kernel_t subcellComputeKernel; // FV Update values
  kernel_t reconstructDGKernel; // recontruct DG solution from cell averages values

 

  lss_t() = default;
  lss_t(platform_t &_platform, mesh_t &_mesh, stab_t _stab, 
              lssSettings_t& _settings) {
    Setup(_platform, _mesh, _stab, _settings);
  }

  //setup
  void Setup(platform_t& platform, mesh_t& mesh, stab_t _stab, 
             lssSettings_t& settings);

  void Run();

  void Advection(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
  void Redistance(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
  void RedistanceFilter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
  void RedistanceArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
  void RedistanceSubcell(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

  void rhsf_subcell(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, 
                    deviceMemory<dfloat>& o_sQ, deviceMemory<dfloat>& o_sRHS, const dfloat T); 
 
  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q, const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);
  void reconstruct(deviceMemory<dfloat>& o_q, const dfloat time);


  void postStep(deviceMemory<dfloat>& o_q,  const dfloat time, const dfloat dt);
  void postStage(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time, const dfloat dt);

  
};


#endif

