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

#ifndef STAB_HPP
#define STAB_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "settings.hpp"
#include "mesh.hpp"
#include "ogs.hpp"

namespace libp {

class stabSettings_t: public settings_t {
public:
  stabSettings_t() = default;
  stabSettings_t(comm_t _comm);
  void report();
};

namespace Stab {
  /*Element types*/
  enum SolverType {
    HJS    =1,
    CNS    =2,
    INS    =3,
  };

  enum DetectorType {
    KLOCKNER =1,
    PERSSON  =2,
  };

  enum StabType {
    FILTER  =1,
    LIMITER =2,
    ARTDIFF =3,
    SUBCELL =4,
  };

} //namespace Stab

class stab_t {
 public:
  platform_t platform;
  properties_t props;
  stabSettings_t settings;
  mesh_t mesh; 
  comm_t comm;

  ogs::halo_t traceHalo;

  int dNfields, sNfields; // # of detector and stabiization fields; 

  /*************************/
  /* Solver Data           */
  /*************************/
  Stab::SolverType solverType;
  Stab::DetectorType detectorType;
  Stab::StabType stabType;

  /*************************/
  /* Detector Data        */
  /*************************/
  memory<dlong> eList;
  deviceMemory<dlong> o_eList; 

  memory<dfloat> qd;
  deviceMemory<dfloat> o_qd; 

  // Memory for artificial Viscosity Activation Function 
  // just in case it is needed!!!!
  memory<dfloat> viscRamp; 
  deviceMemory<dfloat> o_viscRamp; 

  // Klockner Detector
  memory<int> modeMap; 
  deviceMemory<int> o_modeMap;
  memory<dfloat> LSF, BLD; 

  deviceMemory<dfloat> o_LSF, o_BLD; 
  deviceMemory<dfloat> o_invV;

  // Persson Detector
  memory<dfloat> projectNm1; 
  deviceMemory<dfloat> o_projectNm1; 

  kernel_t detectKernel; 
  kernel_t copyFloatKernel; 

  /*****************************/
  /*   FILTER STABILIZATION    */
  /*****************************/
  memory<dfloat> filterM; 
  deviceMemory<dfloat> o_filterM; 

  kernel_t filterKernel; 
 
  /*****************************/
  /*   LIMITER STABILIZATION    */
  /*****************************/






  /*******************************************/
  /*   ARTIFICIAL DIFFUSION STABILIZATION    */
  /*******************************************/
  // Memory for artificial Viscosity Activation Function 
  memory<dfloat> viscScale; 
  deviceMemory<dfloat> o_viscScale; 

  memory<dfloat> visc; 
  deviceMemory<dfloat> o_visc; 

  kernel_t computeViscosityKernel; 








  // kernel_t primitiveToConservativeKernel;
  // kernel_t conservativeToPrimitiveKernel;

  // kernel_t MassMatrixKernel;

  stab_t() = default;
  stab_t(platform_t& _platform, mesh_t &_mesh, stabSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  // mesh setup
  void Setup(platform_t& _platform, mesh_t &_mesh, stabSettings_t& _settings);

  // It is safer o split solver implementations I guess
  void DetectorSetup(){
    switch (solverType) {
      case Stab::HJS:
        DetectorSetupHJS();
        break;
      case Stab::CNS:
        // DetectorSetupCNS();
         LIBP_FORCE_ABORT("STAB is not implemented for CNS yet");
        break;
      case Stab::INS:
        LIBP_FORCE_ABORT("STAB is not implemented for INS yet");
        // DetectorSetupINS();
        break;
    } 
  }

  void DetectorSetupHJS(){
    switch (detectorType) {
      case Stab::KLOCKNER:
         DetectorSetupHJSKlockner(); 
        break;
      case Stab::PERSSON:
         DetectorSetupHJSPersson(); 
        break;
    } 
  }


  void DetectorApply(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (solverType) {
      case Stab::HJS:
        DetectorApplyHJS(o_Q, o_RHS, T);
        break;
      case Stab::CNS:
        // DetectorSetupCNS();
         LIBP_FORCE_ABORT("STAB is not implemented for CNS yet");
        break;
      case Stab::INS:
        LIBP_FORCE_ABORT("STAB is not implemented for INS yet");
        // DetectorSetupINS();
        break;
    } 
  }
  

  void DetectorApplyHJS(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (detectorType) {
      case Stab::KLOCKNER:
         DetectorApplyHJSKlockner(o_Q, o_RHS, T); 
        break;
      case Stab::PERSSON:
         DetectorApplyHJSPersson(o_Q, o_RHS, T); 
        break;
    } 
  }


  // It is safer o split solver implementations I guess
  void StabilizerSetup(){
    switch (solverType) {
      case Stab::HJS:
        StabilizerSetupHJS();
        break;
      case Stab::CNS:
        // StabilizerSetupCNS();
         LIBP_FORCE_ABORT("STABILIZER is not implemented for CNS yet");
        break;
      case Stab::INS:
        // StabilizerSetupINS();
        LIBP_FORCE_ABORT("STABILIZER is not implemented for INS yet");
        // DetectorSetupINS();
        break;
    } 
  }


  // It is safer o split solver implementations I guess
  void StabilizerSetupHJS(){
    switch (stabType) {
      case Stab::FILTER:
        StabilizerSetupHJSFilter();
        break;
      case Stab::LIMITER:
         // StabilizerSetupHJSLimiter();
         LIBP_FORCE_ABORT("Limiter is not implemented yet");
        break;
      case Stab::ARTDIFF:
        StabilizerSetupHJSArtdiff();
        break;
      case Stab::SUBCELL:
        // StabilizerSetupHJSSubcell();
        LIBP_FORCE_ABORT("FV-Subcell is not implemented yet");
        break;
    } 
  }




  // It is safer o split solver implementations I guess
  void StabilizerApply(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (solverType) {
      case Stab::HJS:
        StabilizerApplyHJS(o_Q, o_RHS, T);
        break;
      case Stab::CNS:
         LIBP_FORCE_ABORT("STABILIZER is not implemented for CNS yet");
        break;
      case Stab::INS:
        LIBP_FORCE_ABORT("STABILIZER is not implemented for INS yet");
        break;
    } 
  }


// It is safer o split solver implementations I guess
  void StabilizerApplyHJS(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (stabType) {
      case Stab::FILTER:
        StabilizerApplyHJSFilter(o_Q, o_RHS, T);
        break;
      case Stab::LIMITER:
         LIBP_FORCE_ABORT("Limiter is not implemented yet");
        break;
      case Stab::ARTDIFF:
        StabilizerApplyHJSArtdiff(o_Q, o_RHS, T);         
        break;
      case Stab::SUBCELL:
        LIBP_FORCE_ABORT("FV-Subcell is not implemented yet");
        break;
    } 
  }




   void Report(dfloat time, int tstep); 
   dlong GetElementNumber(deviceMemory<dlong>& eList); 
   void PlotElements(memory<dlong> eList, const std::string fileName);
   void PlotFields(memory<dfloat> Q, const std::string fileName);


   void DetectorSetupHJSKlockner(); 
   void DetectorApplyHJSKlockner(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void DetectorSetupHJSPersson(); 
   void DetectorApplyHJSPersson(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 



   void ModeInfoKlocknerTri2D(int _N, memory<int>& _modeMap);
   void ModeInfoKlocknerQuad2D(int _N, memory<int>& _modeMap);
   void ModeInfoKlocknerTet3D(int _N, memory<int>& _modeMap);
   void ModeInfoKlocknerHex3D(int _N, memory<int>& _modeMap);

   void ModeInfoPerssonTri2D(int _N, memory<dfloat>&  _truncModes);
   void ModeInfoPerssonQuad2D(int _N, memory<dfloat>& _truncModes);
   void ModeInfoPerssonTet3D(int _N, memory<dfloat>&  _truncModes);
   void ModeInfoPerssonHex3D(int _N, memory<dfloat>&  _truncModes);

   void LeastSquaresFitKlockner(int _N, memory<dfloat>& _LSF);
   void BaseLineDecayKlockner(int _N, memory<dfloat>& _BLD);



   void FilterMatrixTri2D (int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);
   void FilterMatrix1D(int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);
   void FilterMatrixTet3D (int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);



   void StabilizerSetupHJSFilter(); 
   void StabilizerApplyHJSFilter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void StabilizerSetupHJSArtdiff(); 
   void StabilizerApplyHJSArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
    
   dfloat ElementViscosityScaleTri2D(dlong e);
   dfloat ElementViscosityScaleQuad2D(dlong e);
   dfloat ElementViscosityScaleTet3D(dlong e);
   dfloat ElementViscosityScaleHex3D(dlong e);


 
 private:
  /*Set the type of mesh*/
  void SetTypes(const Stab::SolverType sType, 
                const Stab::DetectorType dType, 
                const Stab::StabType stType);

 
};

} //namespace libp

#endif

