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

  /***************************/
  int CXID, CYID, CZID, IVID; 
  int FXID, FYID, FZID, SAID; 
  int NXID, NYID, NZID, BCID; 

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

  // Remove Later
  memory<dfloat> efList;
  deviceMemory<dfloat> o_efList; 
  kernel_t copyIntKernel; 

  memory<dfloat> qd;
  deviceMemory<dfloat> o_qd; 

  
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
  kernel_t extractFieldKernel; 

  /*****************************/
  /*   FILTER STABILIZATION    */
  /*****************************/
  memory<dfloat> filterM; 
  deviceMemory<dfloat> o_filterM; 

  // kernel_t filterKernel; 
 
  /*****************************/
  /*   LIMITER STABILIZATION    */
  /*****************************/
  memory<dfloat> qv; 
  deviceMemory<dfloat> o_qv; 

  memory<dfloat> qc; 
  deviceMemory<dfloat> o_qc; 

  memory<dfloat> dq; 
  deviceMemory<dfloat> o_dq; 

  memory<dfloat> dqf; 
  deviceMemory<dfloat> o_dqf; 

  // Project to cell centers
  memory<dfloat> projectC0; 
  deviceMemory<dfloat> o_projectC0; 

  memory<int> vertexNodes;
  deviceMemory<int> o_vertexNodes; 

  // differences with respect to cell-center
  // required for Taylor expansion 
  memory<dfloat> dAVE; 
  
  memory<dfloat> DX; 
  deviceMemory<dfloat>  o_DX; 

   // *********************LIMITER RELATED*****************************//
   void stabSetupLimiter(); 
   void limSetupOperators(); 

   void limGeometricFactorsTri2D(); 
   void limGeometricFactorsQuad2D(); 
   void limGeometricFactorsTet3D(); 
   void limGeometricFactorsHex3D(); 









  /*******************************************/
  /*   ARTIFICIAL DIFFUSION STABILIZATION    */
  /*******************************************/
  // Memory for artificial Viscosity Activation Function 
  memory<dfloat> viscScale; 
  deviceMemory<dfloat> o_viscScale; 

  memory<dfloat> weight; 
  deviceMemory<dfloat> o_weight; 





  // Memory for artificial Viscosity Activation Function 
  // just in case it is needed!!!!
  memory<dfloat> viscActivation; 
  deviceMemory<dfloat> o_viscActivation; 

  memory<dfloat> vertexVisc; 
  deviceMemory<dfloat> o_vertexVisc; 

  memory<dfloat> visc; 
  deviceMemory<dfloat> o_visc; 

  memory<dfloat> projectVisc; 
  deviceMemory<dfloat> o_projectVisc; 

  // smmoth out viscosity
  ogs::ogs_t ogs;
  // memory<dfloat> vweight; 
  // deviceMemory<dfloat> o_vweight; 


  // Vertex Viscosity  
  kernel_t computeViscosityKernel; 
  kernel_t projectViscosityKernel; 

  /*******************************************/
  /*         SUBCELL STABILIZATION           */
  /*******************************************/
  int N, Nsubcells, Nint, Next, Nfields; 
  int Nverts, Nfaces, NfaceVertices, Np; 
  int Nvgeo, Nsgeo; 

  // minor girid connectivity
  memory<int> mEToV, mEToE, mEToF, mFaceNodes; 
  memory<int> faceVertices; 
  memory<dfloat> vr, vs, vt; // local vertex coordinates @ reference element
  memory<dfloat> cr, cs, ct; // center points of subcells @ reference element 
  memory<dfloat> fr, fs, ft; // center points of subcells @ reference element 
  memory<dfloat> mJ; // Jacobian of minor grid



  memory<int> mFToE, mFToF, mDGID; 
  deviceMemory<int> o_mFToE, o_mFToF, o_mDGID, o_mEToV, o_mFaceNodes;   

  // local projection 
  memory<dfloat> PM,  RM, PVM; // volume reconstuction, projection,  vertex  projection
  deviceMemory<dfloat> o_PM, o_RM, o_PVM;  

  memory<dfloat> PFM, RFM, SLIFT; // face projection and reconstruction 
  deviceMemory<dfloat> o_PFM, o_RFM, o_SLIFT; 

  // Connectivity
  memory<int> ielist, eelist; // local connectivity 
  deviceMemory<int> o_ielist, o_eelist; 

  // Glocal Connectivity
  memory<dlong> emapP, fmapP; // global connectivity
  deviceMemory<dlong> o_emapP, o_fmapP; // global connectivity

  // geometric info
  memory<dfloat> vgeo, sgeo; 
  deviceMemory<dfloat> o_vgeo, o_sgeo; 

  deviceMemory<dlong> o_EToE; 

  kernel_t findNeighKernel; 

  stab_t() = default;
  stab_t(platform_t& _platform, mesh_t &_mesh, stabSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  // stab setup
  void Setup(platform_t& _platform, mesh_t &_mesh, stabSettings_t& _settings);

  void detectSetup(){
    switch (detectorType) {
      case Stab::KLOCKNER:
         detectSetupKlockner(); 
        break;
      case Stab::PERSSON:
         detectSetupPersson(); 
        break;
    } 
  }


  void detectApply(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
    switch (detectorType) {
      case Stab::KLOCKNER:
         detectApplyKlockner(o_Q, o_RHS, T); 
        break;
      case Stab::PERSSON:
         detectApplyPersson(o_Q, o_RHS, T); 
        break;
    } 
  }

  // It is safer o split solver implementations I guess
  void stabSetup(){
    switch (stabType) {
      case Stab::FILTER:
        stabSetupFilter();
        break;
      case Stab::LIMITER:
         stabSetupLimiter();
         // LIBP_FORCE_ABORT("Limiter is not implemented yet");
        break;
      case Stab::ARTDIFF:
        stabSetupArtdiff();
        break;
      case Stab::SUBCELL:
        stabSetupSubcell();
        // LIBP_FORCE_ABORT("FV-Subcell is not implemented yet");
        break;
    } 
  }


   //***********************GENERAL IO OPERATIONS****************************//
   void Report(dfloat time, int tstep); 
   dlong GetElementNumber(deviceMemory<dlong>& eList); 
   void PlotElements(memory<dlong> eList, const std::string fileName);
   void PlotFields(memory<dfloat> Q, const std::string fileName);

   // *******************Detector Related Functions***************************//
   void detectSetupKlockner(); 
   void detectApplyKlockner(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void detectSetupPersson(); 
   void detectApplyPersson(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

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


   // *********************FILTER RELATED*****************************************//
   void stabSetupFilter(); 
   void stabApplyFilter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 

   void FilterMatrixTri2D (int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);
   void FilterMatrix1D(int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);
   void FilterMatrixTet3D (int _N, int _Nc, int _s, memory<dfloat>& _filterMatrix);




   // *********************ARTIFICIAL DIFFUSION RELATED*****************************//
   void stabSetupArtdiff(); 
   void stabApplyArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T); 
    
   dfloat ElementViscosityScaleTri2D(dlong e);
   dfloat ElementViscosityScaleQuad2D(dlong e);
   dfloat ElementViscosityScaleTet3D(dlong e);
   dfloat ElementViscosityScaleHex3D(dlong e);


   // *********************SUBCELL RELATED*****************************//
   void stabSetupSubcell(); 
   void stabSetupSubcellTri2D(); 


   void CellLocalConnectTri2D(); 
   void CellGlobalConnectTri2D(); 
   void CellCreateMinorGridTri2D();
   void CellGeometricFactorsTri2D(); 
   void CellSetupOperatorsTri2D(); 
   void CellFindBestMatchTri2D(dfloat x1, dfloat y1, dlong eP, int fP, 
                               memory<int> &elist, memory<dfloat> &x2,  memory<dfloat> &y2, 
                               int &nE, int &nF);


   void CellEquispacedEToVTri2D(const int _N, memory<int>& mEToV);
   void CellWarpBlendEToVTri2D(const int _N, memory<int>& mEToV);

  



 
 private:
  /*Set the type of mesh*/
  void setTypes(const Stab::SolverType sType, 
                const Stab::DetectorType dType, 
                const Stab::StabType stType);

 
};

} //namespace libp

#endif

