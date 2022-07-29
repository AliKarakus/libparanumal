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

//settings for lss solver
lssSettings_t::lssSettings_t(comm_t _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/lssVortex2D.h",
             "Boundary and Initial conditions header");

  newSetting("SOLVER TYPE",
             "ADVECTION",
             "Level set solver type",
             {"ADVECTION", "REDISTANCE"});

  newSetting("ADVECTION TYPE",
             "CUBATURE",
             "Advection integration method",
             {"COLLOCATION", "CUBATURE"});

  newSetting("FLOW TYPE",
             "INCOMPRESSIBLE",
             "Advection flow type",
             {"INCOMPRESSIBLE", "COMPRESSIBLE"});

  newSetting("TIME INTEGRATOR",
             "LSERK4",
             "Time integration method",
             {"AB3", "DOPRI5", "LSERK4"});

  newSetting("CFL NUMBER",
             "1.0",
             "Multiplier for timestep stability bound");

  newSetting("TIME RECONSTRUCTION",
             "ENO2",
             "Time recontruction order",
             {"ENO2", "ENO3"});


/*  newSetting("STABILIZATION",
             "SUBCELL",
             "Stabilization method", 
             {"SUBCELL", "NONE"});*/

/*  newSetting("SUBCELL MINOR GRID",
             "EQUISPACED",
             "Minor Triangulation",
             {"EQUISPACED", "WARPBLEND"});*/

/*  newSetting("INDICATOR TYPE",
             "MDA",
             "Troubled Cell Indicator Type",
             {"MDA", "MDH"});

  newSetting("SUBCELL NUMBER",
             "0",
             "Number of elements per direction");*/

  newSetting("START TIME",
             "0",
             "Start time for time integration");

  newSetting("FINAL TIME",
             "8",
             "End time for time integration");

  newSetting("OUTPUT INTERVAL",
             ".1",
             "Time between printing output data");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FILE NAME",
             "vortex",
             "Output file name");
}

void lssSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "LSS Settings:\n\n";
    reportSetting("SOLVER TYPE");
   /* reportSetting("REDISTANCE SOLVER");*/
    reportSetting("ADVECTION TYPE");
    reportSetting("FLOW TYPE");
    reportSetting("DATA FILE");
    reportSetting("TIME RECONSTRUCTION");
    reportSetting("TIME INTEGRATOR");
    /*reportSetting("STABILIZATION");
    reportSetting("SUBCELL NUMBER");
    reportSetting("SUBCELL MINOR GRID");
    reportS*//*etting("INDICATOR TYPE");*/
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");
  }
}

void lssSettings_t::parseFromFile(platformSettings_t& platformSettings,
                                  meshSettings_t& meshSettings,
                                  const std::string filename) {
  //read all settings from file
  settings_t s(comm);
  s.readSettingsFromFile(filename);

    for(auto it = s.settings.begin(); it != s.settings.end(); ++it) {
    setting_t& set = it->second;
    const std::string name = set.getName();
    const std::string val = set.getVal<std::string>();
    if (platformSettings.hasSetting(name))
      platformSettings.changeSetting(name, val);
    else if (meshSettings.hasSetting(name))
      meshSettings.changeSetting(name, val);
    else if (hasSetting(name)) //self
      changeSetting(name, val);
    else  {
      LIBP_FORCE_ABORT("Unknown setting: [" << name << "] requested");
    }
  }
}