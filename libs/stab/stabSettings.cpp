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

#include "stab.hpp"
#include "parAdogs.hpp"

namespace libp {

stabSettings_t::stabSettings_t(comm_t _comm):
  settings_t(_comm) {

  newSetting("STABILIZATION SOLVER TYPE",
             "1",
             "Solver type for stabilization",
             {"1", "2", "3"});

  newSetting("DETECTOR TYPE",
             "1",
             "Detector type for stabilization",
             {"1", "2"});

  newSetting("STAB TYPE",
             "1",
             "Stabilization type: FILTER:1, LIMITER:2, ART_DIFF:3, SUBCELL:4",
             {"1", "2", "3", "4"});

  newSetting("FILTER CUTOFF",
             "1",
             "Exponential Filter Cutoff Order");

  newSetting("FILTER ORDER",
             "2",
             "Exponential Filter Order (must be even)");


  newSetting("STAB OUTPUT TO FILE",
             "TRUE",
             "Detector Output to File ",
             {"TRUE", "FALSE"});
 
  paradogs::AddSettings(*this);
}

void stabSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "Stabilization Settings:\n\n";
    
    reportSetting("STABILIZATION SOLVER TYPE");
    reportSetting("DETECTOR TYPE");
    reportSetting("STAB TYPE");
    reportSetting("STAB OUTPUT TO FILE");

    reportSetting("FILTER CUTOFF");
    reportSetting("FILTER ORDER");
    reportSetting("STAB OUTPUT TO FILE");

    // if (!compareSetting("MESH FILE","BOX")) {
    //   paradogs::ReportSettings(*this);
    // }
  }
}

} //namespace libp
