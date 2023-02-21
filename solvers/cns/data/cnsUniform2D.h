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

//Mean flow
#define p_RBAR 1.0
#define p_UBAR 0.990268069
#define p_VBAR 0.139173101

//P = RHO * U^2/(Ma^2*gamma)
// define p_PBAR 4.4642857
#define p_PBAR 0.496031746031746

// Ma=0.1
// #define p_PBAR 71.4286
// // Ma=0.3
// #define p_PBAR 7.936507936507937
// Ma=0.4
// #define p_PBAR 4.4642857
// // Ma=0.5
// #define p_PBAR 2.857142857142857
// // Ma=0.6
// #define p_PBAR 1.984126984126984
// // Ma=0.7
// #define p_PBAR 1.457725947521866
// Ma=0.8
// #define p_PBAR 1.116071428571428
// Ma=0.9
// #define p_PBAR 0.881834215167549
// // Ma=1.0
// #define p_PBAR 0.714285714285714
// // Ma=1.1
// #define p_PBAR 0.590318772136954
// Ma=1.2
// #define p_PBAR 0.496031746031746
// Ma=1.3
// #define p_PBAR 0.422654268808115
// Ma=1.4
// #define p_PBAR 0.364431486880467
// Ma=1.5
//#define p_PBAR 0.317460317460317


// Initial conditions (p is ignored for isothermal)
#define cnsInitialConditions2D(gamma, mu, t, x, y, r, u, v, p) \
{                                         \
  *(r) = p_RBAR;           \
  *(u) = p_UBAR;           \
  *(v) = p_VBAR;           \
  *(p) = p_PBAR;           \
}

// Body force
#define cnsBodyForce2D(gamma, mu, t, x, y, r, u, v, p, fx, fy) \
{                                                   \
  *(fx) = 0.0;                                      \
  *(fy) = 0.0;                                      \
}


// Boundary conditions for double check with FLEXI
// Adiabatic Wall, Weak Dirichlet for Inflow / Outflow
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
#define cnsBoundaryConditions2D(bc, gamma, R, CP, CV, mu, \
                                  t, x, y, nx, ny, \
                                  rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
                                  rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
{                                                                     \
const dfloat gammaM1  = gamma - 1.0;                                      \
const dfloat gammaP1  = gamma + 1.0;                                      \
const dfloat igammaM1 = 1.0/ gammaM1;                                     \
const dfloat igammaP1 = 1.0/ gammaP1;                                     \
  if(bc==1){                                                              \
    const dfloat uin  = uM*nx + vM*ny;                                    \
    const dfloat cin  = sqrt(gamma*pM/rM);                                \
    const dfloat min  = uin/cin;                                          \
    const dfloat pfunc= 2.0*gamma*igammaM1;                               \
    const dfloat AR   = 2.0*igammaP1/rM;                                  \
    const dfloat BR   = gammaM1*igammaP1*pM;                              \
    const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*gammaM1*min), pfunc);  \
    const dfloat PR2  = pM+0.5*uin/AR*(uin+sqrt(uin*uin+4.0*AR*(pM+BR))); \
    const dfloat PR   = (uin<=0) ? PR1 : PR2; \
    const dfloat TB   = pM/(rM*R) ;           \
    *(rB) = PR/(TB*R);                        \
    *(uB) = 0.0;                       \
    *(vB) = 0.0;                       \
    *(pB) = PR;                        \
    *(uxB) = uxM;                      \
    *(uyB) = uyM;                      \
    *(vxB) = vxM;                      \
    *(vyB) = vyM;                      \
  } else if(bc==2){                    \
    *(rB) = p_RBAR;                    \
    *(uB) = p_UBAR;                    \
    *(vB) = p_VBAR;                    \
    *(pB) = p_PBAR;                    \
    *(uxB) = uxM;                      \
    *(uyB) = uyM;                      \
    *(vxB) = vxM;                      \
    *(vyB) = vyM;                      \
  } else if(bc==3){                    \
    *(rB) = rM;                    \
    *(uB) = uM;                    \
    *(vB) = vM;                        \
    *(pB) = p_PBAR;                    \
    *(uxB) = 0.0;                      \
    *(uyB) = 0.0;                      \
    *(vxB) = 0.0;                      \
    *(vyB) = 0.0;                      \
  } else if(bc==4||bc==5){             \
    *(rB) = rM;                        \
    *(uB) = uM - (nx*uM+ny*vM)*nx;     \
    *(vB) = vM - (nx*uM+ny*vM)*ny;     \
    *(pB) = pM;                        \
    *(uxB) = uxM;                      \
    *(uyB) = uyM;                      \
    *(vxB) = vxM;                      \
    *(vyB) = vyM;                      \
  }                                    \
}

//**Pressure Outflow **//
// const dfloat c2 = gamma*pM/rM;         \
// const dfloat v2 = uM*uM + vM*vM;       \
// const dfloat pb = (v2 < c2) ? p_PBAR : pM; \             
//     *(rB) = gamma*pb/c2;                \
//     *(uB) = uM;                          \
//     *(vB) = vM;                         \
//     *(pB) = pb;                        \
//     *(uxB) = 0.0;                      \
//     *(uyB) = 0.0;                      \
//     *(vxB) = 0.0;                      \
//     *(vyB) = 0.0;                      \


//** Weak Drichlet Inflow or Outflow **//
//     *(rB) = p_RBAR;                    \
//     *(uB) = p_UBAR;                    \
//     *(vB) = p_VBAR;                    \
//     *(pB) = p_PBAR;                    \
//     *(uxB) = 0.0;                      \
//     *(uyB) = 0.0;                      \
//     *(vxB) = 0.0;                      \
//     *(vyB) = 0.0;                      \



//** Adiabatic / Isothermal Wall **//
// const dfloat gammaM1  = gamma - 1.0;                                      \
// const dfloat gammaP1  = gamma + 1.0;                                      \
// const dfloat igammaM1 = 1.0/ gammaM1;                                     \
// const dfloat igammaP1 = 1.0/ gammaP1;                                     \
//   if(bc==1){                                                              \
//     const dfloat uin  = uM*nx + vM*ny;                                    \
//     const dfloat cin  = sqrt(gamma*pM/rM);                                \
//     const dfloat min  = uin/cin;                                          \
//     const dfloat pfunc= 2.0*gamma*igammaM1;                               \
//     const dfloat AR   = 2.0*igammaP1/rM;                                  \
//     const dfloat BR   = gammaM1*igammaP1*pM;                              \
//     const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*gammaM1*min), pfunc);  \
//     const dfloat PR2  = pM+0.5*uin/AR*(uin+sqrt(uin*uin+4.0*AR*(pM+BR))); \
//     const dfloat PR   = (uin<=0) ? PR1 : PR2; \
//     const dfloat TB   = pM/(rM*R) (If Adiabatic);           \
//     const dfloat TB   = TW (if isothermal) ;           \
//     *(rB) = PR/(TB*R);                        \
//     *(uB) = 0.0;                       \
//     *(vB) = 0.0;                       \
//     *(pB) = PR;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \



//**Subsonic Outflow**//
//**if supersonic use total pressure to compute density**//
// const dfloat c2 = gamma*pM/rM;         \
// const dfloat v2 = uM*uM + vM*vM;       \
// const dfloat un = uM*nx + vM*ny;       \
// const dfloat pb = (v2 >=c2) ? (pM + 0.5*rM*v2) : p_PBAR; \          
// const dfloat uo = (un >=0.0) ? un : sqrt(v2); \          
//     *(rB) = gamma*pb/c2;               \
//     *(uB) = un*nx;                     \
//     *(vB) = un*ny;                     \
//     *(pB) = p_PBAR;                    \
//     *(uxB) = 0.0;                      \
//     *(uyB) = 0.0;                      \
//     *(vxB) = 0.0;                      \
//     *(vyB) = 0.0;                      \


//**Supersonic Outflow**//
//**if supersonic use total pressure to compute density**//       
//     *(rB) = rM;               \
//     *(uB) = uM;                     \
//     *(vB) = vM;                     \
//     *(pB) = pM;                    \
//     *(uxB) = 0.0;                      \
//     *(uyB) = 0.0;                      \
//     *(vxB) = 0.0;                      \
//     *(vyB) = 0.0;                      \




// // Boundary conditions for Subsonic case
// /* Using Adiabatic Wall PR/(TB*R); */     
// /* Pressure is set using Riemann Solver */
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// #define cnsBoundaryConditions2D(bc, gamma, mu, \
//                                   t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                                                     \
// const dfloat bulkMach = sqrt(p_UBAR*p_UBAR + p_VBAR*p_VBAR)/sqrt(gamma*p_PBAR/p_RBAR); \
// const dfloat R        = 1.0/(gamma*bulkMach*bulkMach);                    \
// const dfloat gammaM1  = gamma - 1.0;                                      \
// const dfloat gammaP1  = gamma + 1.0;                                      \
// const dfloat igammaM1 = 1.0/ gammaM1;                                     \
// const dfloat igammaP1 = 1.0/ gammaP1;                                     \
// const dfloat CP       = R*gamma*igammaM1;                                 \
// const dfloat CV       = R*igammaM1;                                       \
//   if(bc==1){                                                              \
//     const dfloat uin  = uM*nx + vM*ny;                                    \
//     const dfloat cin  = sqrt(gamma*pM/rM);                                \
//     const dfloat min  = uin/cin;                                          \
//     const dfloat pfunc= 2.0*gamma*igammaM1;                               \
//     const dfloat AR   = 2.0*igammaP1/rM;                                  \
//     const dfloat BR   = gammaM1*igammaP1*pM;                              \
//     const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*gammaM1*min), pfunc);  \
//     const dfloat PR2  = pM+0.5*uin/AR*(uin+sqrt(uin*uin+4.0*AR*(pM+BR))); \
//     const dfloat PR   = (uin<=0) ? PR1 : PR2; \
//     const dfloat TB   = pM/(rM*R) ;           \
//     *(rB) = PR/(TB*R);                        \
//     *(uB) = 0.0;                       \
//     *(vB) = 0.0;                       \
//     *(pB) = PR;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==2){                    \
//     *(rB) = p_RBAR;                    \
//     *(uB) = p_UBAR;                    \
//     *(vB) = p_VBAR;                    \
//     *(pB) = pM;                    \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==3){                    \
//   const dfloat c2   = gamma*pM/rM;       \
//   const dfloat vm2  = uM*uM +vM*vM;      \
//   const dfloat PB = (vm2 < c2) ? p_PBAR : pM; \
//     *(rB) = gamma*PB/c2;                \
//     *(uB) = uM;                        \
//     *(vB) = vM;                        \
//     *(pB) = PB;                    \
//     *(uxB) = 0.0;                      \
//     *(uyB) = 0.0;                      \
//     *(vxB) = 0.0;                      \
//     *(vyB) = 0.0;                      \
//   } else if(bc==4||bc==5){             \
//     *(rB) = rM;                        \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;     \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;     \
//     *(pB) = pM;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   }                                    \
// }




    // const dfloat uin  = uM*nx + vM*ny;                                    \
    // const dfloat cin  = sqrt(gamma*pM/rM);                                \
    // const dfloat min  = uin/cin;                                          \


// // Boundary conditions using Riemann Type
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// #define cnsBoundaryConditions2D(bc, gamma, R, CP, CV, mu, \
//                                   t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \                                                                   \
// {                                             \
// const dfloat gammaM1  = gamma - 1.0;                                      \
// const dfloat gammaP1  = gamma + 1.0;                                      \
// const dfloat igammaM1 = 1.0/ gammaM1;                                     \
// const dfloat igammaP1 = 1.0/ gammaP1;                                     \
// const dfloat u0n = p_UBAR*nx + p_VBAR*ny;     \
// const dfloat c0n = sqrt(gamma*p_PBAR/p_RBAR); \
// const dfloat uin = uM*nx + vM*ny;             \
// const dfloat cin = sqrt(gamma*pM/rM);         \
// const dfloat min = sqrt(uM*uM + vM*vM)/cin;   \
//   if(bc==1){                                                              \
//     const dfloat pfunc= 2.0*gamma*igammaM1;                               \
//     const dfloat AR   = 2.0*igammaP1/rM;                                  \
//     const dfloat BR   = gammaM1*igammaP1*pM;                              \
//     const dfloat PR1  = pM*pow(max(0.0001, 1.0+0.5*gammaM1*min), pfunc);  \
//     const dfloat PR2  = pM+0.5*uin/AR*(uin+sqrt(uin*uin+4.0*AR*(pM+BR))); \
//     const dfloat PR   = (uin<=0) ? PR1 : PR2; \
//     const dfloat TB   = pM/(rM*R) ;           \
//     *(rB) = PR/(TB*R);                        \
//     *(uB) = 0.0;                       \
//     *(vB) = 0.0;                       \
//     *(pB) = PR;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==2 || bc==3){                    \
//     if(uin<=0.0){                              \
//       if(min >=1.0 ){                           \
//         *(rB) = p_RBAR;                         \
//         *(uB) = p_UBAR;                         \
//         *(vB) = p_VBAR;                         \
//         *(pB) = p_PBAR;                         \
//         *(uxB) = uxM;                           \
//         *(uyB) = uyM;                           \
//         *(vxB) = vxM;                           \
//         *(vyB) = vyM;                           \
//       }else{                                    \
//         const dfloat RM  = u0n - 2.0*c0n*CV;                \
//         const dfloat RP  = uin + 2.0*cin*CV;                \
//         const dfloat vb  = 0.5*(RP + RM);                         \
//         const dfloat cR  = 0.25*gammaM1*(RP-RM);                  \
//         const dfloat uR  = p_UBAR + (vb - u0n)*nx;             \
//         const dfloat vR  = p_VBAR + (vb - u0n)*ny;             \
//         const dfloat sR  = c0n*c0n*pow(p_RBAR,-gammaM1)/gamma;   \
//         const dfloat rR  = pow(cR*cR/(gamma*sR), CV);      \
//         const dfloat pR   = rR * cR * cR / gamma;                 \
//         *(rB) = rR ;                            \
//         *(uB) = uR ;                            \
//         *(vB) = vR ;                            \
//         *(pB) = pR;                             \
//         *(uxB) = uxM;                           \
//         *(uyB) = uyM;                           \
//         *(vxB) = vxM;                           \
//         *(vyB) = vyM;                           \
//       }                                         \
//     }else{                                 \
//       if(min  >= 1.0 ){                     \
//         *(rB) = rM;                        \
//         *(uB) = uM;                        \
//         *(vB) = vM;                        \
//         *(pB) = pM;                        \
//         *(uxB) = 0.0;                      \
//         *(uyB) = 0.0;                      \
//         *(vxB) = 0.0;                      \
//         *(vyB) = 0.0;                      \
//        }else{                              \
//         const dfloat RM  = u0n - 2.0*c0n*CV;                \
//         const dfloat RP  = uin + 2.0*cin*CV;                \
//         const dfloat vb  = 0.5*(RP + RM);                         \
//         const dfloat cR  = 0.25*gammaM1*(RP-RM);                  \
//         const dfloat uR  = uM + (vb - uin)*nx;             \
//         const dfloat vR  = vM + (vb - uin)*ny;             \
//         const dfloat sR  = cin*cin*pow(rM,-gammaM1)/gamma;   \
//         const dfloat rR  =  pow(cR*cR/(gamma*sR), CV);      \
//         const dfloat pR  = rR * cR * cR / gamma;                 \
//         *(rB) = rR ;                       \
//         *(uB) = uR ;                       \
//         *(vB) = vR ;                       \
//         *(pB) = pR;                                               \
//         *(uxB) = 0.0;                           \
//         *(uyB) = 0.0;                           \
//         *(vxB) = 0.0;                           \
//         *(vyB) = 0.0;                           \
//        }                                                          \
//      }                                                            \
//   }else if(bc==4||bc==5){             \
//     *(rB) = rM;                        \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;     \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;     \
//     *(pB) = pM;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   }                                    \
// }







// // // Boundary conditions for Supersonic case
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// #define cnsBoundaryConditions2D(bc, gamma, mu, \
//                                   t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                      \
// const dfloat umn = uM*nx + vM*ny;      \
// const dfloat amn = sqrt(gamma*pM/rM);  \
// const dfloat man = sqrt(uM*uM + vM*vM)/amn;  \
//   if(bc==1){                           \
//     *(rB) = p_RBAR;                    \
//     *(uB) = 0.f;                       \
//     *(vB) = 0.f;                       \
//     *(pB) = p_PBAR;                    \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==2){                    \
//     *(rB) = p_RBAR;                       \
//     *(uB) = p_UBAR;                       \
//     *(vB) = p_VBAR;                       \
//     *(pB) = p_PBAR;                       \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==3){                    \
//     *(rB) = rM;                        \
//     *(uB) = uM;                        \
//     *(vB) = vM;                        \
//     *(pB) = pM;                        \
//     *(uxB) = 0.0;                      \
//     *(uyB) = 0.0;                      \
//     *(vxB) = 0.0;                      \
//     *(vyB) = 0.0;                      \
//   } else if(bc==4||bc==5){             \
//     *(rB) = rM;                        \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;     \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;     \
//     *(pB) = pM;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   }                                    \
// }







// dfloat RP=0., RM = 0.;                         \
// if(uin <=  0){                                  \
//   if(min>=1.0){                                  \
//     RP = u0n + 2.0*c0n*igammaM1;                \
//   }else{                                        \
//     RP = uin + 2.0*cin*igammaM1;                \
//   }                                             \
// }else{                                           \
//   RP = uin + 2.0*cin*igammaM1;                  \
//   if(min>=1.0){                                 \
//     RM = uin - 2.0*cin*igammaM1;                \
//   }else{                                        \
//     RM = u0n - 2.0*c0n*igammaM1;                \
//   }                                             \
} 

// // Boundary conditions
// // const dfloat m0n = sqrt(p_UBAR*p_UBAR + p_VBAR*p_VBAR)/c0n;  \
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// #define cnsBoundaryConditions2D(bc, gamma, mu, \
//                                   t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                      \
// const dfloat gammaM1 = gamma -1.0;            \
// const dfloat igammaM1 = 1.0/ gammaM1;          \
// const dfloat uin = uM*nx + vM*ny;              \
// const dfloat cin = sqrt(gamma*pM/rM);           \
// const dfloat min = sqrt(uM*uM + vM*vM)/cin;     \
// const dfloat u0n = p_UBAR*nx + p_VBAR*ny;       \
// const dfloat c0n = sqrt(gamma*p_PBAR/p_RBAR);   \
// dfloat RP  = uin + 2.0*cin*igammaM1;            \
// dfloat RM  = u0n - 2.0*c0n*igammaM1;           \
// RM         = (min >= 1.0 && uin>=0) ? (uin - 2.0*cin*igammaM1) : RM; \
// RP         = (min >= 1.0 && uin <0) ? (u0n + 2.0*c0n*igammaM1) : RP; \                                              \
// const dfloat cb  = 0.25*gammaM1*(RP-RM);                \
// const dfloat vb  = 0.5*(RP + RM);                       \
// const dfloat sin = cin*cin/(gamma*pow(rM,gammaM1));     \
// const dfloat s0n = c0n*c0n/(gamma*pow(p_RBAR,gammaM1)); \
// const dfloat sb =  (uin >=0.0) ? sin : s0n; \
// const dfloat rb =  pow(cb*cb/(gamma*sb), igammaM1); \
// const dfloat pb =  rb*cb*cb/gamma; \
//   if(bc==1){                           \
//     *(rB) = p_RBAR;                    \
//     *(uB) = 0.f;                       \
//     *(vB) = 0.f;                       \
//     *(pB) = p_PBAR;                    \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==2 || bc==3){           \
//     if(uin <= 0.0){                     \
//         *(rB) = rb;                     \
//         *(uB) = p_UBAR + (vb - u0n)*nx;                      \
//         *(vB) = p_VBAR + (vb - u0n)*ny;                      \
//         *(pB) = pb;                    \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//     }else{                                 \
//         *(rB) = rb;                \
//         *(uB) = uM + (vb - uin)*nx;                      \
//         *(vB) = vM + (vb - uin)*ny;                      \
//         *(pB) = pb;                    \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//      }                                     \
//   }else if(bc==4||bc==5){             \
//     *(rB) = rM;                        \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;     \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;     \
//     *(pB) = pM;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   }                                    \
// }

// #define cnsBoundaryConditions2D(bc, gamma, mu, \
//                                   t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                      \
// const dfloat uin = uM*nx + vM*ny;      \
// const dfloat cin = sqrt(gamma*pM/rM);  \
// const dfloat min = sqrt(uM*uM + vM*vM)/cin;  \
// const dfloat u0n = p_UBAR*nx + p_VBAR*ny;      \
// const dfloat c0n = sqrt(gamma*p_PBAR/p_RBAR);  \
// const dfloat m0n = sqrt(p_UBAR*p_UBAR + p_VBAR*p_VBAR)/c0n;  \
// const dfloat sin = cin*cin/(gamma*pow(rM,gamma-1.0));      \
// const dfloat s0n = c0n*c0n/(gamma*pow(p_RBAR,gamma-1.0)); \
// const dfloat sn =  uin> 0 ? sin : s0n; \
//   if(bc==1){                           \
//     *(rB) = p_RBAR;                    \
//     *(uB) = 0.f;                       \
//     *(vB) = 0.f;                       \
//     *(pB) = p_PBAR;                    \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==2 || bc==3){           \
//     if(umn < 0.0){                     \
//       if(man > 1.0 ){                     \
//         *(rB) = p_RBAR;                    \
//         *(uB) = p_UBAR;                    \
//         *(vB) = p_VBAR;                    \
//         *(pB) = p_PBAR;                    \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//       }else{                               \
//         *(rB) = p_RBAR ;                   \
//         *(uB) = p_UBAR ;                   \
//         *(vB) = p_VBAR ;                   \
//         *(pB) = pM;                        \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//       }                                    \
//     }else{                                 \
//       if(man  > 1.0 ){                     \
//         *(rB) = rM;                        \
//         *(uB) = uM;                        \
//         *(vB) = vM;                        \
//         *(pB) = pM;                        \
//         *(uxB) = 0.0;                      \
//         *(uyB) = 0.0;                      \
//         *(vxB) = 0.0;                      \
//         *(vyB) = 0.0;                      \
//        }else{                              \
//         *(rB) = rM ;                       \
//         *(uB) = uM ;                       \
//         *(vB) = vM ;                       \
//         *(pB) = p_PBAR;                                                   \
//         *(uxB) = 0.0;                                                     \
//         *(uyB) = 0.0;                                                     \
//         *(vxB) = 0.0;                                                     \
//         *(vyB) = 0.0;                                                     \
//        }                                                                  \
//      }                                                                    \
//   }else if(bc==4||bc==5){             \
//     *(rB) = rM;                        \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;     \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;     \
//     *(pB) = pM;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   }                                    \
// }

// // Boundary conditions
// /* wall 1, inflow 2, outflow 3, x-slip 4, y-slip 5 */
// #define cnsBoundaryConditions2D(bc, gamma, mu, \
//                                   t, x, y, nx, ny, \
//                                   rM, uM, vM, pM, uxM, uyM, vxM, vyM, \
//                                   rB, uB, vB, pB, uxB, uyB, vxB, vyB) \
// {                                      \
// const dfloat umn = uM*nx + vM*ny;      \
// const dfloat amn = sqrt(gamma*pM/rM);  \
// const dfloat man = sqrt(uM*uM + vM*vM)/amn;  \
//   if(bc==1){                           \
//     *(rB) = p_RBAR;                    \
//     *(uB) = 0.f;                       \
//     *(vB) = 0.f;                       \
//     *(pB) = p_PBAR;                    \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   } else if(bc==2 || bc==3){           \
//     if(umn < 0.0){                     \
//       if(man > 1.0 ){                     \
//         *(rB) = p_RBAR;                    \
//         *(uB) = p_UBAR;                    \
//         *(vB) = p_VBAR;                    \
//         *(pB) = p_PBAR;                    \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//       }else{                               \
//         const dfloat pR = 0.5*(pM + p_PBAR - rM*amn*(nx*(uM-p_UBAR) + ny*(vM - p_VBAR)));    \
//         *(rB) = p_RBAR + (pR-p_PBAR)/(amn*amn);                                             \
//         *(uB) = p_UBAR + nx*(pR - p_PBAR)/(amn*rM);                                         \
//         *(vB) = p_VBAR + ny*(pR - p_PBAR)/(amn*rM);                                         \
//         *(pB) = pR;                        \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//       }                                    \
//     }else{                                 \
//       if(man  > 1.0 ){                     \
//         *(rB) = rM;                        \
//         *(uB) = uM;                        \
//         *(vB) = vM;                        \
//         *(pB) = pM;                        \
//         *(uxB) = uxM;                      \
//         *(uyB) = uyM;                      \
//         *(vxB) = vxM;                      \
//         *(vyB) = vyM;                      \
//        }else{                                                             \
//         *(rB) = rM + (p_PBAR - pM) / (amn * amn);                         \
//         *(uB) = uM - nx* (p_PBAR - pM) / (amn * rM);                      \
//         *(vB) = vM - ny* (p_PBAR - pM) / (amn * rM);                      \
//         *(pB) = p_PBAR;                                                   \
//         *(uxB) = uxM;                                                     \
//         *(uyB) = uyM;                                                     \
//         *(vxB) = vxM;                                                     \
//         *(vyB) = vyM;                                                     \
//        }                                                                  \
//      }                                                                    \
//   }else if(bc==4||bc==5){             \
//     *(rB) = rM;                        \
//     *(uB) = uM - (nx*uM+ny*vM)*nx;     \
//     *(vB) = vM - (nx*uM+ny*vM)*ny;     \
//     *(pB) = pM;                        \
//     *(uxB) = uxM;                      \
//     *(uyB) = uyM;                      \
//     *(vxB) = vxM;                      \
//     *(vyB) = vyM;                      \
//   }                                    \
// }

