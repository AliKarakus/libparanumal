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


/// Level-Set function
#define lssInitialConditions2D(t, x, y, q) \
{    \
  const dfloat x1c  =  0.5,  y1c =  0.5, r1  = 0.25; \
  const dfloat x2c  = -0.5,  y2c = -0.5, r2  = 0.25; \
  const dfloat x3c  = -0.5,  y3c =  0.5, r3  = 0.25; \
  const dfloat x4c  =  0.5,  y4c = -0.5, r4  = 0.25; \
  const dfloat x5c  =  1.0,  y5c =  0.0, r5  = 0.25; \
  const dfloat x6c  =  0.0,  y6c =  1.0, r6  = 0.25; \
  const dfloat x7c  = -1.0,  y7c =  0.0, r7  = 0.25; \
  const dfloat x8c  =  0.0,  y8c = -1.0, r8  = 0.25; \
  const dfloat x9c  =  1.0,  y9c =  1.0, r9  = 0.25; \
  const dfloat x10c =  1.0,  y10c = -1.0, r10  = 0.25; \
  const dfloat x11c = -1.0,  y11c = -1.0, r11  = 0.25; \
  const dfloat x12c = -1.0,  y12c =  1.0, r12  = 0.25; \
  const dfloat test1 =  sqrt((x-x1c)*(x-x1c) + (y-y1c)*(y-y1c)) - r1;\
  const dfloat test2 =  sqrt((x-x2c)*(x-x2c) + (y-y2c)*(y-y2c)) - r2;\
  const dfloat test3 =  sqrt((x-x3c)*(x-x3c) + (y-y3c)*(y-y3c)) - r3;\
  const dfloat test4 =  sqrt((x-x4c)*(x-x4c) + (y-y4c)*(y-y4c)) - r4;\
  const dfloat test5 =  sqrt((x-x5c)*(x-x5c) + (y-y5c)*(y-y5c)) - r5;\
  const dfloat test6 =  sqrt((x-x6c)*(x-x6c) + (y-y6c)*(y-y6c)) - r6;\
  const dfloat test7 =  sqrt((x-x7c)*(x-x7c) + (y-y7c)*(y-y7c)) - r7;\
  const dfloat test8 =  sqrt((x-x8c)*(x-x8c) + (y-y8c)*(y-y8c)) - r8;\
  const dfloat test9 =  sqrt((x-x9c)*(x-x9c) + (y-y9c)*(y-y9c)) - r9;\
  const dfloat test10 =  sqrt((x-x10c)*(x-x10c) + (y-y10c)*(y-y10c)) - r10;\
  const dfloat test11 =  sqrt((x-x11c)*(x-x11c) + (y-y11c)*(y-y11c)) - r11;\
  const dfloat test12 =  sqrt((x-x12c)*(x-x12c) + (y-y12c)*(y-y12c)) - r12;\
  const dfloat t1 = min(min(min(min(min(min(min(test1, test2),test3),test4),test5),test6),test7),test8); \
  *q = pow( (x-1.0)*(x-1.0) + (y-1.0)*(y-1.0) + 0.1, 1.0)*min(min(min(min(t1, test9), test10),test11),test12); \
}

  // const dfloat test2 =-2.0*(sqrt(x*x + y*y) -(0.30 - 0.075*sin(4*atan((y-0.80)/x)))); \
// Level-Set function
#define lssExactSolution2D(t, x, y, q) \
{                                       \
}

 
// LS Advective field
#define lssAdvectionField2D(t, x, y, q, u, v) \
{                                       \
}


// Boundary conditions
/* wall 1, outflow 2 */
#define lssDirichletConditions2D(bc, t, x, y, nx, ny, qM, qB) \
{                                       \
  if(bc==1){                            \
    *(qB) = qM;                        \
  } else if(bc==2){                     \
    *(qB) = qM;                         \
  }                                     \
}