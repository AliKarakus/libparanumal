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

// Level-Set function Rotated Square
#define lssInitialConditions2D(t, x, y, q) \
{                                       \
  const dfloat xc = 0.875f;           \
  const dfloat yc = 0.5f;         \
  const dfloat A  = 1.0f;          \
  const dfloat B  = 0.5f;          \
  const dfloat scale = pow( (x-xc)*(x-xc) + (y-yc)*(y-yc) + 0.1, 1.0); \
  (*q) = scale*(sqrt(x*x/(A*A) + y*y/(B*B)) - 1.0) ; \
}

  // const dfloat yc = 2.0f/4.0f;         \

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