/*
 * This file is part of the coherent_dmd_cimulator.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "cudaUtils.h"
#include <cuComplex.h>

extern "C" __global__ void calcSingleMirror(int n, double m,
        double ax, double ay, double az, double alpha, float *finalField) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    
    double out[3];
    calcSphericalOut(i, out);
    double bx = out[0];
    double by = out[1];
    double bz = out[2];
    
    double s2 = sqrtf(2);
    double ca = __cosf(alpha);
    double sa = __sinf(alpha);

    double r0 = lambda*lambda;
    double r1 = M_PI*M_PI;
    double r2 = ax+ay-bx-by+(ax-ay-bx+by)*ca-s2*(az-bz)*sa;
    double r3 = -ax-ay+bx+by+(ax-ay-bx+by)*ca-s2*(az-bz)*sa;
    double r = r0/r1/r2/r3;

    //double arg0 = 0;
    double arg1 = (2*ax*m+2*ay*m-2*bx*m-2*by*m)*M_PI/lambda;
    double arg2 = (ax*m+ay*m-bx*m-by*m+(ax-ay-bx+by)*m*ca-s2*(az-bz)*m*sa)*M_PI/lambda;
    double arg3 = (ax*m+ay*m-bx*m-by*m-(ax-ay-bx+by)*m*ca+s2*(az-bz)*m*sa)*M_PI/lambda;

    double re0 = 1;
    double im0 = 0;
    double re1 = __cosf(arg1);
    double im1 = __sinf(arg1);
    double re2 = __cosf(arg2);
    double im2 = __sinf(arg2);
    double re3 = __cosf(arg3);
    double im3 = __sinf(arg3);

    double nx = 1/s2*sa;
    double ny = -nx;
    double nz = sqrtf(1-sa*sa);
    double intesityFactor = fabsf(ax*nx+ay*ny+az*nz);

    finalField[2*i] = static_cast<float>(intesityFactor * r * (re0 + re1 - re2 - re3));
    finalField[2*i + 1] = static_cast<float>(intesityFactor * r * (im0 + im1 - im2 - im3));
    //if (i == 994119) printf("GPU: %f %f %f %f %f \n", -ax-ay+bx+by, (ax-ay-bx+by)*ca, s2*(az-bz)*sa, r3, r);
}