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

extern "C" __global__ void calcDeltaPeaks(int n, double m, double ax, double ay,
        float *finalField) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    int mx = nrX - 1;
    int my = nrY - 1;
    
    double out[3];
    calcSphericalOut(i, out);
    double bx = out[0];
    double by = out[1];
    
    double arg0;
    if (ax == bx) {
        arg0 = -((2*by*m*my*M_PI)/lambda);
    } else if (ay == by) {
        arg0 = -((2*bx*m*mx*M_PI)/lambda);
    } else if (ax == bx && ay == by) {
        finalField[2*i] = nrX*nrY;
        finalField[2*i+1] = 0;
        return;
    } else arg0 = -((2*m*(bx*mx+by*my)*M_PI)/lambda);
    
    double arg1 = (2*ax*m*(1+mx)*M_PI)/lambda;
    double arg2 = (2*bx*m*(1+mx)*M_PI)/lambda;
    double arg3 = (2*ay*m*(1+my)*M_PI)/lambda;
    double arg4 = (2*by*m*(1+my)*M_PI)/lambda;
    double arg5 = (2*ax*m*M_PI)/lambda;
    double arg6 = (2*bx*m*M_PI)/lambda;
    double arg7 = (2*ay*m*M_PI)/lambda;
    double arg8 = (2*by*m*M_PI)/lambda;
    
    cuDoubleComplex z0 = make_cuDoubleComplex(cos(arg0), sin(arg0));
    cuDoubleComplex z1;
    cuDoubleComplex z2;
    cuDoubleComplex z3;
    cuDoubleComplex z4;
    
    if (ax == bx) {
        z1 = make_cuDoubleComplex(nrX, 0);
        z2 = make_cuDoubleComplex(
                cos(arg3)-cos(arg4), sin(arg3)-sin(arg4));
        z3 = make_cuDoubleComplex(1, 0);
        z4 = make_cuDoubleComplex(
                cos(arg7)-cos(arg8), sin(arg7)-sin(arg8));
    } else if (ay == by) {
        z1 = make_cuDoubleComplex(
                cos(arg1)-cos(arg2), sin(arg1)-sin(arg2));
        z2 = make_cuDoubleComplex(nrY, 0);
        z3 = make_cuDoubleComplex(
                cos(arg5)-cos(arg6), sin(arg5)-sin(arg6));
        z4 = make_cuDoubleComplex(1, 0);
    } else {
        z1 = make_cuDoubleComplex(
                cos(arg1)-cos(arg2), sin(arg1)-sin(arg2));
        z2 = make_cuDoubleComplex(
                cos(arg3)-cos(arg4), sin(arg3)-sin(arg4));
        z3 = make_cuDoubleComplex(
                cos(arg5)-cos(arg6), sin(arg5)-sin(arg6));
        z4 = make_cuDoubleComplex(
                cos(arg7)-cos(arg8), sin(arg7)-sin(arg8));
    }
    
    z0 = cuCmul(z0, z1);
    z0 = cuCmul(z0, z2);
    z3 = cuCmul(z3, z4);
    z0 = cuCdiv(z0, z3);
    
    finalField[2*i] = static_cast<float>(z0.x);
    finalField[2*i + 1] = static_cast<float>(z0.y);
}