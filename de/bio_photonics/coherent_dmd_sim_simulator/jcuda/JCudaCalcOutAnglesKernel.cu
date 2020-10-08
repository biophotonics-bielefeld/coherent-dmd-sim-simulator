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

/** cuda kernel for the dmd ray tracer
 * calculates the field for all out angles
 * @param n number of out angles/pixels
 * @param mStart index of first mirror for this kernel
 * @param mEnd index of (last mirror + 1) for this kernel
 * @param *tiltStates pointer to the array of tiltStates (1&0 as true/false)
 * @param *mirrorTrue pointer to the array of the reference field for true mirrors
 * @param *mirrorFalse pointer to the array of the reference field for false mirrors
 * @param *inOffsetPathLengths pointer to the array which contains the pathlength
 * between the incoming wave and each mirror
 * @param *beamProfile pointer to the array of in beam intensitys for each mirror
 * @param *finalField pointer to the array for the resulting field for each out angle
 */

extern "C" __global__ void calcOutAngles(int n, int mStart, int mEnd, int *tiltStates,
		float *mirrorTrue, float *mirrorFalse, float *inOffsetPathLengths, float *beamProfile, float *finalField) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    
    float out[3];
    calcSphericalOut(i, out);
    //calcCartesianOut(i, out);
    
    for (int m = mStart; m < mEnd; m++) {
        float* referenceMirror = tiltStates[m] ? mirrorTrue : mirrorFalse;
        float referenceFieldRe = referenceMirror[2 * i];
        float referenceFieldIm = referenceMirror[2 * i + 1];
        float referenceFieldAbs = calcComplexAbs(referenceFieldRe, referenceFieldIm);
        float referenceFieldArg = atan2(referenceFieldIm, referenceFieldRe);
        
        int my = m / nrX;
        int mx = m % nrX;
        float dmdPosition[3];
        calcDmdPosition(mx, my, dmdPosition);
        
        float initInPathLength = inOffsetPathLengths[m];
        float outPathLength = calcOutPathLength(out, dmdPosition);
        float additionalPl = initInPathLength + outPathLength;
        float additionalPhase = (additionalPl / lambda) * 2 * M_PI;
        float r = referenceFieldAbs * beamProfile[m];
        float p = referenceFieldArg + additionalPhase;
        float finalFieldRe = r * cos(p);
        float finalFieldIm = r * sin(p);
        finalField[2*i] += finalFieldRe;
        finalField[2*i + 1] += finalFieldIm;
    }
}