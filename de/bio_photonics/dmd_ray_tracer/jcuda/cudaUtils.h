#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

__constant__ int constantInts[4];
#define nrX constantInts[0] // nr of mirrors in x
#define nrY constantInts[1] // nr of mirrors in y
#define pMax constantInts[2] // nr of out angles in phi/x
#define tMax constantInts[3] //nr of out angles in theta/y

#define nrM (nrX*nrY) // nr of mirrors

__constant__ float constantFloats[13];
#define thetaMin constantFloats[0] // theta start value in rad
#define thetaMax constantFloats[1] // theta end value in rad
#define phiMin constantFloats[2] // phi start value in rad
#define phiMax constantFloats[3] // phi end value in rad
#define outStepSize constantFloats[4] // out step size in rad
#define mirrorSize constantFloats[5] // edge lenth of a single mirror
#define mirrorGap constantFloats[6] // size of gap between single mirrors
#define tilt constantFloats[7] // tilt angle of mirror in degree (not in use)
#define lambda constantFloats[8] // wavelenth of in beam laser in microns
#define beamDiameter constantFloats[9] // diameter of the gaussian in beam (not in use)
#define refPosition &constantFloats[10] // pointer to the reference position vector, equals calcDmdPosition(int mx = 0, int my = 0, float* dmdPosition)

#define dmdWidth (nrX * (mirrorSize + mirrorGap) - mirrorGap) // width of the whole dmd
#define dmdHeight (nrY * (mirrorSize + mirrorGap) - mirrorGap) // height of the whole dmd

/**
 * calculates a single normalized spherical vector
 * @param phi phi/x angle
 * @param theta theta/x angle
 * @param *out pointer to the out vector
 */
__device__ __inline__ void calcSphericalOutVector(float phi, float theta, float *out) {
    float tp = __tanf(phi);
    float tt = __tanf(theta);
    out[2] = sqrtf(1.0 / (tp*tp + tt*tt + 1));
    out[0] = out[2]*tp;
    out[1] = out[2]*tt;
}

/**
 * calculates out angles in spherical coordinates
 * @param i index of out angle/pixel
 * @param *out pointer to the result vector
 */
__device__ __inline__ void calcSphericalOut(int i, float *out) {
    int th = i / pMax;
    int ph = i % pMax;
    float theta = thetaMin + th * outStepSize;
    float phi = phiMin + ph * outStepSize;
    calcSphericalOutVector(phi, theta, out);
}

__device__ __inline__ void calcSphericalOutVector(double phi, double theta, double *out) {
    double tp = tan(phi);
    double tt = tan(theta);
    out[2] = sqrt(1.0 / (tp*tp + tt*tt + 1));
    out[0] = out[2]*tp;
    out[1] = out[2]*tt;
}

__device__ __inline__ void calcSphericalOut(int i, double *out) {
    int th = i / pMax;
    int ph = i % pMax;
    double theta = thetaMin + th * outStepSize;
    double phi = phiMin + ph * outStepSize;
    calcSphericalOutVector(phi, theta, out);
}

/**
 * calculates a single normalized cartesian vector (not in use)
 * @param x x coordinate
 * @param y y coordinate
 * @param z z coordinate
 * @param *out pointer to the out vector
 */
__device__ __inline__ void calcCartesianOutAngles(float x, float y, float z, float *out) {
    float length = sqrt(x*x+y*y+z*z);
    out[0] = x/length;
    out[1] = y/length;
    out[2] = z/length;
}

/**
 * calculates out angles in cartesian coordinates
 * @param i index of out angle/pixel
 * @param *out pointer to the result vector
 */
__device__ void calcCartesianOut(int i, float *out) {
    int xStep = i % pMax;
    int yStep = i / pMax;
    float xStart = phiMin / M_PI * 180.0;
    float yStart = thetaMin / M_PI * 180.0;
    float stepSize = outStepSize /M_PI * 180.0;
    float x = xStart + xStep * stepSize;
    float y = yStart + yStep * stepSize;
    float z = 20.0;
    calcCartesianOutAngles(x, y, z, out);
}

/**
 * calculates the cartesian vector to the position of a single mirror on the dmd
 * @param mx index of the x position
 * @param my index of the y position
 * @param *dmdPosition pointer to the result vector
 */
__device__ __inline__ void calcDmdPosition(int mx, int my, float* dmdPosition) {
    dmdPosition[0] = (mirrorSize + mirrorGap) * mx - dmdWidth / 2;
    dmdPosition[1] = (mirrorSize + mirrorGap) * my - dmdHeight / 2;
    dmdPosition[2] = 0;
}

/**
 * calculates the path length between dmd position and the surface
 * which is perpendicular to the out vector
 * @param *out pointer to the out vector
 * @param *dmdPosition pointer to the vector of the position on the dmd 
 */
__device__ __inline__ float calcOutPathLength(float* out, float *dmdPosition) {
	float* referencePosition = refPosition;
    float referencePl = -(referencePosition[0] * out[0] +
            referencePosition[1] * out[1] + referencePosition[2] * out[2]);
    
    float currentPl = -(dmdPosition[0] * out[0] + dmdPosition[1] * out[1] +
            dmdPosition[2] * out[2]);
    return currentPl - referencePl;
}

/**
 * calculates the absolute value of a complex number
 * @param re real part
 * @param im imaginary part
 */
__device__ __inline__ float calcComplexAbs(float re, float im) {
    return sqrt(re*re + im*im);
}
