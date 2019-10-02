/*
 * Copyright (C) 2019 m.lachetta
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
package de.bio_photonics.coherent_dmd_sim_simulator;

import de.bio_photonics.dmd_ray_tracer.jcuda.JCudaEngine;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;
import java.io.IOException;

/**
 * class for dmd simulations
 * @author m.lachetta
 */
public class DmdSimulator {
    static final JCudaEngine gpuEngine;
    static String gpuModuleName;
    static String gpuFunctionName;
    
    int nrX, nrY, pMax, tMax;
    double lambda, beamDiameter, mirrorSize, gap, thetaMinD, thetaMinR,
            thetaMaxD, thetaMaxR, phiMinD, phiMinR, phiMaxD, phiMaxR,
            outStepSizeD, outStepSizeR, tiltD;
    Dmd dmd;
    Vector inBeam;
    Vector referencePosition;
    
    boolean[][] tiltStates;
    Vector[][] dmdPositions;
    double[][] gaussians;
    double[][] inOffsetPathLengths;
    
    Vector[][] outAngles;
    Complex[][] finalField;
    Complex[][] currentField;
    Complex[][] mirrorTrue;
    Complex[][] mirrorFalse;
    
    int trueThStart;
    int trueThEnd;
    int truePhStart;
    int truePhEnd;
    int falseThStart;
    int falseThEnd;
    int falsePhStart;
    int falsePhEnd;
    
    boolean gpuActive;
    int jcPixelsPerKernel;
    JCudaEngine.IntArray jcTiltState;
    JCudaEngine.FloatArray jcGaussians;
    JCudaEngine.FloatArray jcField;
    
    static {
        JCudaEngine.setExceptionsEnabled(true);
        gpuEngine = new JCudaEngine();
        gpuModuleName = "JCudaCalcOutAnglesKernel";
        String cuFilePath ="de\\bio_photonics\\dmd_ray_tracer\\jcuda\\JCudaCalcOutAnglesKernel.cu";
        gpuFunctionName = "calcOutAngles";
        try {
            gpuEngine.loadModule(gpuModuleName, cuFilePath);
            gpuEngine.loadModule("JCudaCalcDeltaPeaks", "de\\bio_photonics\\dmd_ray_tracer\\jcuda\\JCudaCalcDeltaPeaksKernel.cu");
            gpuEngine.loadModule("JCudaCalcSingleMirror", "de\\bio_photonics\\dmd_ray_tracer\\jcuda\\JCudaCalcSingleMirrorKernel.cu");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        
        gpuEngine.loadFunktion(gpuModuleName, gpuFunctionName);
        gpuEngine.loadFunktion("JCudaCalcDeltaPeaks", "calcDeltaPeaks");
        gpuEngine.loadFunktion("JCudaCalcSingleMirror", "calcSingleMirror");
    }
    
    /**
     * class to handle meta data
     */
    static class MetaData {
        String outDir;
        int[] lambdas;
        boolean gpuActive;
        int beamDiameter;
        
        int nrX;
        int nrY;
        
        double latticeConstant;
        double fillFactor;
        
        double tiltAngle;
        
        double phiOutStart;
        double phiOutEnd;
        double thetaOutStart;
        double thetaOutEnd;
        double outStepSize;
        
        double phiInStart;
        double phiInEnd;
        double thetaInStart;
        double thetaInEnd;
        double inStepSize;
        
        Image bmp;
        
        @Override
        public String toString() {
            String lambdasString = "";
            for (int lambda : lambdas) lambdasString += (lambda + " ");
            String retString = "outDir: " + outDir +
                    "\ngpuActive: " + gpuActive +
                    "\nbeamDiameter: " + beamDiameter +
                    "\nlambdas: " + lambdasString +
                    "\nnrX: " + nrX +
                    "\nnrY: " + nrY +
                    "\nmirrorSize: " + latticeConstant +
                    "\nfillFactor: " + fillFactor +
                    "\ntiltAngle: " + tiltAngle +
                    "\nphiOutStart: " + phiOutStart +
                    "\nphiOutEnd: " + phiOutEnd +
                    "\nthetaOutStart: " + thetaOutStart +
                    "\nthetaOutEnd: " + thetaOutEnd +
                    "\noutStepSize: " + outStepSize +
                    "\nphiInStart: " + phiInStart +
                    "\nphiInEnd: " + phiInEnd +
                    "\nthetaInStart: " + thetaInStart +
                    "\nthetaInEnd: " + thetaInEnd +
                    "\ninStepSize: " + inStepSize;
            return retString;
        }
    }
    
    /**
     * 
     * @param dmd dmd object
     * @param tiltAngle tilt angle of the mirrors in degree
     * @param lambda wavelength of the incident beam in µm
     * @param beamDiameter in µm
     * @param phMin start out angle of x in degree
     * @param phMax end out angle of x in degree
     * @param thMin start out angle of y in degree
     * @param thMax end out angle of y in degree
     * @param stepSize step size of out angles
     * @param flipStates bitmap (1/0;true/false) image for flip states of the mirrors on the dmd
     * @param gpuActive running ray tracer with CUDA support?
     */
    DmdSimulator(Dmd dmd, double tiltAngle, double lambda,
            double beamDiameter, double phMin, double phMax, double thMin,
            double thMax, double stepSize, Image flipStates, boolean gpuActive) {
        this.dmd = dmd;
        this.nrX = dmd.nrX;
        this.nrY = dmd.nrY;
        this.mirrorSize = dmd.mirrorWidth;
        if (mirrorSize != dmd.mirrorHeight)
            throw new RuntimeException("Only squared mirrors are supported");
        this.gap = dmd.gapX;
        if (gap != dmd.gapY) throw new RuntimeException("Only equal x y gaps are ");
        
        this.tiltD = tiltAngle;
        
        this.lambda = lambda;
        this.beamDiameter = beamDiameter;
        this.phiMinD = phMin;
        this.phiMaxD = phMax;
        this.thetaMinD = thMin;
        this.thetaMaxD = thMax;
        this.outStepSizeD = stepSize;
        
        // switch from degree to rad
        thetaMinR = thetaMinD * Math.PI / 180;
        thetaMaxR = thetaMaxD * Math.PI / 180;
        phiMinR = phiMinD * Math.PI / 180;
        phiMaxR = phiMaxD * Math.PI / 180;
        outStepSizeR = outStepSizeD * Math.PI / 180;
        
        // calculateing range of for loops over out angles
        tMax = (int) ((thetaMaxD - thetaMinD) / outStepSizeD);
        pMax = (int) ((phiMaxD - phiMinD) / outStepSizeD);
        
        // position on the dmd for which the reference fields will be calculated
        referencePosition = dmd.getCoordinates(0, 0, tiltAngle, 0, 0);
        //referenceFalse = dmd.getCoordinates(0, 0, -tiltAngle, 0, 0);
        
        // pre calculations over the mirrors
        this.tiltStates = new boolean[nrY][nrX];
        dmdPositions = new Vector[nrY][nrX];
        gaussians = new double[nrY][nrX];
        inOffsetPathLengths = new double[nrY][nrX];
        for (int my = 0; my < nrY; my++) {
            for (int mx = 0; mx < nrX; mx++) {
                tiltStates[my][mx] = flipStates.get(mx, my) != 0;
                double angleOfTilt = tiltStates[my][mx] ? tiltD : -tiltD;
                dmdPositions[my][mx] = dmd.getCoordinates(mx, my, angleOfTilt, 0, 0);
                double x = dmd.getX(mx, my, tiltAngle, 0, 0);
                double y = dmd.getY(mx, my, tiltAngle, 0, 0);
                double r = Math.sqrt(x*x+y*y);
                gaussians[my][mx] = gaussian(r, beamDiameter);
            }
        }
        
        // pre calculations over the out angles
        mirrorTrue = new Complex[tMax][pMax];
        mirrorFalse = new Complex[tMax][pMax];
        finalField = new Complex[tMax][pMax];
        for (int th = 0; th < tMax; th++) {
            for (int ph = 0; ph < pMax; ph++) {
                mirrorTrue[th][ph] = new Complex(0, 0);
                mirrorFalse[th][ph] = new Complex(0, 0);
                finalField[th][ph] = new Complex(0, 0);
            }
        }
        outAngles = calcSphericalOutAngles(phiMinR, pMax, thetaMinR, tMax, outStepSizeR);
        //outAngles = calcCartesianOutAngles(phiMinD, pMax, thetaMinD, tMax, outStepSizeD, 12.56);
        
        // init CUDA GPU
        this.gpuActive = gpuActive;
        if (gpuActive) initCudaMemory();
        
        System.out.println("Finished constructor");
    }
    
    /**
     * copies data from CPU to GPU for CUDA
     */
    private void initCudaMemory() {
        int[] constantInts = new int[4];
        constantInts[0] = nrX;
        constantInts[1] = nrY;
        constantInts[2] = pMax;
        constantInts[3] = tMax;
        gpuEngine.writeConstant(gpuModuleName, "constantInts", constantInts);
        gpuEngine.writeConstant("JCudaCalcDeltaPeaks", "constantInts", constantInts);
        gpuEngine.writeConstant("JCudaCalcSingleMirror", "constantInts", constantInts);
        
        float[] constantFloats = new float[13];
        constantFloats[0] = (float) thetaMinR;
        constantFloats[1] = (float) thetaMaxR;
        constantFloats[2] = (float) phiMinR;
        constantFloats[3] = (float) phiMaxR;
        constantFloats[4] = (float) outStepSizeR;
        constantFloats[5] = (float) mirrorSize;
        constantFloats[6] = (float) gap;
        constantFloats[7] = (float) tiltD;
        constantFloats[8] = (float) lambda;
        constantFloats[9] = (float) beamDiameter;
        constantFloats[10] = (float) referencePosition.getX();
        constantFloats[11] = (float) referencePosition.getY();
        constantFloats[12] = (float) referencePosition.getZ();
        gpuEngine.writeConstant(gpuModuleName, "constantFloats", constantFloats);
        gpuEngine.writeConstant("JCudaCalcDeltaPeaks", "constantFloats", constantFloats);
        gpuEngine.writeConstant("JCudaCalcSingleMirror", "constantFloats", constantFloats);
        
        for (int y = 0; y < tiltStates.length; y++) {
            for (int x = 0; x < tiltStates[0].length; x++) {
            }
        }
        jcTiltState = gpuEngine.createIntArray(tiltStates);
        jcGaussians = gpuEngine.createFloatArray(gaussians);
        jcField = gpuEngine.createFloatArray(finalField);
        jcTiltState.hostToDevice();
        jcGaussians.hostToDevice();
        jcField.hostToDevice();
    }
    
    /**
     * creates spherical out vectors (not standard spherical)
     * @param phiMin start out angle of phi in degree
     * @param nrPhi amount of phi steps
     * @param thetaMin start out angle of phi in degree
     * @param nrTheta amount of theta steps
     * @param stepSize degree per step
     * @return 2D array of vectors
     */
    private static Vector[][] calcSphericalOutAngles(double phiMin, int nrPhi,
            double thetaMin, int nrTheta, double stepSize) {
        Vector[][] angles = new Vector[nrTheta][nrPhi];
        for (int th = 0; th < nrTheta; th++) {
            double theta = thetaMin + th * stepSize;
            for (int ph = 0; ph < nrPhi; ph++) {
                double phi = phiMin + ph * stepSize;
                angles[th][ph] = new Vector(phi, theta);
            }
        }
        return angles;
    }
    
    /**
     * creates normalized cartesian out vectors
     * @param xStart start x value in free to choose units
     * @param xSteps amount of x steps
     * @param yStart start y value in free to choose units
     * @param ySteps amount of y steps
     * @param stepSize free to choose unit per step
     * @param zDistance distance between dmd and screen
     * @return 2D array of vectors
     */
    private static Vector[][] calcCartesianOutAngles(double xStart, int xSteps,
            double yStart, int ySteps, double stepSize, double zDistance) {
        Vector[][] angles = new Vector[ySteps][xSteps];
        for (int y = 0; y < ySteps; y++) {
            double yy = yStart + y * stepSize;
            for (int x = 0; x < xSteps; x++) {
                double xx = xStart + x * stepSize;
                angles[y][x] = new Vector(xx, yy, zDistance);
                angles[y][x].normalize();
            }
        }
        return angles;
    }
    
    /**
     * sets a new incident beam angle for this ray tracer
     * @param inBeam normalized in beam vector
     */
    void setInBeam(Vector inBeam) {
        this.inBeam = inBeam;
        for (int my = 0; my < nrY; my++) {
            for (int mx = 0; mx < nrX; mx++) {
                double tiltState = tiltStates[my][mx] ? tiltD : -tiltD;
                inOffsetPathLengths[my][mx] = calcInOffsetPathLength(mx, my, tiltState);
            }
        }
    }
    
    /**
     * calculates the electric complex field for a specified out angle in an
     * analytical manner (look into dmd.nb for details)
     * @param out vector of the out angle
     * @param alpha tilt angle of the mirror in rad
     * @return the complex field of the out angle
     */
    Complex calcAnalyticSingleMirrorOutAngle(Vector out, double alpha) {
        double m = mirrorSize;
        double ax = inBeam.getX();
        double ay = inBeam.getY();
        double az = inBeam.getZ();
        double bx = out.getX();
        double by = out.getY();
        double bz = out.getZ();
        
        double s2 = Math.sqrt(2);
        double ca = Math.cos(alpha);
        double sa = Math.sin(alpha);
        
        double r0 = lambda*lambda;
        double r1 = Math.PI*Math.PI;
        double r2 = ax+ay-bx-by+(ax-ay-bx+by)*ca-s2*(az-bz)*sa;
        double r3 = -ax-ay+bx+by+(ax-ay-bx+by)*ca-s2*(az-bz)*sa;
        double r = r0/r1/r2/r3;
        
        double argFactor = Math.PI/lambda;
        double arg0 = 0;
        double arg1 = (2*ax*m+2*ay*m-2*bx*m-2*by*m)*argFactor;
        double arg2 = (ax*m+ay*m-bx*m-by*m+(ax-ay-bx+by)*m*ca-s2*(az-bz)*m*sa)*argFactor;
        double arg3 = (ax*m+ay*m-bx*m-by*m-(ax-ay-bx+by)*m*ca+s2*(az-bz)*m*sa)*argFactor;
        
        double re0 = 1;
        double im0 = 0;
        double re1 = Math.cos(arg1);
        double im1 = Math.sin(arg1);
        double re2 = Math.cos(arg2);
        double im2 = Math.sin(arg2);
        double re3 = Math.cos(arg3);
        double im3 = Math.sin(arg3);
        
        double nx = 1/s2*sa;
        double ny = -nx;
        double nz = Math.sqrt(1-sa*sa);
        Vector n = new Vector(nx,ny,nz);
        double intesityFactor = Math.abs(inBeam.times(n));
        
        double re = intesityFactor * r * (re0 + re1 - re2 - re3);
        double im = intesityFactor * r * (im0 + im1 - im2 - im3);
        
        return new Complex(re, im);
    }
    
    
    Complex[][] calcAllDeltaPeaks() {
        double pi = Math.PI;
        double m = mirrorSize+gap;
        double ax = inBeam.getX();
        double ay = inBeam.getY();
        int mx = nrX - 1;
        int my = nrY - 1;
        
        Complex[][] field = new Complex[tMax][pMax];
        for (int th = 0; th < tMax; th++) {
            double theta = thetaMinR + th * outStepSizeR;
            for (int ph = 0; ph < pMax; ph++) {
                double phi = phiMinR + ph * outStepSizeR;
                Vector out = new Vector(phi, theta);
                double bx = out.getX();
                double by = out.getY();
                
                double arg0 = -((2*m*(bx*mx+by*my)*pi)/lambda);
                double arg1 = (2*ax*m*(1+mx)*pi)/lambda;
                double arg2 = (2*bx*m*(1+mx)*pi)/lambda;
                double arg3 = (2*ay*m*(1+my)*pi)/lambda;
                double arg4 = (2*by*m*(1+my)*pi)/lambda;
                double arg5 = (2*ax*m*pi)/lambda;
                double arg6 = (2*bx*m*pi)/lambda;
                double arg7 = (2*ay*m*pi)/lambda;
                double arg8 = (2*by*m*pi)/lambda;
                
                if (ax == bx) {
                    arg0 = -((2*by*m*my*pi)/lambda);
                } else if (ay == by) {
                    arg0 = -((2*bx*m*mx*pi)/lambda);
                } else if (ax == bx && ay == by) {
                    field[th][ph] = new Complex(nrX*nrY, 0);
                    break;
                }
                
                double re0 = Math.cos(arg0);
                double im0 = Math.sin(arg0);
                double re1 = Math.cos(arg1);
                double im1 = Math.sin(arg1);
                double re2 = Math.cos(arg2);
                double im2 = Math.sin(arg2);
                double re3 = Math.cos(arg3);
                double im3 = Math.sin(arg3);
                double re4 = Math.cos(arg4);
                double im4 = Math.sin(arg4);
                double re5 = Math.cos(arg5);
                double im5 = Math.sin(arg5);
                double re6 = Math.cos(arg6);
                double im6 = Math.sin(arg6);
                double re7 = Math.cos(arg7);
                double im7 = Math.sin(arg7);
                double re8 = Math.cos(arg8);
                double im8 = Math.sin(arg8);
                
                Complex z0 = new Complex(re0, im0);
                Complex z1 = new Complex(re1-re2, im1-im2);
                Complex z2 = new Complex(re3-re4, im3-im4);
                Complex z3 = new Complex(re5-re6, im5-im6);
                Complex z4 = new Complex(re7-re8, im7-im8);
                
                if (ax == bx) {
                    z1 = new Complex(nrX, 0);
                    z3 = new Complex(1, 0);
                } else if (ay == by) {
                    z2 = new Complex(nrY, 0);
                    z4 = new Complex(1, 0);
                }
                
                z0.multi(z1);
                z0.multi(z2);
                z3.multi(z4);
                z0.divide(z3);
                
                field[th][ph] = z0;
            }
        }
        return field;
    }
    
    /**
     * calculates the complex fields for true and false tilt states for all out angles
     * @return intensity distribution for a true and false mirror
     */
    Image[] computeAnalyticSingleMirrors() {
        Image viewTrue = new Image(pMax, tMax);
        Image viewFalse = new Image(pMax, tMax);
        viewTrue.setTitle((int) (lambda*1000) + "-true");
        viewFalse.setTitle((int) (lambda*1000) + "-false");
        double alpha = tiltD / 180. * Math.PI;
        
        for (int th = 0; th < tMax; th++) {
            for (int ph = 0; ph < pMax; ph++) {
                mirrorTrue[th][ph] = calcAnalyticSingleMirrorOutAngle(outAngles[th][ph], alpha);
                mirrorFalse[th][ph] = calcAnalyticSingleMirrorOutAngle(outAngles[th][ph], -alpha);
                viewTrue.set(ph, th, (float) mirrorTrue[th][ph].square());
                viewFalse.set(ph, th, (float) mirrorFalse[th][ph].square());
            }
        }
        
        if (!gpuActive) setCpuInterestingArea(7);
        Image[] trueFalse = {viewTrue, viewFalse};
        return trueFalse;
    }
    
    /**
     * sets a region of interest around the maximum in the reference field for cpu calculations
     * @param rangeOfInterest in degree
     */
    protected void setCpuInterestingArea(int rangeOfInterest) {
        double trueMax = -1;
        double falseMax = -1;
        int trueMaxTh = -1, trueMaxPh = -1, falseMaxTh = -1, falseMaxPh = -1;
        for (int th = 0; th < tMax; th++) {
            for (int ph = 0; ph < pMax; ph++) {
                double t = mirrorTrue[th][ph].square();
                double f = mirrorFalse[th][ph].square();
                if (t > trueMax) {
                    trueMax = t;
                    trueMaxTh = th;
                    trueMaxPh = ph;
                }
                if (f > falseMax) {
                    falseMax = f;
                    falseMaxTh = th;
                    falseMaxPh = ph;
                }
            }
        }
        
        double range = rangeOfInterest / 180. * Math.PI;
        int steps = (int) (range / outStepSizeR) + 1;
        trueThStart = trueMaxTh - steps;
        trueThEnd = trueMaxTh + steps;
        truePhStart = trueMaxPh - steps;
        truePhEnd = trueMaxPh + steps;
        falseThStart = falseMaxTh - steps;
        falseThEnd = falseMaxTh + steps;
        falsePhStart = falseMaxPh - steps;
        falsePhEnd = falseMaxPh + steps;
        
        trueThStart = checkThRange(trueThStart);
        trueThEnd = checkThRange(trueThEnd);
        truePhStart = checkPhRange(truePhStart);
        truePhEnd = checkPhRange(truePhEnd);
        falseThStart = checkThRange(falseThStart);
        falseThEnd = checkThRange(falseThEnd);
        falsePhStart = checkPhRange(falsePhStart);
        falsePhEnd = checkPhRange(falsePhEnd);
    }
    
    protected int checkThRange(int thValue) {
        if (thValue < 0) return 0;
        else if (thValue > tMax) return tMax;
        else return thValue;
    }
    
    protected int checkPhRange(int phValue) {
        if (phValue < 0) return 0;
        else if (phValue > pMax) return pMax;
        else return phValue;
    }
    
    /**
     * calculates the offset b
     * @param mx
     * @param my
     * @param flipState
     * @return 
     */
    double calcInOffsetPathLength(int mx, int my, double flipState) {
        Vector currentPosition = dmd.getCoordinates(mx, my, flipState, 0, 0); // set reference postion to 0 0
        //double referencePl = referencePosition.times(inBeam);
        double currentPl = currentPosition.times(inBeam);
        return currentPl; // - referencePl
    }
    
    /**
     * calculates the path length between dmd and the surface which is
     * perpendicular to the out vector
     * @param mx index of mirror in x
     * @param my index of mirror in y
     * @param out out beam vector
     * @return out path length
     */
    double calcOutOffsetPathLength(int mx, int my, Vector out) {
        //Vector currentPosition = dmd.getCoordinates(mx, my, tiltAngles[my][mx], 0, 0);
        double referencePl = -referencePosition.times(out);
        double currentPl = -dmdPositions[my][mx].times(out);
        return currentPl - referencePl;
    }
    
    /**
     * calculates the field of a specific out angle and saves it into {@link #finalField}
     * @param ph step index of phi
     * @param th step index of theta
     */
    void calcCpuOutAngle(int ph, int th) {
        for (int my = 0; my < nrY; my++) {
            for (int mx = 0; mx < nrX; mx++) {
                //double tiltAngle = tiltAngles[my][mx];
                //Vector referencePosition = tiltStates[my][mx] ? referenceTrue : referenceFalse;
                Complex[][] referenceMirror = tiltStates[my][mx] ? mirrorTrue : mirrorFalse;
                double initInPathLength = inOffsetPathLengths[my][mx];
                Vector out = outAngles[th][ph];
                double outPathLength = calcOutOffsetPathLength(mx, my, out);
                double additionalPl = initInPathLength + outPathLength;
                double additionalPhase = (additionalPl / lambda) * 2 * Math.PI;
                Complex referenceField = referenceMirror[th][ph];
                double centerDistance = gaussians[my][mx];
                double r = referenceField.abs() * gaussian(centerDistance, beamDiameter);
                double p = referenceField.arg() + additionalPhase;
                Complex field = new Complex(r*Math.cos(p), r*Math.sin(p));
                finalField[th][ph].add(field);
            }
        }
    }
    
    /**
     * calculates on CPU in the interesting area {@link #setCpuInterestingArea(int)}
     * the field of all out angles
     * @see #calcCpuOutAngle(int, int) 
     */
    void computeCpuOutAngles() {
        for (int th = 0; th < tMax; th++) {
            for (int ph = 0; ph < pMax; ph++) {
                finalField[th][ph] = new Complex(0, 0);
            }
        }
        for (int th = trueThStart; th < trueThEnd; th++) {
            for (int ph = truePhStart; ph < truePhEnd; ph++) {
                calcCpuOutAngle(th, ph);
            }
        }
        for (int th = falseThStart; th < falseThEnd; th++) {
            for (int ph = falsePhStart; ph < falsePhEnd; ph++) {
                calcCpuOutAngle(th, ph);
            }
        }
    }
    
    /**
     * calculates on GPU via CUDA the field of all out angles
     */
    void computeGpuOutAngles() {
        
        // manage gpu memory
        JCudaEngine.FloatArray jcMirrorTrue = gpuEngine.createFloatArray(mirrorTrue);
        JCudaEngine.FloatArray jcMirrorFalse = gpuEngine.createFloatArray(mirrorFalse);
        JCudaEngine.FloatArray jcInOffsetPathLengths = gpuEngine.createFloatArray(inOffsetPathLengths);
        jcMirrorTrue.hostToDevice();
        jcMirrorFalse.hostToDevice();
        jcInOffsetPathLengths.hostToDevice();
        JCudaEngine.FloatArray jcFinalField = gpuEngine.createFloatArray(finalField);
        
        // calculations before launching CUDA kernels
        long refTime = 64 * 64 * (long) (900 * 900);
        long nrMirrors = nrX * nrY;
        long nrPixels = pMax * tMax;
        long thisTime = nrMirrors * nrPixels;
        long nrKernels = thisTime / refTime;
        nrKernels = nrKernels < 1 ? 1 : nrKernels;
        int mirrorsPerKernel = (int) (nrMirrors / nrKernels);
        int nrP = (int) nrPixels;
        if (nrP != nrPixels) throw new RuntimeException("to many pixels");
        
        // launching CUDA kernels
        int blockSizeX = 256;
        for (int k = 0; k < nrKernels; k++) {
            System.out.println(k + "/" + nrKernels);
            int mpk = -1;
            if (k == nrKernels - 1) mpk = (int) (mirrorsPerKernel + nrMirrors % nrKernels);
            else mpk = mirrorsPerKernel;
            int mStart = mirrorsPerKernel * k;
            int mEnd = mStart + mpk;
            if (mpk <= 0) throw new RuntimeException("mirrorsPerKernel error");
            gpuEngine.launchKernel(gpuModuleName, gpuFunctionName, blockSizeX, nrP, mStart, mEnd, jcTiltState,
                jcMirrorTrue, jcMirrorFalse, jcInOffsetPathLengths, jcGaussians, jcFinalField);
            
        }
        
        // coping results from gpu to cpu
        jcFinalField.deviceToHost();
        finalField = jcFinalField.toComplex2d(tMax, pMax);
    }
    
    Complex[][] computeGpuDeltaPeaks() {
        int nrPixels = pMax * tMax;
        int blockSizeX = 256;
        gpuEngine.launchKernel("JCudaCalcDeltaPeaks", "calcDeltaPeaks",
                blockSizeX, nrPixels, mirrorSize+gap, inBeam.getX(),
                inBeam.getY(), jcField);
        
        // coping results from gpu to cpu
        jcField.deviceToHost();
        return jcField.toComplex2d(tMax, pMax);
    }
    
    Complex[][] computeGpuSingleMirror(boolean mirrorState) {
        int nrPixels = pMax * tMax;
        int blockSizeX = 256;
        double tiltR = tiltD / 180 * Math.PI;
        double alpha = mirrorState ? tiltR : -tiltR;
        gpuEngine.launchKernel("JCudaCalcSingleMirror", "calcSingleMirror",
                blockSizeX, nrPixels, mirrorSize, inBeam.getX(),
                inBeam.getY(), inBeam.getZ(), alpha, jcField);
        
        // coping results from gpu to cpu
        jcField.deviceToHost();
        return jcField.toComplex2d(tMax, pMax);
    }
    
    Image[] computeGpuBothSingleMirrors() {
        computeAnalyticSingleMirrors();
        //mirrorTrue = computeGpuSingleMirror(true);
        //mirrorFalse = computeGpuSingleMirror(false);
        Image viewTrue = buildIntensityImage(mirrorTrue);
        Image viewFalse = buildIntensityImage(mirrorFalse);
        Image[] trueFalse = {viewTrue, viewFalse};
        return trueFalse;
    }
    
    /**
     * creates an intensity image from a 2D complex field array
     * @param field 2D field array
     * @return intensity image
     */
    static Image buildIntensityImage(Complex[][] field) {
        int width = field[0].length;
        int height = field.length;
        Image view = new Image(width, height);
        //view.setTitle((int) (lambda*1000) + "-intensity");
        for (int th = 0; th < height; th++) {
            for (int ph = 0; ph < width; ph++) {
                view.set(ph, th, (float) field[th][ph].square());
            }
        }
        return view;
    }
    
    /**
     * creates an phase image from a 2D complex field array
     * @param field 2D field array
     * @return phase image
     */
    static Image buildPhaseImage(Complex[][] field) {
        int width = field[0].length;
        int height = field.length;
        Image view = new Image(width, height);
        //view.setTitle((int) (lambda*1000) + "-phase");
        for (int th = 0; th < height; th++) {
            for (int ph = 0; ph < width; ph++) {
                view.set(ph, th, (float) field[th][ph].arg());
            }
        }
        return view;
    }
    
    /**
     * simulates the field for a specified incident beam vector
     * @param inBeam incident beam vector
     * @return image array of
     * [0,1,2,3]=[intensity, true reference field, false reference field, phase]
     */
    Image[] simulateFieldAccurate(Vector inBeam) {
        setInBeam(inBeam);
        Image[] trueFalse = gpuActive ? computeGpuBothSingleMirrors() : computeAnalyticSingleMirrors();
        if (gpuActive) computeGpuOutAngles();
        else computeCpuOutAngles();
        
        Image intesity = buildIntensityImage(finalField);
        Image phase = buildPhaseImage(finalField);
        Image[] retArray = {intesity, trueFalse[0], trueFalse[1], phase};
        return retArray;
    }
    
    Image[] simulateFieldCoarse(Vector inBeam, boolean flipState) {
        setInBeam(inBeam);
        Image trueFalse = buildIntensityImage(computeGpuSingleMirror(flipState));
        //Image cpuTrueFalse = computeAnalyticSingleMirrors()[1];
        Image deltaPeaks = buildIntensityImage(computeGpuDeltaPeaks());
        //Image cpuDeltaPeaks = buildIntensityImage(calcAllDeltaPeaks());
        return new Image[] {trueFalse, deltaPeaks};
    }
    
    /**
     * calculates the value of a gaussian distribution
     * https://en.wikipedia.org/wiki/Normal_distribution
     * @param x argument of the gaussian function
     * @param m is the mean or expectation of the distribution (and also its median and mode)
     * @param fwhm fill width half maximum
     * @return value of the gaussian function at x
     */
    static double gaussian(double x, double m, double fwhm) {
        double sigma = fwhm/2.3548;
        sigma = sigma*1.41421; // to get gaussian intensity instead of field amplitude
        return 1/(sigma*Math.sqrt(2*Math.PI))*Math.exp(-0.5*Math.pow(((x-m)/sigma),2));
    }
    
    /**
     * calculates the value of a gaussian distribution
     * https://en.wikipedia.org/wiki/Normal_distribution
     * @param x argument of the gaussian function
     * @param fwhm fill width half maximum
     * @return value of the gaussian function at x
     */
    static double gaussian(double x, double fwhm) {
        return gaussian(x,0,fwhm);
    }

    /**
     * main method to simulate images for one single color. The resulting images
     * can be used to find the best incident beam vector
     * @param meta all necessary meta data
     * @param tilt simulate for true or false tilt states of all mirrors
     * @param saveImgs saveing all simulated images? needs a lot of ram
     */
    static void simulateColorSlow(MetaData meta, boolean tilt, boolean saveImgs) {
        
        // reading the meta data        
        String outDir = meta.outDir;
        double beamDiameter = meta.beamDiameter;
        int nrX = meta.nrX;
        int nrY = meta.nrY;
        double mirrorSize = Math.sqrt(meta.latticeConstant*meta.latticeConstant*meta.fillFactor);
        double gap = meta.latticeConstant-mirrorSize;
        double tiltAngle = meta.tiltAngle;
        
        double phiMin = meta.phiOutStart;
        double phiMax = meta.phiOutEnd;
        double thetaMin = meta.thetaOutStart;
        double thetaMax = meta.thetaOutEnd;
        double outStepSize = meta.outStepSize;
        
        
        double thetaInStart = meta.thetaInStart;
        double thetaInEnd = meta.thetaInEnd;
        double phiInStart = meta.phiInStart;
        double phiInEnd = meta.phiInEnd;
        double inStepSize = meta.inStepSize;
        
        int lambda = meta.lambdas[0];
        boolean gpuActive = meta.gpuActive;
        //Image bmp = meta.bmp;
        Image bmp = new Image(meta.nrX, meta.nrY);
        //bmp.setAll(tilt ? 1 : 0);
        bmp = meta.bmp;
        double lambdaUm = lambda * 0.001;
        
        // pre calculations
        double zeroPhiOut = -phiMin / outStepSize;
        double zeroThetaOut = -thetaMin / outStepSize;
        int width = (int) ((phiMax - phiMin) / outStepSize);
        int height = (int) ((thetaMax - thetaMin) / outStepSize);
        int thetaInSteps = (int) ((thetaInEnd - thetaInStart) / inStepSize);
        int phiInSteps = (int) ((phiInEnd - phiInStart) / inStepSize);
        int inSteps = thetaInSteps * phiInSteps;
        double maxIntensity = -Double.MAX_VALUE;
        int maxIntTh = -1, maxIntPh = -1;
        double minDistance = Double.MAX_VALUE;
        int minDisTh = -1, minDisPh = -1;
        int inCounter = 0;
        
        // init dmd and ray tracer
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        DmdSimulator ort = new DmdSimulator(dmd, tiltAngle, lambdaUm,
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp, gpuActive);
        
        // init images for simulations
        Image intensity = new Image(width, height);
        Image phase = new Image(width, height);
        Image intensityMax = new Image(phiInSteps, thetaInSteps);
        Image trueFalse = new Image(width, height);
        Image intensityDistances = new Image(phiInSteps, thetaInSteps);
        Image trueFalseDistances = new Image(phiInSteps, thetaInSteps);
        intensity.setTitle(String.valueOf(lambda) + "_intensity");
        phase.setTitle(String.valueOf(lambda) + "_phase");
        intensityMax.setTitle(lambda + "_intensityMax");
        trueFalse.setTitle("_trueFalse");
        intensityDistances.setTitle(lambda + "_" +  tilt + "_intensity_distances");
        trueFalseDistances.setTitle(lambda + "_" +  tilt + "_trueFalse_distances");
        intensity.show();
        phase.show();
        intensityMax.show();
        trueFalse.show();
        intensityDistances.show();
        trueFalseDistances.show();
        ImageStack intensityStack = new ImageStack(width, height);
        ImageStack phaseStack = new ImageStack(width, height);
        ImageStack trueFalseStack = new ImageStack(width, height);
        
        // simulating over all in angles
        for (int th = 0; th < thetaInSteps; th++) {
            double thetaIn = thetaInStart + th * inStepSize;
            for (int ph = 0; ph < phiInSteps; ph++) {
                double phiIn = phiInStart + ph * inStepSize;
                long startTime = System.currentTimeMillis();
                
                // set in beam
                Vector inBeam = new Vector(phiIn / 180. * Math.PI, thetaIn / 180. * Math.PI);
                inBeam.times(-1);
                
                // the actual simulation
                Image[] simulateField = ort.simulateFieldAccurate(inBeam);
                
                // add simulated images to stacks for saving
                if (saveImgs) {
                    intensityStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_intensity", simulateField[0].getFloatProcessor());
                    phaseStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_phase", simulateField[3].getFloatProcessor());
                    trueFalseStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_" + tilt, simulateField[tilt ? 1 : 2].getFloatProcessor());
                }
                
                // set simulated images to corresponding images for screen
                intensity.set(simulateField[0]);
                phase.set(simulateField[3]);
                trueFalse.set(simulateField[tilt ? 1 : 2]);
                
                // calculates values for result images for finding 1 color optimum
                int[] max0 = intensity.findMax();
                int[] max1 = trueFalse.findMax();
                double intensityDistance = Math.sqrt(Math.pow(max1[0]-zeroPhiOut, 2)
                        + Math.pow(max1[1]-zeroThetaOut, 2));
                intensityDistances.set(ph, th, (float) intensityDistance);
                double trueFaleIntensityDistance = Math.sqrt(Math.pow(max1[0]-max0[0], 2)
                        + Math.pow(max1[1]-max0[1], 2));
                if (max0[2] > maxIntensity) {
                    maxIntensity = max0[2];
                    maxIntTh = th;
                    maxIntPh = ph;
                }
                if (trueFaleIntensityDistance < minDistance) {
                    minDistance = trueFaleIntensityDistance;
                    minDisTh = th;
                    minDisPh = ph;
                }
                
                // set calculates values for 1 color optimum into corresponding images
                intensityMax.set(ph, th, max0[2]);
                trueFalseDistances.set(ph, th, (float) trueFaleIntensityDistance);
                
                // updateing images on the screen
                intensity.repaint();
                phase.repaint();
                intensityMax.repaint();
                trueFalse.repaint();
                intensityDistances.repaint();
                trueFalseDistances.repaint();
                
                inCounter++;
                
                // calculation for additional information in console
                double calcTime = ((System.currentTimeMillis() - startTime) / 1000. / 60.);
                int expectedTime = (int) ((inSteps - inCounter) * calcTime);
                double progres = (int) Math.round(
                        (float)((inCounter / (double) inSteps)) * 10000.0) * 0.01;
                double memmory = (int) Math.round((float) (Runtime.getRuntime().totalMemory()
                        / (double) Runtime.getRuntime().maxMemory()) * 10000.0) * 0.01;
                System.out.println("phiIn " + phiIn + " ; thetaIn " + thetaIn
                        + " ; intensityMax " + max0[2] + " ; distance "
                        + trueFaleIntensityDistance + " ; min " + expectedTime
                        + " ; progress " + progres + "%" + " ; memmory " + memmory + "%");
            }
        }
        
        // saveing results on the hard disc
        intensityMax.saveAsTiff(
                outDir + lambda + "_" + tilt + "_intensityMax.tif", meta);
        intensityDistances.saveAsTiff(
                outDir + lambda + "_" + tilt + "_intensity_distances.tif", meta);
        trueFalseDistances.saveAsTiff(
                outDir + lambda + "_" + tilt + "_trueFalse_distances.tif", meta);
        if (saveImgs) {
            ImagePlus intensityImagePlus = new ImagePlus(
                    lambda + "_" + tilt + "_intensity.tif", intensityStack);
            if (intensityImagePlus.getImageStackSize() > 1)
                intensityImagePlus = HyperStackConverter.toHyperStack(
                        intensityImagePlus, 1, thetaInSteps, phiInSteps,
                        "xyztc", "composite");
            intensityImagePlus.setProperty("Info", meta.toString());
            new FileSaver(intensityImagePlus).saveAsTiff(
                    outDir + lambda + "_intensity.tif");
            ImagePlus phaseImagePlus = new ImagePlus(
                    lambda + "_" + tilt + "_phase.tif", phaseStack);
            if (phaseImagePlus.getImageStackSize() > 1) phaseImagePlus =
                    HyperStackConverter.toHyperStack(phaseImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            phaseImagePlus.setProperty("Info", meta.toString());
            new FileSaver(phaseImagePlus).saveAsTiff(
                    outDir + lambda + "_phase.tif");
            ImagePlus trueFalseImagePlus = new ImagePlus(
                    lambda + "_" + tilt + "_" + tilt + ".tif", trueFalseStack);
            if (trueFalseImagePlus.getImageStackSize() > 1) trueFalseImagePlus =
                    HyperStackConverter.toHyperStack(trueFalseImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            trueFalseImagePlus.setProperty("Info", meta.toString());
            new FileSaver(trueFalseImagePlus).saveAsTiff(
                    outDir + lambda + "_" + tilt + ".tif");
        }
        
        // cloeses all images
        intensity.close();
        phase.close();
        trueFalse.close();
        intensityMax.close();
        intensityDistances.close();
        trueFalseDistances.close();
    }
    
    /**
     * simulates the dmd faster but can not simulate patterns on the dmd
     * only all mirror in true or false state
     * @param meta meta data
     * @param lambda wavelength in nm
     * @param tiltState tilt states for all mirrors
     */
    static void simulateColorFast(MetaData meta, int lambda, boolean tiltState) {
        
        // reading the meta data        
        String outDir = meta.outDir;
        double beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);;
        int nrX = meta.nrX;
        int nrY = meta.nrY;
        double mirrorSize = Math.sqrt(meta.latticeConstant*meta.latticeConstant*meta.fillFactor);
        double gap = meta.latticeConstant-mirrorSize;
        double tiltAngle = meta.tiltAngle;
        
        double phiMin = meta.phiOutStart;
        double phiMax = meta.phiOutEnd;
        double thetaMin = meta.thetaOutStart;
        double thetaMax = meta.thetaOutEnd;
        double outStepSize = meta.outStepSize;
        
        
        double thetaInStart = meta.thetaInStart;
        double thetaInEnd = meta.thetaInEnd;
        double phiInStart = meta.phiInStart;
        double phiInEnd = meta.phiInEnd;
        double inStepSize = meta.inStepSize;
        
        boolean gpuActive = meta.gpuActive;
        double lambdaUm = lambda * 0.001;
        
        // pre calculations
        int width = (int) ((phiMax - phiMin) / outStepSize);
        int height = (int) ((thetaMax - thetaMin) / outStepSize);
        int thetaInSteps = (int) ((thetaInEnd - thetaInStart) / inStepSize);
        int phiInSteps = (int) ((phiInEnd - phiInStart) / inStepSize);
        int inSteps = thetaInSteps * phiInSteps;
        int inCounter = 0;
        
        // init dmd and ray tracer
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        DmdSimulator ort = new DmdSimulator(dmd, tiltAngle, lambdaUm,
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, new Image(meta.nrX, meta.nrY), gpuActive);
        
        // init images for simulations
        Image gratingPeaks = new Image(width, height);
        //Image intensityOn = new Image(width, height);
        Image intensityOff = new Image(width, height);
        //Image envelopeOn = new Image(width, height);
        Image envelopeOff = new Image(width, height);
        //Image intensityMaxXOn = new Image(phiInSteps, thetaInSteps);
        Image intensityMaxXOff = new Image(phiInSteps, thetaInSteps);
        //Image intensityMaxYOn = new Image(phiInSteps, thetaInSteps);
        Image intensityMaxYOff = new Image(phiInSteps, thetaInSteps);
        //Image enveopePeakDistanceOn = new Image(phiInSteps, thetaInSteps);
        Image enveopePeakDistanceOff = new Image(phiInSteps, thetaInSteps);
        
        gratingPeaks.setTitle(String.valueOf(lambda) + "_grating-peaks");
        //intensityOn.setTitle(String.valueOf(lambda) + "_on-intensity");
        intensityOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity");
        //envelopeOn.setTitle(String.valueOf(lambda) + "_on-evelope");
        envelopeOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-evelope");
        //intensityMaxXOn.setTitle(String.valueOf(lambda) + "_on-intensity-max-x");
        intensityMaxXOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity-max-x");
        //intensityMaxYOn.setTitle(String.valueOf(lambda) + "_on-intensity-max-y");
        intensityMaxYOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity-max-y");
        //enveopePeakDistanceOn.setTitle(String.valueOf(lambda) + "_on-eveope-peak-distance");
        enveopePeakDistanceOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-eveope-peak-distance");

        gratingPeaks.show();
        //intensityOn.show();
        intensityOff.show();
        //envelopeOn.show();
        envelopeOff.show();
        //intensityMaxXOn.show();
        intensityMaxXOff.show();
        //intensityMaxYOn.show();
        intensityMaxYOff.show();
        //enveopePeakDistanceOn.show();
        enveopePeakDistanceOff.show();
        
        // simulating over all in angles
        for (int th = 0; th < thetaInSteps; th++) {
            double thetaIn = thetaInStart + th * inStepSize;
            for (int ph = 0; ph < phiInSteps; ph++) {
                double phiIn = phiInStart + ph * inStepSize;
                long startTime = System.currentTimeMillis();
                
                // set in beam
                Vector inBeam = new Vector(phiIn / 180. * Math.PI, thetaIn / 180. * Math.PI);
                inBeam.times(-1);
                
                // the actual simulation
                ort.setInBeam(inBeam);
                gratingPeaks.set(buildIntensityImage(ort.computeGpuDeltaPeaks()));
                //envelopeOn.set(buildIntensityImage(ort.computeGpuSingleMirror(true)));
                envelopeOff.set(buildIntensityImage(ort.computeGpuSingleMirror(false)));
                //intensityOn.set(Image.multiply(envelopeOn, gratingPeaks));
                intensityOff.set(Image.multiply(envelopeOff, gratingPeaks));
                
                //int[] envelopeOnMax = envelopeOn.findMax();
                int[] envelopeOffMax = envelopeOff.findMax();
                //int[] intensityOnMax = intensityOn.findMax();
                int[] intensityOffMax = intensityOff.findMax();
                //double intensityEnvelopeDistanceOn = Math.sqrt((Math.pow(envelopeOnMax[0]-intensityOnMax[0], 2) + Math.pow(envelopeOnMax[1]-intensityOnMax[1], 2)))*outStepSize;
                double intensityEnvelopeDistanceOff = Math.sqrt((Math.pow(envelopeOffMax[0]-intensityOffMax[0], 2) + Math.pow(envelopeOffMax[1]-intensityOffMax[1], 2)));
                
                // set calculates values for 1 color optimum into corresponding images
                //intensityMaxXOn.set(ph, th, (float) (intensityOnMax[0]*outStepSize+phiMin));
                //intensityMaxYOn.set(ph, th, (float) (intensityOnMax[1]*outStepSize+thetaMin));
                intensityMaxXOff.set(ph, th, (float) (intensityOffMax[0]*outStepSize+phiMin));
                intensityMaxYOff.set(ph, th, (float) (intensityOffMax[1]*outStepSize+thetaMin));
                //enveopePeakDistanceOn.set(ph, th, (float) intensityEnvelopeDistanceOn);
                enveopePeakDistanceOff.set(ph, th, (float) intensityEnvelopeDistanceOff);
                
                // updateing images on the screen
                gratingPeaks.repaint();
                //envelopeOn.repaint();
                envelopeOff.repaint();
                //intensityOn.repaint();
                intensityOff.repaint();
                //intensityMaxXOn.repaint();
                //intensityMaxYOn.repaint();
                intensityMaxXOff.repaint();
                intensityMaxYOff.repaint();
                //enveopePeakDistanceOn.repaint();
                enveopePeakDistanceOff.repaint();
                
                inCounter++;
                
                // calculation for additional information in console
                if (inCounter % 100 == 0) {
                    double calcTime = ((System.currentTimeMillis() - startTime) / 1000. / 60.);
                    int expectedTime = (int) ((inSteps - inCounter) * calcTime);
                    double progres = (int) Math.round(
                            (float)((inCounter / (double) inSteps)) * 10000.0) * 0.01;
                    double memmory = (int) Math.round((float) (Runtime.getRuntime().totalMemory()
                            / (double) Runtime.getRuntime().maxMemory()) * 10000.0) * 0.01;
                    System.out.println("phiIn " + phiIn + " ; thetaIn " + thetaIn
                            + " ; min " + expectedTime + " ; progress " + progres
                            + "%" + " ; memmory " + memmory + "%" + " ; fps " + 1./calcTime/60.);
                }
            }
        }
        
        // saveing results on the hard disc
        //intensityMaxXOn.saveAsTiff(outDir + intensityMaxXOn.getTitle() + ".tif", meta);
        //intensityMaxYOn.saveAsTiff(outDir + intensityMaxYOn.getTitle() + ".tif", meta);
        intensityMaxXOff.saveAsTiff(outDir + intensityMaxXOff.getTitle() + ".tif", meta);
        intensityMaxYOff.saveAsTiff(outDir + intensityMaxYOff.getTitle() + ".tif", meta);
        //enveopePeakDistanceOn.saveAsTiff(outDir + enveopePeakDistanceOn.getTitle() + ".tif", meta);
        enveopePeakDistanceOff.saveAsTiff(outDir + enveopePeakDistanceOff.getTitle() + ".tif", meta);
        
        // cloeses all images
        gratingPeaks.close();
        //envelopeOn.close();
        envelopeOff.close();
        //intensityOn.close();
        intensityOff.close();
        //intensityMaxXOn.close();
        //intensityMaxYOn.close();
        intensityMaxXOff.close();
        intensityMaxYOff.close();
        //enveopePeakDistanceOn.close();
        enveopePeakDistanceOff.close();
    }
    
    static double[] simulateColorFastDiagonal(MetaData meta, int lambda, boolean tiltState) {
        
        // reading the meta data        
        String outDir = meta.outDir;
        double beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);;
        int nrX = meta.nrX;
        int nrY = meta.nrY;
        double mirrorSize = Math.sqrt(meta.latticeConstant*meta.latticeConstant*meta.fillFactor);
        double gap = meta.latticeConstant-mirrorSize;
        double tiltAngle = meta.tiltAngle;
        
        double phiMin = meta.phiOutStart;
        double phiMax = meta.phiOutEnd;
        double thetaMin = meta.thetaOutStart;
        double thetaMax = meta.thetaOutEnd;
        double outStepSize = meta.outStepSize;
        
        double phiInStart = meta.phiInStart;
        double phiInEnd = meta.phiInEnd;
        double inStepSize = meta.inStepSize;
        
        boolean gpuActive = meta.gpuActive;
        double lambdaUm = lambda * 0.001;
        
        // pre calculations
        int width = (int) ((phiMax - phiMin) / outStepSize);
        int height = (int) ((thetaMax - thetaMin) / outStepSize);
        int phiInSteps = (int) ((phiInEnd - phiInStart) / inStepSize);
        
        // init dmd and ray tracer
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        DmdSimulator ort = new DmdSimulator(dmd, tiltAngle, lambdaUm,
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, new Image(meta.nrX, meta.nrY), gpuActive);
        
        // init images for simulations
        Image gratingPeaks = new Image(width, height);
        //Image intensityOn = new Image(width, height);
        Image intensiry = new Image(width, height);
        //Image envelopeOn = new Image(width, height);
        Image envelope = new Image(width, height);

        
        gratingPeaks.setTitle(String.valueOf(lambda) + "_grating-peaks");
        //intensityOn.setTitle(String.valueOf(lambda) + "_on-intensity");
        intensiry.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity");
        //envelopeOn.setTitle(String.valueOf(lambda) + "_on-evelope");
        envelope.setTitle(String.valueOf(lambda) + "_" + tiltState + "-evelope");

        gratingPeaks.show();
        //intensityOn.show();
        intensiry.show();
        //envelopeOn.show();
        envelope.show();
        
        //double[] intensityMaxX = new double[phiInSteps];
        //double[] intensityMaxY = new double[phiInSteps];
        double[] envelopePeakDistance = new double[phiInSteps];;
        
        // simulating over all in angles
        for (int ph = 0; ph < phiInSteps; ph++) {
            double phiIn = phiInStart + ph * inStepSize;

            // set in beam
            Vector inBeam = new Vector(phiIn / 180. * Math.PI, -phiIn / 180. * Math.PI);
            inBeam.times(-1);

            // the actual simulation
            ort.setInBeam(inBeam);
            gratingPeaks.set(buildIntensityImage(ort.computeGpuDeltaPeaks()));
            //envelopeOn.set(buildIntensityImage(ort.computeGpuSingleMirror(true)));
            envelope.set(buildIntensityImage(ort.computeGpuSingleMirror(false)));
            //intensityOn.set(Image.multiply(envelopeOn, gratingPeaks));
            intensiry.set(Image.multiply(envelope, gratingPeaks));

            //int[] envelopeOnMax = envelopeOn.findMax();
            int[] envelopeOffMax = envelope.findMax();
            //int[] intensityOnMax = intensityOn.findMax();
            int[] intensityOffMax = intensiry.findMax();
            //double intensityEnvelopeDistanceOn = Math.sqrt((Math.pow(envelopeOnMax[0]-intensityOnMax[0], 2) + Math.pow(envelopeOnMax[1]-intensityOnMax[1], 2)))*outStepSize;
            double intensityEnvelopeDistanceOff = Math.sqrt((Math.pow(envelopeOffMax[0]-intensityOffMax[0], 2) + Math.pow(envelopeOffMax[1]-intensityOffMax[1], 2)));

            // set calculates values for 1 color optimum into corresponding images
            //intensityMaxXOn.set(ph, th, (float) (intensityOnMax[0]*outStepSize+phiMin));
            //intensityMaxYOn.set(ph, th, (float) (intensityOnMax[1]*outStepSize+thetaMin));
            //intensityMaxX[ph] = intensityOffMax[0]*outStepSize+phiMin;
            //intensityMaxY[ph] = intensityOffMax[1]*outStepSize+thetaMin;
            envelopePeakDistance[ph] = intensityEnvelopeDistanceOff;

            // updateing images on the screen
            gratingPeaks.repaint();
            //envelopeOn.repaint();
            envelope.repaint();
            //intensityOn.repaint();
            intensiry.repaint();
        }
        
        // cloeses all images
        gratingPeaks.close();
        //envelopeOn.close();
        envelope.close();
        //intensityOn.close();
        intensiry.close();
        
        return envelopePeakDistance;
    }
    
    /**
     * main method to run this ray tracer
     *
     * @param args "[0:directory for saving] [1:true or false for gpuActive]" +
     * "[2:lambda0] [3:lambda1] ..." + "\n or \n" + "[0:directory for saving]
     * [1:true or false for gpuActive]" + "[1:beamDiameter in µm] [2:nrX
     * mirrors] [3:nrYmirrors] " + "[4:single mirrorSize in µm] [5:dmd fill
     * factor] " + "[6:tiltAngle] [7:phiOutStart] [8:phiOutEnd] " +
     * "[9:thetaOutStart] [10:thetaOutEnd] [11:outStepSize] " + "[12:phiInStart]
     * [13:phiInEnd] [14:thetaInStart] " + "[15:thetaInEnd] [16:inStepSize]" +
     * "[17:flipState as bitMap path] [18:lambda0] [19:lambda1] ...\n" + "all
     * angles in degree"
     */
    public static void main(String args[]) {
        
        // creating meta data object
        MetaData meta = new MetaData();
        try {
            if (args.length >= 2 && args.length < 19) {
                meta.outDir = args[0] + "/";
                meta.gpuActive = Boolean.valueOf(args[1]);
                
                meta.lambdas = new int[args.length - 1];
                for (int i = 0; i < args.length-2; i++) meta.lambdas[i] = Integer.parseInt(args[i+2]);

                meta.nrX = 100;
                meta.nrY = 100;

                meta.latticeConstant = 7.56;
                meta.fillFactor = 0.92;
                meta.tiltAngle = 11.7;
                
                meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

                meta.phiOutStart = -2.5;
                meta.phiOutEnd = 22.5;
                meta.thetaOutStart = -22.5;
                meta.thetaOutEnd = 2.5;
                meta.outStepSize = 0.005;

                meta.phiInStart = -17.4;
                meta.phiInEnd = -16.4;
                meta.thetaInStart = 17.4;
                meta.thetaInEnd = 18.4;
                meta.inStepSize = 1;

                meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\DLP6500_1,35_1,75_33_wl532_ang2_pha0.bmp");
                //meta.bmp = new Image(meta.nrX, meta.nrY);
                //meta.bmp.setAll(1);
            } else if (args.length >= 19) {
                meta.outDir = args[0] + "/";
                meta.gpuActive = Boolean.valueOf(args[1]);
                meta.lambdas = new int[args.length - 18];
                for (int i = 0; i < args.length-18; i++) meta.lambdas[i] = Integer.parseInt(args[i+18]);

                meta.nrX = Integer.parseInt(args[2]);
                meta.nrY = Integer.parseInt(args[3]);

                meta.latticeConstant = Double.parseDouble(args[4]);
                meta.fillFactor = Double.parseDouble(args[5]);
                meta.tiltAngle = Double.parseDouble(args[6]);

                meta.phiOutStart = Double.parseDouble(args[7]);
                meta.phiOutEnd = Double.parseDouble(args[8]);
                meta.thetaOutStart = Double.parseDouble(args[9]);
                meta.thetaOutEnd = Double.parseDouble(args[10]);
                meta.outStepSize = Double.parseDouble(args[11]);

                meta.phiInStart = Double.parseDouble(args[12]);
                meta.phiInEnd = Double.parseDouble(args[13]);
                meta.thetaInStart = Double.parseDouble(args[14]);
                meta.thetaInEnd = Double.parseDouble(args[15]);
                meta.inStepSize = Double.parseDouble(args[16]);

                try {
                    int flipState = Integer.parseInt(args[17]);
                    meta.bmp = new Image(meta.nrX, meta.nrY);
                    meta.bmp.setAll(flipState);
                } catch (NumberFormatException ex) {
                    meta.bmp = Image.readBitmap(args[17]);
                }
            } else {
                String error = "Usage: [0:directory for saving] [1:true or false for gpuActive]"
                        + "[2:lambda0]  [3:lambda1] ..."
                        + "\n or \n"
                        + "[0:directory for saving] [1:true or false for gpuActive]"
                        + "[1:beamDiameter in µm] [2:nrX mirrors] [3:nrYmirrors] "
                        + "[4:single mirrorSize in µm] [5:dmd fill factor] "
                        + "[6:tiltAngle] [7:phiOutStart] [8:phiOutEnd] "
                        + "[9:thetaOutStart] [10:thetaOutEnd] [11:outStepSize] "
                        + "[12:phiInStart] [13:phiInEnd] [14:thetaInStart] "
                        + "[15:thetaInEnd] [16:inStepSize]"
                        + "[17:flipState as bitMap path] [18:lambda0]  [19:lambda1] ...\n"
                        + "all angles in degree";
                System.err.println(error);
                throw new RuntimeException(error);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        
        // start ray tracer
        simulateColorSlow(meta, false, true);
    }
}
