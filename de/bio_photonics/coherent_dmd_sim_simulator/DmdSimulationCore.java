/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.coherent_dmd_sim_simulator;

import ij.IJ;
import java.util.stream.IntStream;

/**
 *
 * @author Mario
 */
public class DmdSimulationCore {
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
    
    /**
     * class to handle meta data
     */
    public static class MetaData {
        public String outDir;
        public int[] lambdas;
        public boolean gpuActive;
        public int beamDiameter;
        
        public int nrX;
        public int nrY;
        
        public double latticeConstant;
        public double fillFactor;
        
        public double tiltAngle;
        
        public double phiOutStart;
        public double phiOutEnd;
        public double thetaOutStart;
        public double thetaOutEnd;
        public double outStepSize;
        
        public double phiInStart;
        public double phiInEnd;
        public double thetaInStart;
        public double thetaInEnd;
        public double inStepSize;
        
        public Image bmp;
        
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
     */
    DmdSimulationCore(Dmd dmd, double tiltAngle, double lambda,
            double beamDiameter, double phMin, double phMax, double thMin,
            double thMax, double stepSize, Image flipStates) {
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
        
        init(flipStates);
    }
    
    DmdSimulationCore(MetaData meta) {
        this.nrX = meta.nrX;
        this.nrY = meta.nrY;
        this.mirrorSize = Math.sqrt(meta.latticeConstant*meta.latticeConstant*meta.fillFactor);
        this.gap = meta.latticeConstant - mirrorSize;
        
        this.tiltD = meta.tiltAngle;
        
        this.dmd = new Dmd(nrX, nrY, meta.latticeConstant, gap);
        
        this.lambda = meta.lambdas[0];
        this.beamDiameter = meta.beamDiameter;
        this.phiMinD = meta.phiOutStart;
        this.phiMaxD = meta.phiOutEnd;
        this.thetaMinD = meta.thetaOutStart;
        this.thetaMaxD = meta.thetaOutEnd;
        this.outStepSizeD = meta.outStepSize;
        
        init(meta.bmp);
    }
    
    private void init(Image flipStates) {
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
        referencePosition = dmd.getCoordinates(0, 0, tiltD, 0, 0);
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
                double x = dmd.getX(mx, my, tiltD, 0, 0);
                double y = dmd.getY(mx, my, tiltD, 0, 0);
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
     * sets a region of interest around the maximum in the reference field for calculations
     * @param rangeOfInterest in degree
     */
    /*
    private void setInterestingArea(int rangeOfInterest) {
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
    */
    
    private int checkThRange(int thValue) {
        if (thValue < 0) return 0;
        else if (thValue > tMax) return tMax;
        else return thValue;
    }
    
    private int checkPhRange(int phValue) {
        if (phValue < 0) return 0;
        else if (phValue > pMax) return pMax;
        else return phValue;
    }
    
    /**
     * calculates the complex field of the native grating of the DMD for all out angles
     * @return 
     */
    protected Complex[][] calcAnalyticDeltaPeaks() {
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
    protected Complex[][] calcAnalyticSingleMirror(boolean tiltState) {
        Complex[][] field = new Complex[tMax][pMax];
        double alpha = tiltD / 180. * Math.PI;
        alpha = tiltState ? alpha : -alpha;
        
        for (int th = 0; th < tMax; th++) {
            for (int ph = 0; ph < pMax; ph++) {
                field[th][ph] = calcAnalyticSingleMirrorOutAngle(outAngles[th][ph], alpha);
            }
        }
        return field;
    }
    
    /**
     * calculates the electric complex field for a specified out angle in an
     * analytical manner (look into dmd.nb for details)
     * @param out vector of the out angle
     * @param alpha tilt angle of the mirror in rad
     * @return the complex field of the out angle
     */
    private Complex calcAnalyticSingleMirrorOutAngle(Vector out, double alpha) {
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
    
    /**
     * calculates in the interesting area {@link #setInterestingArea(int)}
     * the field of all out angles
     * @return 
     * @see #calcCpuOutAngle(int, int) 
     */
    protected Complex[][] calcOutAngles() {
        //System.out.println(trueThStart + " " + trueThEnd + " " + truePhStart + " " + truePhEnd);
        Complex[][] field = new Complex[tMax][pMax];
        for (int th = 0; th < tMax; th++) {
            final int finalTh = th;
            IntStream pStream = IntStream.range(0, pMax).parallel();
            pStream.forEach(ph -> field[finalTh][ph] = calcOutAngle(ph, finalTh));
//            for (int ph = 0; ph < pMax; ph++) {
//                //field[th][ph] = new Complex(0, 0);
//                field[th][ph] = calcOutAngle(ph, th);
//            }
            IJ.log("Progress: " + (th + 1) + "/" + tMax);
        }
        /*
        for (int th = trueThStart; th < trueThEnd; th++) {
            for (int ph = truePhStart; ph < truePhEnd; ph++) {
                field[th][ph].add(calcOutAngle(ph, th));
            }
        }
        for (int th = falseThStart; th < falseThEnd; th++) {
            for (int ph = falsePhStart; ph < falsePhEnd; ph++) {
                field[th][ph].add(calcOutAngle(ph, th));
            }
        }
        */
        return field;
    }
    
    /**
     * calculates the field of a specific out angle and saves it into {@link #finalField}
     * @param ph step index of phi
     * @param th step index of theta
     */
    private Complex calcOutAngle(int ph, int th) {
        Complex field = new Complex(0, 0);
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
                double r = referenceField.abs() * gaussians[my][mx];
                double p = referenceField.arg() + additionalPhase;
                field.add(new Complex(r*Math.cos(p), r*Math.sin(p)));
            }
        }
        return field;
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
     * calculates the offset b
     * @param mx
     * @param my
     * @param flipState
     * @return 
     */
    private double calcInOffsetPathLength(int mx, int my, double flipState) {
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
    private double calcOutOffsetPathLength(int mx, int my, Vector out) {
        //Vector currentPosition = dmd.getCoordinates(mx, my, tiltAngles[my][mx], 0, 0);
        double referencePl = -referencePosition.times(out);
        double currentPl = -dmdPositions[my][mx].times(out);
        return currentPl - referencePl;
    }
    
    /**
     * creates an intensity image from a 2D complex field array
     * @param field 2D field array
     * @return intensity image
     */
    public static Image buildIntensityImage(Complex[][] field) {
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
    public static Image buildPhaseImage(Complex[][] field) {
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
    public Image[] simulateFieldAccurate(Vector inBeam) {
        setInBeam(inBeam);
        mirrorTrue = calcAnalyticSingleMirror(true);
        mirrorFalse = calcAnalyticSingleMirror(false);
        //setInterestingArea(interestingAreaInDegrees);
        Complex[][] field = calcOutAngles();
        
        Image fieldIntensity = buildIntensityImage(field);
        Image trueIntensity = buildIntensityImage(mirrorTrue);
        Image falseIntensity = buildIntensityImage(mirrorFalse);
        Image fieldPhase = buildPhaseImage(field);
        Image[] retArray = {fieldIntensity, trueIntensity, falseIntensity, fieldPhase};
        return retArray;
    }
    
    public Image[] simulateFieldCoarse(Vector inBeam, boolean flipState) {
        setInBeam(inBeam);
        Complex[][] trueFalse = calcAnalyticSingleMirror(flipState);
        Complex[][] deltaPeaks = calcAnalyticDeltaPeaks();
        Image trueFalseIntensity = buildIntensityImage(trueFalse);
        Image deltaPeaksIntensity = buildIntensityImage(deltaPeaks);
        return new Image[] {trueFalseIntensity, deltaPeaksIntensity};
    }
    
    
}
