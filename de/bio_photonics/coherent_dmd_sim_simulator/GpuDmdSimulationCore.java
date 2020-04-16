/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.coherent_dmd_sim_simulator;

import de.bio_photonics.dmd_ray_tracer.jcuda.JCudaEngine;
import java.io.IOException;

/**
 *
 * @author Mario
 */
public class GpuDmdSimulationCore extends DmdSimulationCore {
    static boolean gpuInit;
    static JCudaEngine gpuEngine;
    static String gpuModuleName;
    static String gpuFunctionName;
    int jcPixelsPerKernel;
    JCudaEngine.IntArray jcTiltState;
    JCudaEngine.FloatArray jcGaussians;
    JCudaEngine.FloatArray jcField;
    
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
    GpuDmdSimulationCore(Dmd dmd, double tiltAngle, double lambda,
            double beamDiameter, double phMin, double phMax, double thMin,
            double thMax, double stepSize, Image flipStates) {
        
        super(dmd, tiltAngle, lambda, beamDiameter, phMin, phMax, thMin, thMax,
                stepSize, flipStates);
        
        init();
    }
    
    GpuDmdSimulationCore(MetaData meta) {
        super(meta);
        init();
    }
    
    private void init() {
        if (gpuInit) return;
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
            throw new RuntimeException("loading GPU modules failed");
        }
        
        gpuEngine.loadFunktion(gpuModuleName, gpuFunctionName);
        gpuEngine.loadFunktion("JCudaCalcDeltaPeaks", "calcDeltaPeaks");
        gpuEngine.loadFunktion("JCudaCalcSingleMirror", "calcSingleMirror");
        gpuInit = true;
        
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
     * calculates in the interesting area {@link #setInterestingArea(int)}
     * the field of all out angles
     * @see #calcCpuOutAngle(int, int) 
     */
    @Override
    protected Complex[][] calcOutAngles() {
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
        return finalField;
    }
    
    /**
     * calculates the complex field of the native grating of the DMD for all out angles
     * @return 
     */
    @Override
    protected Complex[][] calcAnalyticDeltaPeaks() {
        int nrPixels = pMax * tMax;
        int blockSizeX = 256;
        gpuEngine.launchKernel("JCudaCalcDeltaPeaks", "calcDeltaPeaks",
                blockSizeX, nrPixels, mirrorSize+gap, inBeam.getX(),
                inBeam.getY(), jcField);
        
        // coping results from gpu to cpu
        jcField.deviceToHost();
        return jcField.toComplex2d(tMax, pMax);
    }
    
    /**
     * calculates the complex fields for true and false tilt states for all out angles
     * @return intensity distribution for a true and false mirror
     */
    @Override
    protected Complex[][] calcAnalyticSingleMirror(boolean tiltState) {
        int nrPixels = pMax * tMax;
        int blockSizeX = 256;
        double tiltR = tiltD / 180 * Math.PI;
        double alpha = tiltState ? tiltR : -tiltR;
        gpuEngine.launchKernel("JCudaCalcSingleMirror", "calcSingleMirror",
                blockSizeX, nrPixels, mirrorSize, inBeam.getX(),
                inBeam.getY(), inBeam.getZ(), alpha, jcField);
        
        // coping results from gpu to cpu
        jcField.deviceToHost();
        return jcField.toComplex2d(tMax, pMax);
    }
}
