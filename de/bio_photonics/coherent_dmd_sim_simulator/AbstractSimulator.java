/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.coherent_dmd_sim_simulator;

import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore.MetaData;

/**
 *
 * @author Mario
 */
public abstract class AbstractSimulator {
    MetaData meta;
    // meta data
    final String outDir;
    final double beamDiameter;
    final int nrX;
    final int nrY;
    final double mirrorSize;
    final double gap;
    final double tiltAngle;

    final double phiMin;
    final double phiMax;
    final double thetaMin;
    final double thetaMax;
    final double outStepSize;

    final double thetaInStart;
    final double thetaInEnd;
    final double phiInStart;
    final double phiInEnd;
    final double inStepSize;

    final int lambda;
    boolean gpuActive;
    Image bmp;
    double lambdaUm;
    
    // pre calculations
    final double zeroPhiOut;
    final double zeroThetaOut;
    final int width;
    final int height;
    final int thetaInSteps;
    final int phiInSteps;
    final int inSteps;
    
    double maxIntensity = -Double.MAX_VALUE;
    int maxIntTh = -1, maxIntPh = -1;
    double minDistance = Double.MAX_VALUE;
    int minDisTh = -1, minDisPh = -1;
    int inCounter = 0;
    
    DmdSimulationCore dsc;
    
    AbstractSimulator (MetaData meta) {
        this.meta = meta;
        // reading the meta data
        outDir = meta.outDir + "/";
        beamDiameter = meta.beamDiameter;
        nrX = meta.nrX;
        nrY = meta.nrY;
        mirrorSize = Math.sqrt(meta.latticeConstant*meta.latticeConstant*meta.fillFactor);
        gap = meta.latticeConstant-mirrorSize;
        tiltAngle = meta.tiltAngle;
        
        phiMin = meta.phiOutStart;
        phiMax = meta.phiOutEnd;
        thetaMin = meta.thetaOutStart;
        thetaMax = meta.thetaOutEnd;
        outStepSize = meta.outStepSize;
        
        thetaInStart = meta.thetaInStart;
        thetaInEnd = meta.thetaInEnd;
        phiInStart = meta.phiInStart;
        phiInEnd = meta.phiInEnd;
        inStepSize = meta.inStepSize;
        
        lambda = meta.lambdas[0];
        gpuActive = meta.gpuActive;
        bmp = meta.bmp;
        lambdaUm = lambda * 0.001;
        
        // pre calculations
        zeroPhiOut = -phiMin / outStepSize;
        zeroThetaOut = -thetaMin / outStepSize;
        width = (int) ((phiMax - phiMin) / outStepSize);
        height = (int) ((thetaMax - thetaMin) / outStepSize);
        thetaInSteps = (int) ((thetaInEnd - thetaInStart) / inStepSize);
        phiInSteps = (int) ((phiInEnd - phiInStart) / inStepSize);
        inSteps = thetaInSteps * phiInSteps;
        maxIntensity = -Double.MAX_VALUE;
        maxIntTh = -1; maxIntPh = -1;
        minDistance = Double.MAX_VALUE;
        minDisTh = -1; minDisPh = -1;
        inCounter = 0;
        
        // init dmd and simulation core
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        if (gpuActive) dsc = new GpuDmdSimulationCore(dmd, tiltAngle, lambdaUm, 
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp);
        else dsc = new DmdSimulationCore(dmd, tiltAngle, lambdaUm, 
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp);
    }
    
    abstract void simulate();
}
