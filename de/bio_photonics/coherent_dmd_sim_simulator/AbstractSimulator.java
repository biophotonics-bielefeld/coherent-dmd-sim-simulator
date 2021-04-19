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
package de.bio_photonics.coherent_dmd_sim_simulator;

import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore.MetaData;

/**
 * class which provides a general simuilation framework for simulating DMDs
 * with coherent Light
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

    final double lambda;
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
    //int inCounter = 0;
    
    DmdSimulationCore dsc;
    
    /**
     * casic constructor for all simulators stores meta data and generates the
     * dmd model
     * @param meta 
     */
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
        //inCounter = 0;
        
        // init dmd and simulation core
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        if (gpuActive) dsc = new GpuDmdSimulationCore(dmd, tiltAngle, lambdaUm, 
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp);
        else dsc = new DmdSimulationCore(dmd, tiltAngle, lambdaUm, 
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp);
    }
    
    protected final void reloadDmDSimulationCore(double lambdaUm) {
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        if (gpuActive) dsc = new GpuDmdSimulationCore(dmd, tiltAngle, lambdaUm, 
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp);
        else dsc = new DmdSimulationCore(dmd, tiltAngle, lambdaUm, 
                beamDiameter, phiMin, phiMax, thetaMin, thetaMax, outStepSize, bmp);
    }
    
    /**
     * method to simulate with the desired approach
     */
    abstract public void simulate();
}
