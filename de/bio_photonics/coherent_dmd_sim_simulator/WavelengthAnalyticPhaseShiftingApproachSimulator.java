/*
 * Copyright (C) 2021 Mario
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

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;

/**
 *
 * @author Mario
 */
public class WavelengthAnalyticPhaseShiftingApproachSimulator extends AbstractSimulator {
    
    public WavelengthAnalyticPhaseShiftingApproachSimulator(DmdSimulationCore.MetaData meta) {
        super(meta);
    }
    
    @Override
    public void simulate() {
        // init images for simulations
        Image intensity = new Image(width, height);
        Image phase = new Image(width, height);

        Image trueSingleMirrorIntensity = new Image(width, height);
        Image falseSingleMirrorIntensity = new Image(width, height);
        intensity.setTitle(String.valueOf(meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1]) + "_intensity");
        phase.setTitle(String.valueOf(meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1]) + "_phase");
        trueSingleMirrorIntensity.setTitle(String.valueOf(meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1]) + "_trueSingleMirrorIntensity");
        falseSingleMirrorIntensity.setTitle(String.valueOf(meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1]) + "_falseSingleMirrorIntensity");
        intensity.show();
        phase.show();
        trueSingleMirrorIntensity.show();
        falseSingleMirrorIntensity.show();
        ImageStack intensityStack = new ImageStack(width, height);
        ImageStack phaseStack = new ImageStack(width, height);
        ImageStack trueStack = new ImageStack(width, height);
        ImageStack falseStack = new ImageStack(width, height);

        // simulating over all in angles
        for (double lambdaActive : meta.lambdas) {
            reloadDmDSimulationCore(lambdaActive * 0.001);
            System.out.println("Set lambda to: " + lambdaActive);
            double thetaIn = thetaInStart;
            double phiIn = phiInStart;

            // set in beam
            Vector inBeam = new Vector(phiIn / 180. * Math.PI, thetaIn / 180. * Math.PI);
            inBeam.times(-1);

            // the actual simulation
            Image[] simulateField = dsc.simulateAnalyticPhaseShiftingApproach(inBeam);

            // add simulated images to stacks for saving
            intensityStack.addSlice(lambdaActive + "_" + phiIn + "_" + thetaIn
                    + "_intensity", simulateField[0].getFloatProcessor());
            phaseStack.addSlice(lambdaActive + "_" + phiIn + "_" + thetaIn
                    + "_phase", simulateField[3].getFloatProcessor());
            trueStack.addSlice(lambdaActive + "_" + phiIn + "_" + thetaIn
                    + "_true", simulateField[1].getFloatProcessor());
            falseStack.addSlice(lambdaActive + "_" + phiIn + "_" + thetaIn
                    + "_false", simulateField[2].getFloatProcessor());

            // set simulated images to corresponding images for screen
            intensity.set(simulateField[0]);
            phase.set(simulateField[3]);
            trueSingleMirrorIntensity.set(simulateField[1]);
            falseSingleMirrorIntensity.set(simulateField[2]);
//                
            // updateing images on the screen
            intensity.repaint();
            phase.repaint();
            trueSingleMirrorIntensity.repaint();
            falseSingleMirrorIntensity.repaint();
        }
        ImagePlus intensityImagePlus = new ImagePlus(
                meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_intensity.tif", intensityStack);
        if (intensityImagePlus.getImageStackSize() > 1)
            intensityImagePlus = HyperStackConverter.toHyperStack(
                    intensityImagePlus, 1, 1, meta.lambdas.length,
                    "xyztc", "composite");
        intensityImagePlus.setProperty("Info", meta.toString());
        new FileSaver(intensityImagePlus).saveAsTiff(
                outDir + meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_intensity.tif");
        ImagePlus phaseImagePlus = new ImagePlus(
                meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_phase.tif", phaseStack);
        if (phaseImagePlus.getImageStackSize() > 1) phaseImagePlus =
                HyperStackConverter.toHyperStack(phaseImagePlus, 1,
                        1, meta.lambdas.length, "xyztc", "composite");
        phaseImagePlus.setProperty("Info", meta.toString());
        new FileSaver(phaseImagePlus).saveAsTiff(
                outDir + meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_phase.tif");
        ImagePlus trueImagePlus = new ImagePlus(
                meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_true.tif", trueStack);
        if (trueImagePlus.getImageStackSize() > 1) trueImagePlus =
                HyperStackConverter.toHyperStack(trueImagePlus, 1,
                        1, meta.lambdas.length, "xyztc", "composite");
        trueImagePlus.setProperty("Info", meta.toString());
        new FileSaver(trueImagePlus).saveAsTiff(outDir + meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_true.tif");
        ImagePlus falseImagePlus = new ImagePlus(
                meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_false.tif", falseStack);
        if (falseImagePlus.getImageStackSize() > 1) falseImagePlus =
                HyperStackConverter.toHyperStack(falseImagePlus, 1,
                        1, meta.lambdas.length, "xyztc", "composite");
        falseImagePlus.setProperty("Info", meta.toString());
        new FileSaver(falseImagePlus).saveAsTiff(outDir + meta.lambdas[0] + "_" + meta.lambdas[meta.lambdas.length-1] + "_false.tif");

        // cloeses all images
        intensity.close();
        phase.close();
        trueSingleMirrorIntensity.close();
        falseSingleMirrorIntensity.close();
    }
    
    /**
     * starts the analytic phase shifting approach, values in this method need
     * to be adjusted for the desired system conditions
     * @param args 
     */
    public static void main(String[] args) {
        // creating meta data object
        DmdSimulationCore.MetaData meta = new DmdSimulationCore.MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = true;
        
        meta.lambdas = new double[]{628.6215, 628.81568, 629.00986, 629.20404,
            629.39822, 629.5924, 629.78658, 629.98076, 630.17494, 630.36912,
            630.5633, 630.75748, 630.95166, 631.14584, 631.34002, 631.5342,
            631.72838, 631.92256, 632.11674, 632.31092, 632.5051, 632.69928,
            632.89346, 633.08764, 633.28182, 633.476, 633.67018, 633.86436,
            634.05854, 634.25272, 634.4469, 634.64108, 634.83526};

        meta.nrX = 350;
        meta.nrY = 350;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 12.0;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -14.2-0.651361;
        meta.phiOutEnd = -14.2+0.651361+0.00167877;
        meta.thetaOutStart = 14.2-0.651361;
        meta.thetaOutEnd = 14.2+0.651361+0.00167877;
        meta.outStepSize = 0.00167877;

        meta.phiInStart = 34.05;
        meta.phiInEnd = 35.05;
        meta.thetaInStart = -34.05;
        meta.thetaInEnd = -33.05;
        meta.inStepSize = 1.0;

        meta.bmp = Image.readBitmap("D:\\dmd-simulator-images\\cSIM-60x-165-133-473-631\\DLP6500_1,33_1,65_33_wl631_ang2_pha0.bmp");
        //meta.bmp = new Image(meta.nrX, meta.nrY);
        
        
        WavelengthAnalyticPhaseShiftingApproachSimulator was = new WavelengthAnalyticPhaseShiftingApproachSimulator(meta);
        long timeStart = System.currentTimeMillis();
        was.simulate();
        System.out.println("Time in seconds: " + ((System.currentTimeMillis() - timeStart) * 0.001));
    }
}
