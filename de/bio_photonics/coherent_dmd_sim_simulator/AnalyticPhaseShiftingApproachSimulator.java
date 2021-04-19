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
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;

/**
 * class for implementing the analytic phase shifting approach
 * @author Mario
 */
public class AnalyticPhaseShiftingApproachSimulator extends AbstractSimulator {
    
    boolean saveImgs;
    
    public AnalyticPhaseShiftingApproachSimulator(MetaData meta, boolean saveImgs) {
        super(meta);
        this.saveImgs = saveImgs;
    }
    
    @Override
    public void simulate() {
        // init images for simulations
        Image intensity = new Image(width, height);
        Image phase = new Image(width, height);
        
        Image trueSingleMirrorIntensity = new Image(width, height);
        Image falseSingleMirrorIntensity = new Image(width, height);
        intensity.setTitle(String.valueOf(lambda) + "_intensity");
        phase.setTitle(String.valueOf(lambda) + "_phase");
        trueSingleMirrorIntensity.setTitle(String.valueOf(lambda) + "_trueSingleMirrorIntensity");
        falseSingleMirrorIntensity.setTitle(String.valueOf(lambda) + "_falseSingleMirrorIntensity");
        intensity.show();
        phase.show();
        trueSingleMirrorIntensity.show();
        falseSingleMirrorIntensity.show();
        ImageStack intensityStack = new ImageStack(width, height);
        ImageStack phaseStack = new ImageStack(width, height);
        ImageStack trueStack = new ImageStack(width, height);
        ImageStack falseStack = new ImageStack(width, height);
        
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
                Image[] simulateField = dsc.simulateAnalyticPhaseShiftingApproach(inBeam);
                
                // add simulated images to stacks for saving
                if (saveImgs) {
                    intensityStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_intensity", simulateField[0].getFloatProcessor());
                    phaseStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_phase", simulateField[3].getFloatProcessor());
                    trueStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_true", simulateField[1].getFloatProcessor());
                    falseStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_false", simulateField[2].getFloatProcessor());
                }
                
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
        }
        if (saveImgs) {
            ImagePlus intensityImagePlus = new ImagePlus(
                    lambda + "_intensity.tif", intensityStack);
            if (intensityImagePlus.getImageStackSize() > 1)
                intensityImagePlus = HyperStackConverter.toHyperStack(
                        intensityImagePlus, 1, thetaInSteps, phiInSteps,
                        "xyztc", "composite");
            intensityImagePlus.setProperty("Info", meta.toString());
            new FileSaver(intensityImagePlus).saveAsTiff(
                    outDir + lambda + "_intensity.tif");
            ImagePlus phaseImagePlus = new ImagePlus(
                    lambda + "_phase.tif", phaseStack);
            if (phaseImagePlus.getImageStackSize() > 1) phaseImagePlus =
                    HyperStackConverter.toHyperStack(phaseImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            phaseImagePlus.setProperty("Info", meta.toString());
            new FileSaver(phaseImagePlus).saveAsTiff(
                    outDir + lambda + "_phase.tif");
            ImagePlus trueImagePlus = new ImagePlus(
                    lambda + "_true.tif", trueStack);
            if (trueImagePlus.getImageStackSize() > 1) trueImagePlus =
                    HyperStackConverter.toHyperStack(trueImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            trueImagePlus.setProperty("Info", meta.toString());
            new FileSaver(trueImagePlus).saveAsTiff(outDir + lambda + "_true.tif");
            ImagePlus falseImagePlus = new ImagePlus(
                    lambda + "_false.tif", falseStack);
            if (falseImagePlus.getImageStackSize() > 1) falseImagePlus =
                    HyperStackConverter.toHyperStack(falseImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            falseImagePlus.setProperty("Info", meta.toString());
            new FileSaver(falseImagePlus).saveAsTiff(outDir + lambda + "_false.tif");
        }
        
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
        MetaData meta = new MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = true;
        
        meta.lambdas = new double[]{631};

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

        meta.bmp = Image.readBitmap("D:\\dmd-simulator-images\\cSIM-60x-165-133-473-631\\DLP6500_1,33_1,65_33_wl631_ang0_pha0.bmp");
        //meta.bmp = new Image(meta.nrX, meta.nrY);
        
        
        AnalyticPhaseShiftingApproachSimulator as = new AnalyticPhaseShiftingApproachSimulator(meta, true);
        long timeStart = System.currentTimeMillis();
        as.simulate();
        System.out.println("Time in seconds: " + ((System.currentTimeMillis() - timeStart) * 0.001));
    }
    
}
