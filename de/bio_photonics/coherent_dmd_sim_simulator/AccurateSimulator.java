/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.coherent_dmd_sim_simulator;

import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore.MetaData;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;

/**
 *
 * @author Mario
 */
public class AccurateSimulator extends AbstractSimulator {
    
    boolean tiltState, saveImgs;
    
    AccurateSimulator(MetaData meta, boolean tiltState, boolean saveImgs) {
        super(meta);
        this.tiltState = tiltState;
        this.saveImgs = saveImgs;
    }
    
    @Override
    void simulate() {
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
        intensityDistances.setTitle(lambda + "_" +  tiltState + "_intensity_distances");
        trueFalseDistances.setTitle(lambda + "_" +  tiltState + "_trueFalse_distances");
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
                Image[] simulateField = dsc.simulateFieldAccurate(inBeam);
                
                // add simulated images to stacks for saving
                if (saveImgs) {
                    intensityStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_intensity", simulateField[0].getFloatProcessor());
                    phaseStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_phase", simulateField[3].getFloatProcessor());
                    trueFalseStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_" + tiltState, simulateField[tiltState ? 1 : 2].getFloatProcessor());
                }
                
                // set simulated images to corresponding images for screen
                intensity.set(simulateField[0]);
                phase.set(simulateField[3]);
                trueFalse.set(simulateField[tiltState ? 1 : 2]);
                
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
        intensityMax.saveAsTiff(outDir + lambda + "_" + tiltState + "_intensityMax.tif", meta);
        intensityDistances.saveAsTiff(outDir + lambda + "_" + tiltState + "_intensity_distances.tif", meta);
        trueFalseDistances.saveAsTiff(outDir + lambda + "_" + tiltState + "_trueFalse_distances.tif", meta);
        if (saveImgs) {
            ImagePlus intensityImagePlus = new ImagePlus(
                    lambda + "_" + tiltState + "_intensity.tif", intensityStack);
            if (intensityImagePlus.getImageStackSize() > 1)
                intensityImagePlus = HyperStackConverter.toHyperStack(
                        intensityImagePlus, 1, thetaInSteps, phiInSteps,
                        "xyztc", "composite");
            intensityImagePlus.setProperty("Info", meta.toString());
            new FileSaver(intensityImagePlus).saveAsTiff(
                    outDir + lambda + "_intensity.tif");
            ImagePlus phaseImagePlus = new ImagePlus(
                    lambda + "_" + tiltState + "_phase.tif", phaseStack);
            if (phaseImagePlus.getImageStackSize() > 1) phaseImagePlus =
                    HyperStackConverter.toHyperStack(phaseImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            phaseImagePlus.setProperty("Info", meta.toString());
            new FileSaver(phaseImagePlus).saveAsTiff(
                    outDir + lambda + "_phase.tif");
            ImagePlus trueFalseImagePlus = new ImagePlus(
                    lambda + "_" + tiltState + "_" + tiltState + ".tif", trueFalseStack);
            if (trueFalseImagePlus.getImageStackSize() > 1) trueFalseImagePlus =
                    HyperStackConverter.toHyperStack(trueFalseImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            trueFalseImagePlus.setProperty("Info", meta.toString());
            new FileSaver(trueFalseImagePlus).saveAsTiff(outDir + lambda + "_" + tiltState + ".tif");
        }
        
        // cloeses all images
        intensity.close();
        phase.close();
        trueFalse.close();
        intensityMax.close();
        intensityDistances.close();
        trueFalseDistances.close();
    }
    
}
