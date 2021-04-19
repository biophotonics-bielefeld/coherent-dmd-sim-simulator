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
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.plugin.HyperStackConverter;

/**
 * class which implements the grating approach
 * @author Mario
 */
public class GratingApproachSimulator extends AbstractSimulator {
    
    boolean tiltState, saveRaws;
    
    public GratingApproachSimulator(MetaData meta, boolean tiltState, boolean saveRaws) {
        super(meta);
        this.tiltState = tiltState;
        this.saveRaws = saveRaws;
    }
    
    @Override
    public void simulate() {
        // init images for simulations
        Image gratingPeaks = new Image(width, height);
        Image intensityOff = new Image(width, height);
        Image envelopeOff = new Image(width, height);
        Image intensityMaxXOff = new Image(phiInSteps, thetaInSteps);
        Image intensityMaxYOff = new Image(phiInSteps, thetaInSteps);
        Image enveopePeakDistanceOff = new Image(phiInSteps, thetaInSteps);
        
        gratingPeaks.setTitle(String.valueOf(lambda) + "_grating-peaks");
        intensityOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity");
        envelopeOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-evelope");
        intensityMaxXOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity-max-x");
        intensityMaxYOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity-max-y");
        enveopePeakDistanceOff.setTitle(String.valueOf(lambda) + "_" + tiltState + "-eveope-peak-distance");
        
        ImageStack gratingPeaksStack = new ImageStack(width, height);
        ImageStack intensityStack = new ImageStack(width, height);
        ImageStack envelopeStack = new ImageStack(width, height);

        gratingPeaks.show();
        intensityOff.show();
        envelopeOff.show();
        intensityMaxXOff.show();
        intensityMaxYOff.show();
        enveopePeakDistanceOff.show();
        
        // simulating over all in angles
        int step = 0;
        int nrSteps = thetaInSteps * phiInSteps;
        for (int th = 0; th < thetaInSteps; th++) {
            double thetaIn = thetaInStart + th * inStepSize;
            for (int ph = 0; ph < phiInSteps; ph++) {
                double phiIn = phiInStart + ph * inStepSize;
                IJ.log("Progress: " + step++ + "/" + nrSteps);
                
                // set in beam
                Vector inBeam = new Vector(phiIn / 180. * Math.PI, thetaIn / 180. * Math.PI);
                inBeam.times(-1);
                
                // the actual simulation
                Image[] diffractionImages = dsc.simulateGratingApproach(inBeam, tiltState);
                envelopeOff.set(diffractionImages[0]);
                gratingPeaks.set(diffractionImages[1]);
                intensityOff.set(Image.multiply(envelopeOff, gratingPeaks));
                
                if (saveRaws) {
                    gratingPeaksStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_gratingPeaks", gratingPeaks.getFloatProcessor());
                    envelopeStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_envelope", envelopeOff.getFloatProcessor());
                    intensityStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_intensity", intensityOff.getFloatProcessor());
                }
                
                int[] envelopeOffMax = envelopeOff.findMax();
                int[] intensityOffMax = intensityOff.findMax();
                double intensityEnvelopeDistanceOff = Math.sqrt((Math.pow(envelopeOffMax[0]-intensityOffMax[0], 2) + Math.pow(envelopeOffMax[1]-intensityOffMax[1], 2)));
                
                // set calculates values for 1 color optimum into corresponding images
                intensityMaxXOff.set(ph, th, (float) (intensityOffMax[0]*outStepSize+phiMin));
                intensityMaxYOff.set(ph, th, (float) (intensityOffMax[1]*outStepSize+thetaMin));
                enveopePeakDistanceOff.set(ph, th, (float) intensityEnvelopeDistanceOff);
                
                // updateing images on the screen
                gratingPeaks.repaint();
                envelopeOff.repaint();
                intensityOff.repaint();
                intensityMaxXOff.repaint();
                intensityMaxYOff.repaint();
                enveopePeakDistanceOff.repaint();
            }
        }
        
        // saveing results on the hard disc
        intensityMaxXOff.saveAsTiff(outDir + intensityMaxXOff.getTitle() + ".tif", meta);
        intensityMaxYOff.saveAsTiff(outDir + intensityMaxYOff.getTitle() + ".tif", meta);
        enveopePeakDistanceOff.saveAsTiff(outDir + enveopePeakDistanceOff.getTitle() + ".tif", meta);
        
        if (saveRaws) {
            ImagePlus gratingPeakImagePlus = new ImagePlus(
                    lambda + "_coarse_gratingPeaks.tif", gratingPeaksStack);
            if (gratingPeakImagePlus.getImageStackSize() > 1) gratingPeakImagePlus =
                    HyperStackConverter.toHyperStack(gratingPeakImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            gratingPeakImagePlus.setProperty("Info", meta.toString());
            new FileSaver(gratingPeakImagePlus).saveAsTiff(
                    outDir + lambda + "_coarse_gratingPeaks.tif");
            ImagePlus envelopeImagePlus = new ImagePlus(
                    lambda + "_coarse_envelope.tif", envelopeStack);
            if (envelopeImagePlus.getImageStackSize() > 1) envelopeImagePlus =
                    HyperStackConverter.toHyperStack(envelopeImagePlus, 1,
                            thetaInSteps, phiInSteps, "xyztc", "composite");
            envelopeImagePlus.setProperty("Info", meta.toString());
            new FileSaver(envelopeImagePlus).saveAsTiff(
                    outDir + lambda + "_coarse_envelope.tif");
            ImagePlus intensityImagePlus = new ImagePlus(
                    lambda + "_coarse_intensity.tif", intensityStack);
            if (intensityImagePlus.getImageStackSize() > 1)
                intensityImagePlus = HyperStackConverter.toHyperStack(
                        intensityImagePlus, 1, thetaInSteps, phiInSteps,
                        "xyztc", "composite");
            intensityImagePlus.setProperty("Info", meta.toString());
            new FileSaver(intensityImagePlus).saveAsTiff(
                    outDir + lambda + "_coarse_intensity.tif");
        }
        
        // cloeses all images
        gratingPeaks.close();
        envelopeOff.close();
        intensityOff.close();
        intensityMaxXOff.close();
        intensityMaxYOff.close();
        enveopePeakDistanceOff.close();
    }
    
    /**
     * starts the grating approach, values in this method need
     * to be adjusted for the desired system conditions
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // creating meta data object
        MetaData meta = new MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = true;
        
        
        meta.lambdas = new double[]{633};

        meta.nrX = 100;
        meta.nrY = 100;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 12.0;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -25;
        meta.phiOutEnd = -11;
        meta.thetaOutStart = 11;
        meta.thetaOutEnd = 25;
        meta.outStepSize = 7.0/2460.0;

        meta.phiInStart = 0;
        meta.phiInEnd = 1;
        meta.thetaInStart = 0;
        meta.thetaInEnd = 1;
        meta.inStepSize = 1;

        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        boolean save = true;
        GratingApproachSimulator cs = new GratingApproachSimulator(meta, false, save);
        long timeStart = System.currentTimeMillis();
        cs.simulate();
        System.out.println("Time in seconds: " + ((System.currentTimeMillis() - timeStart) * 0.001));
    }
    
}
