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
 *
 * @author Mario
 */
public class CoarseSimulator extends AbstractSimulator {
    
    boolean tiltState, saveRaws;
    
    public CoarseSimulator(MetaData meta, boolean tiltState, boolean saveRaws) {
        super(meta);
        this.tiltState = tiltState;
        this.saveRaws = saveRaws;
    }
    
    @Override
    public void simulate() {
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
        
        ImageStack gratingPeaksStack = new ImageStack(width, height);
        ImageStack intensityStack = new ImageStack(width, height);
        ImageStack envelopeStack = new ImageStack(width, height);

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
        int step = 0;
        int nrSteps = thetaInSteps * phiInSteps;
        for (int th = 0; th < thetaInSteps; th++) {
            double thetaIn = thetaInStart + th * inStepSize;
            for (int ph = 0; ph < phiInSteps; ph++) {
                double phiIn = phiInStart + ph * inStepSize;
                IJ.log("Progress: " + step++ + "/" + nrSteps);
                long startTime = System.currentTimeMillis();
                
                // set in beam
                Vector inBeam = new Vector(phiIn / 180. * Math.PI, thetaIn / 180. * Math.PI);
                inBeam.times(-1);
                
                // the actual simulation
                dsc.setInBeam(inBeam);
                gratingPeaks.set(DmdSimulationCore.buildIntensityImage(dsc.calcAnalyticDeltaPeaks()));
                //envelopeOn.set(buildIntensityImage(ort.computeGpuSingleMirror(true)));
                envelopeOff.set(DmdSimulationCore.buildIntensityImage(dsc.calcAnalyticSingleMirror(tiltState)));
                //intensityOn.set(Image.multiply(envelopeOn, gratingPeaks));
                intensityOff.set(Image.multiply(envelopeOff, gratingPeaks));
                
                if (saveRaws) {
                    gratingPeaksStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_gratingPeaks", gratingPeaks.getFloatProcessor());
                    envelopeStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_envelope", envelopeOff.getFloatProcessor());
                    intensityStack.addSlice(lambda + "_" + phiIn + "_" + thetaIn
                            + "_intensity", intensityOff.getFloatProcessor());
                }
                
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
                    IJ.log("phiIn " + phiIn + " ; thetaIn " + thetaIn
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
    
    public static void main(String[] args) {
        // creating meta data object
        MetaData meta = new MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = true;
        
        int lambdaStart = 532;
        int lambdaEnd = 700;
        int lambdaStepSize = 100;
        int nrLambdas = (lambdaEnd - lambdaStart) / lambdaStepSize + 1;
        meta.lambdas = new int[(lambdaEnd - lambdaStart) / lambdaStepSize + 1];
        for (int i = 0; i < nrLambdas; i++) {
            meta.lambdas[i] = lambdaStart + i*lambdaStepSize;
            System.out.println(i + " " + meta.lambdas[i]);
        }

        meta.nrX = 50;
        meta.nrY = 50;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 12.0;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -80;
        meta.phiOutEnd = 80;
        meta.thetaOutStart = -80;
        meta.thetaOutEnd = 80;
        meta.outStepSize = 0.1;

        meta.phiInStart = -60;
        meta.phiInEnd = 60;
        meta.thetaInStart = -60;
        meta.thetaInEnd = 60;
        meta.inStepSize = 0.2;

        //meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\SLM_0,40_1,75_33_wl532_ang0_pha0.bmp");
        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        //DmdSimulationCore dsc = new DmdSimulationCore(meta);
        
        boolean save = false;
        CoarseSimulator cs = new CoarseSimulator(meta, false, save);
        long timeStart = System.currentTimeMillis();
        cs.simulate();
        System.out.println("Time in seconds: " + ((System.currentTimeMillis() - timeStart) * 0.001));
    }
    
}
