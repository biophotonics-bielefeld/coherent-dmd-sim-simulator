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
    
    //int interestingAreaInDegrees;
    boolean saveImgs;
    
    public AccurateSimulator(MetaData meta, boolean saveImgs) {
        super(meta);
//        this.interestingAreaInDegrees = interestingAreaInDegrees;
//        if (this.interestingAreaInDegrees <= 0) 
//            throw new RuntimeException("interestingAreaInDegrees has to be lager than 0"); 
        this.saveImgs = saveImgs;
    }
    
    @Override
    public void simulate() {
        // init images for simulations
        Image intensity = new Image(width, height);
        Image phase = new Image(width, height);
        
        //Image intensityMax = new Image(phiInSteps, thetaInSteps);
        Image trueSingleMirrorIntensity = new Image(width, height);
        Image falseSingleMirrorIntensity = new Image(width, height);
        //Image intensityDistances = new Image(phiInSteps, thetaInSteps);
        //Image trueFalseDistances = new Image(phiInSteps, thetaInSteps);
        intensity.setTitle(String.valueOf(lambda) + "_intensity");
        phase.setTitle(String.valueOf(lambda) + "_phase");
        //intensityMax.setTitle(lambda + "_intensityMax");
        trueSingleMirrorIntensity.setTitle(String.valueOf(lambda) + "_trueSingleMirrorIntensity");
        falseSingleMirrorIntensity.setTitle(String.valueOf(lambda) + "_falseSingleMirrorIntensity");
        //intensityDistances.setTitle(lambda + "_" +  tiltState + "_intensity_distances");
        //trueFalseDistances.setTitle(lambda + "_" +  tiltState + "_trueFalse_distances");
        intensity.show();
        phase.show();
        //intensityMax.show();
        trueSingleMirrorIntensity.show();
        falseSingleMirrorIntensity.show();
        //intensityDistances.show();
        //trueFalseDistances.show();
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
                Image[] simulateField = dsc.simulateFieldAccurate(inBeam);
                
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
                
//                // calculates values for result images for finding 1 color optimum
//                int[] max0 = intensity.findMax();
//                int[] max1 = trueFalse.findMax();
//                double intensityDistance = Math.sqrt(Math.pow(max1[0]-zeroPhiOut, 2)
//                        + Math.pow(max1[1]-zeroThetaOut, 2));
//                intensityDistances.set(ph, th, (float) intensityDistance);
//                double trueFaleIntensityDistance = Math.sqrt(Math.pow(max1[0]-max0[0], 2)
//                        + Math.pow(max1[1]-max0[1], 2));
//                if (max0[2] > maxIntensity) {
//                    maxIntensity = max0[2];
//                    maxIntTh = th;
//                    maxIntPh = ph;
//                }
//                if (trueFaleIntensityDistance < minDistance) {
//                    minDistance = trueFaleIntensityDistance;
//                    minDisTh = th;
//                    minDisPh = ph;
//                }
//                
//                // set calculates values for 1 color optimum into corresponding images
//                intensityMax.set(ph, th, max0[2]);
//                trueFalseDistances.set(ph, th, (float) trueFaleIntensityDistance);
//                
                // updateing images on the screen
                intensity.repaint();
                phase.repaint();
                //intensityMax.repaint();
                trueSingleMirrorIntensity.repaint();
                falseSingleMirrorIntensity.repaint();
                //intensityDistances.repaint();
                //trueFalseDistances.repaint();
                
                inCounter++;
                
                // calculation for additional information in console
                double calcTime = ((System.currentTimeMillis() - startTime) / 1000. / 60.);
                int expectedTime = (int) ((inSteps - inCounter) * calcTime);
                double progres = (int) Math.round(
                        (float)((inCounter / (double) inSteps)) * 10000.0) * 0.01;
                double memmory = (int) Math.round((float) (Runtime.getRuntime().totalMemory()
                        / (double) Runtime.getRuntime().maxMemory()) * 10000.0) * 0.01;
                System.out.println("phiIn " + phiIn + " ; thetaIn " + thetaIn
                        /* + " ; intensityMax " + max0[2] + " ; distance "
                        + trueFaleIntensityDistance */ + " ; min " + expectedTime
                        + " ; progress " + progres + "%" + " ; memmory " + memmory + "%");
            }
        }
        
//        // saveing results on the hard disc
//        intensityMax.saveAsTiff(outDir + lambda + "_" + tiltState + "_intensityMax.tif", meta);
//        intensityDistances.saveAsTiff(outDir + lambda + "_" + tiltState + "_intensity_distances.tif", meta);
//        trueFalseDistances.saveAsTiff(outDir + lambda + "_" + tiltState + "_trueFalse_distances.tif", meta);
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
//        intensityMax.close();
//        intensityDistances.close();
//        trueFalseDistances.close();
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

        meta.phiOutStart = -15;
        meta.phiOutEnd = 15;
        meta.thetaOutStart = -15;
        meta.thetaOutEnd = 15;
        meta.outStepSize = 0.02;

        meta.phiInStart = -21;
        meta.phiInEnd = -20;
        meta.thetaInStart = 21;
        meta.thetaInEnd = 22;
        meta.inStepSize = 1.0;

        //meta.bmp = Image.readBitmap("D:\\dmd-simulator-images\\interesting patterns\\circles-50.bmp");
        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        //DmdSimulationCore dsc = new DmdSimulationCore(meta);
        
        
        AccurateSimulator as = new AccurateSimulator(meta, true);
        long timeStart = System.currentTimeMillis();
        as.simulate();
        System.out.println("Time in seconds: " + ((System.currentTimeMillis() - timeStart) * 0.001));
    }
    
}
