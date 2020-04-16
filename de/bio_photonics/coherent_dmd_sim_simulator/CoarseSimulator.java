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
public class CoarseSimulator extends AbstractSimulator {
    
    boolean tiltState;
    
    CoarseSimulator(MetaData meta, boolean tiltState) {
        super(meta);
        this.tiltState = tiltState;
    }
    
    @Override
    void simulate() {
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
        for (int th = 0; th < thetaInSteps; th++) {
            double thetaIn = thetaInStart + th * inStepSize;
            for (int ph = 0; ph < phiInSteps; ph++) {
                double phiIn = phiInStart + ph * inStepSize;
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
                    System.out.println("phiIn " + phiIn + " ; thetaIn " + thetaIn
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
    
}
