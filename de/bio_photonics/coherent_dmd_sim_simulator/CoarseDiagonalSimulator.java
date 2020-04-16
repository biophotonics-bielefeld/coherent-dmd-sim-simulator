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
public class CoarseDiagonalSimulator extends AbstractSimulator {
    
    int lambdaStart, lambdaEnd, lambdaStepSize;
    boolean tiltState;
    
    CoarseDiagonalSimulator(MetaData meta, int lambdaStart, int lambdaEnd, int lambdaStepSize, boolean tiltState) {
        super(meta);
        this.lambdaStart = lambdaStart;
        this.lambdaEnd = lambdaEnd;
        this.lambdaStepSize = lambdaStepSize;
        this.tiltState = tiltState;
    }
    
    @Override
    void simulate() {
        int diagonalInSteps = (int) ((meta.phiInEnd - meta.phiInStart) / meta.inStepSize);
        
        int counter = 0;
        Image epdByLambda = new Image(diagonalInSteps, meta.lambdas.length);
        epdByLambda.setTitle(lambdaStart + "_" + lambdaEnd + "_epd_" + (int) (meta.tiltAngle*10));
        epdByLambda.show();

        int y = 0;
        for(int lambda : meta.lambdas) {
            System.out.println(lambda);
            long startTime = System.currentTimeMillis();
            double[] epdDiagonal = simulateDiagonal();
            for(int x = 0; x < epdDiagonal.length; x++) {
                System.out.println(x + " " + (lambda-lambdaStart) + " " + epdDiagonal[x]);
                epdByLambda.set(x, y, (float) epdDiagonal[x]);
            }
            epdByLambda.repaint();
            long endTime = System.currentTimeMillis();
            int timeExpected = (int) ((endTime - startTime)*0.001/60 * (meta.lambdas.length - counter++));
            y++;
            System.out.println(timeExpected + " minutes left");
        }
        
        
            epdByLambda.saveAsTiff(meta.outDir + epdByLambda.getTitle() + ".tif", meta);
            epdByLambda.close();
    }
    
    double[] simulateDiagonal() {
        // init images for simulations
        Image gratingPeaks = new Image(width, height);
        //Image intensityOn = new Image(width, height);
        Image intensity = new Image(width, height);
        //Image envelopeOn = new Image(width, height);
        Image envelope = new Image(width, height);
        //ImageStack intensityStack = new ImageStack(width, height);
        
        gratingPeaks.setTitle(String.valueOf(lambda) + "_grating-peaks");
        //intensityOn.setTitle(String.valueOf(lambda) + "_on-intensity");
        intensity.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity");
        //envelopeOn.setTitle(String.valueOf(lambda) + "_on-evelope");
        envelope.setTitle(String.valueOf(lambda) + "_" + tiltState + "-evelope");

        gratingPeaks.show();
        //intensityOn.show();
        intensity.show();
        //envelopeOn.show();
        envelope.show();
        
        //double[] intensityMaxX = new double[phiInSteps];
        //double[] intensityMaxY = new double[phiInSteps];
        double[] envelopePeakDistance = new double[phiInSteps];;
        
        // simulating over all in angles
        for (int ph = 0; ph < phiInSteps; ph++) {
            double phiIn = phiInStart + ph * inStepSize;

            // set in beam
            Vector inBeam = new Vector(phiIn / 180. * Math.PI, -phiIn / 180. * Math.PI);
            inBeam.times(-1);

            // the actual simulation
            dsc.setInBeam(inBeam);
            gratingPeaks.set(DmdSimulationCore.buildIntensityImage(dsc.calcAnalyticDeltaPeaks()));
            //envelopeOn.set(buildIntensityImage(ort.computeGpuSingleMirror(true)));
            envelope.set(DmdSimulationCore.buildIntensityImage(dsc.calcAnalyticSingleMirror(tiltState)));
            //intensityOn.set(Image.multiply(envelopeOn, gratingPeaks));
            intensity.set(Image.multiply(envelope, gratingPeaks));
            //intensityStack.addSlice(intensity.getFloatProcessor());

            //int[] envelopeOnMax = envelopeOn.findMax();
            int[] envelopeOffMax = envelope.findMax();
            //int[] intensityOnMax = intensityOn.findMax();
            int[] intensityOffMax = intensity.findMax();
            //double intensityEnvelopeDistanceOn = Math.sqrt((Math.pow(envelopeOnMax[0]-intensityOnMax[0], 2) + Math.pow(envelopeOnMax[1]-intensityOnMax[1], 2)))*outStepSize;
            double intensityEnvelopeDistanceOff = Math.sqrt((Math.pow(envelopeOffMax[0]-intensityOffMax[0], 2) + Math.pow(envelopeOffMax[1]-intensityOffMax[1], 2)));

            // set calculates values for 1 color optimum into corresponding images
            //intensityMaxXOn.set(ph, th, (float) (intensityOnMax[0]*outStepSize+phiMin));
            //intensityMaxYOn.set(ph, th, (float) (intensityOnMax[1]*outStepSize+thetaMin));
            //intensityMaxX[ph] = intensityOffMax[0]*outStepSize+phiMin;
            //intensityMaxY[ph] = intensityOffMax[1]*outStepSize+thetaMin;
            envelopePeakDistance[ph] = intensityEnvelopeDistanceOff;

            // updateing images on the screen
            gratingPeaks.repaint();
            //envelopeOn.repaint();
            envelope.repaint();
            //intensityOn.repaint();
            intensity.repaint();
        }
        
        // cloeses all images
        gratingPeaks.close();
        //envelopeOn.close();
        envelope.close();
        //intensityOn.close();
        intensity.close();
        //ImagePlus intensities = new ImagePlus("Intensities", intensityStack);
        //new FileSaver(intensities).saveAsTiff(outDir + "/diagonals/" + lambda + ".tif");
        return envelopePeakDistance;
    }
    
    public static void main(String[] args) {
        // creating meta data object
        MetaData meta = new MetaData();
        meta.outDir = "D:\\dmd-simulator-images";
        meta.gpuActive = false;
        
        int lambdaStart = 400;
        int lambdaEnd = 700;
        int lambdaStepSize = 100;
        int nrLambdas = (lambdaEnd - lambdaStart) / lambdaStepSize + 1;
        meta.lambdas = new int[(lambdaEnd - lambdaStart) / lambdaStepSize + 1];
        for (int i = 0; i < nrLambdas; i++) {
            meta.lambdas[i] = lambdaStart + i*lambdaStepSize;
            System.out.println(i + " " + meta.lambdas[i]);
        }

        meta.nrX = 5;
        meta.nrY = 5;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 12.0;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -60;
        meta.phiOutEnd = 60;
        meta.thetaOutStart = -60;
        meta.thetaOutEnd = 60;
        meta.outStepSize = 0.2;

        meta.phiInStart = -45;
        meta.phiInEnd = +45;
        //meta.thetaInStart = -60;
        //meta.thetaInEnd = 60;
        meta.inStepSize = 5.0;

        //meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\SLM_0,40_1,75_33_wl532_ang0_pha0.bmp");
        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        //DmdSimulationCore dsc = new DmdSimulationCore(meta);
        
        
        CoarseDiagonalSimulator cds = new CoarseDiagonalSimulator(meta, lambdaStart, lambdaEnd, lambdaStepSize, false);
        cds.simulate();
    }
}
