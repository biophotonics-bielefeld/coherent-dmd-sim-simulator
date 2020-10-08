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

/**
 *
 * @author Mario
 */
public class CoarseDiagonalSimulator extends AbstractSimulator {
    
    boolean tiltState;
    
    CoarseDiagonalSimulator(MetaData meta, boolean tiltState) {
        super(meta);
        this.tiltState = tiltState;
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
            double phiIn = (phiInStart + ph * inStepSize) / 180. * Math.PI;
            phiIn = Math.atan(1/Math.sqrt(2)*Math.tan(phiIn));
            // set in beam
            Vector inBeam = new Vector(phiIn, -phiIn);
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
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = true;
        
        int lambdaStart = 400;
        int lambdaEnd = 700;
        int lambdaStepSize = 1;
        int nrLambdas = (lambdaEnd - lambdaStart) / lambdaStepSize + 1;
        meta.lambdas = new int[(lambdaEnd - lambdaStart) / lambdaStepSize + 1];
        for (int i = 0; i < nrLambdas; i++) {
            meta.lambdas[i] = lambdaStart + i*lambdaStepSize;
            System.out.println(i + " " + meta.lambdas[i]);
        }

        meta.nrX = 20;
        meta.nrY = 20;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 12.0;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -90;
        meta.phiOutEnd = 90;
        meta.thetaOutStart = -90;
        meta.thetaOutEnd = 90;
        meta.outStepSize = 0.1;

        meta.phiInStart = -90;
        meta.phiInEnd = 90;
        //meta.thetaInStart = -60;
        //meta.thetaInEnd = 60;
        meta.inStepSize = 0.2;

        //meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\SLM_0,40_1,75_33_wl532_ang0_pha0.bmp");
        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        //DmdSimulationCore dsc = new DmdSimulationCore(meta);
        
        
//        CoarseDiagonalSimulator cds = new CoarseDiagonalSimulator(meta, lambdaStart, lambdaEnd, lambdaStepSize, false);
//        cds.simulate();
        
        int diagonalInSteps = (int) ((meta.phiInEnd - meta.phiInStart) / meta.inStepSize);
        
        int counter = 0;
        Image epdByLambda = new Image(diagonalInSteps, meta.lambdas.length);
        epdByLambda.setTitle(lambdaStart + "_" + lambdaEnd + "_epd_" + (int) (meta.tiltAngle*10));
        epdByLambda.show();

        int y = 0;
        for(int wavelength : meta.lambdas) {
            System.out.println(wavelength);
            long startTime = System.currentTimeMillis();
            meta.lambdas[0] = wavelength;
            CoarseDiagonalSimulator cds = new CoarseDiagonalSimulator(meta, false);
            double[] epdDiagonal = cds.simulateDiagonal();
            for(int x = 0; x < epdDiagonal.length; x++) {
                System.out.println(x + " " + (wavelength-lambdaStart) + " " + epdDiagonal[x]);
                epdByLambda.set(x, y, (float) epdDiagonal[x]);
            }
            epdByLambda.repaint();
            long endTime = System.currentTimeMillis();
            int timeExpected = (int) ((endTime - startTime)*0.001/60 * (meta.lambdas.length - counter++));
            y++;
            IJ.log(timeExpected + " minutes left");
        }
        
        
            epdByLambda.saveAsTiff(meta.outDir + epdByLambda.getTitle() + ".tif", meta);
            epdByLambda.close();
    }

    @Override
    public void simulate() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
