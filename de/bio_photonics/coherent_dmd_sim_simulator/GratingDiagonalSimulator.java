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
 * class which implements the grating approach along the diagonal of the dmd
 * @author Mario
 */
public class GratingDiagonalSimulator extends AbstractSimulator {
    
    boolean tiltState;
    
    GratingDiagonalSimulator(MetaData meta, boolean tiltState) {
        super(meta);
        this.tiltState = tiltState;
    }
    
    /**
     * simulates diffraction patterns along the diagonal of the dmd and
     * measures the distance between the center of the envelope and the
     * brightes native grating peak in pixels
     * @return envelope to peak distances in pixels for the desired input angles
     */
    double[] simulateDiagonal() {
        // init images for simulations
        Image gratingPeaks = new Image(width, height);
        Image intensity = new Image(width, height);
        Image envelope = new Image(width, height);
        
        gratingPeaks.setTitle(String.valueOf(lambda) + "_grating-peaks");
        intensity.setTitle(String.valueOf(lambda) + "_" + tiltState + "-intensity");
        envelope.setTitle(String.valueOf(lambda) + "_" + tiltState + "-evelope");

        gratingPeaks.show();
        intensity.show();
        envelope.show();
        
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
            envelope.set(DmdSimulationCore.buildIntensityImage(dsc.calcAnalyticSingleMirror(tiltState)));
            intensity.set(Image.multiply(envelope, gratingPeaks));

            int[] envelopeOffMax = envelope.findMax();
            int[] intensityOffMax = intensity.findMax();
            double intensityEnvelopeDistanceOff = Math.sqrt((Math.pow(envelopeOffMax[0]-intensityOffMax[0], 2) + Math.pow(envelopeOffMax[1]-intensityOffMax[1], 2)));

            envelopePeakDistance[ph] = intensityEnvelopeDistanceOff;

            // updateing images on the screen
            gratingPeaks.repaint();
            envelope.repaint();
            intensity.repaint();
        }
        
        // cloeses all images
        gratingPeaks.close();
        envelope.close();
        intensity.close();
        return envelopePeakDistance;
    }
    
    /**
     * starts the grating approach along the diagonal, values in this method need
     * to be adjusted for the desired system conditions
     * @param args the command line arguments
     */
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
        meta.inStepSize = 0.2;

        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        int diagonalInSteps = (int) ((meta.phiInEnd - meta.phiInStart) / meta.inStepSize);
        
        int counter = 0;
        Image epdByLambda = new Image(diagonalInSteps, meta.lambdas.length);
        epdByLambda.setTitle(lambdaStart + "_" + lambdaEnd + "_epd_" + (int) (meta.tiltAngle*10));
        epdByLambda.show();

        int y = 0;
        for(int wavelength : meta.lambdas) {
            long startTime = System.currentTimeMillis();
            meta.lambdas[0] = wavelength;
            GratingDiagonalSimulator cds = new GratingDiagonalSimulator(meta, false);
            double[] epdDiagonal = cds.simulateDiagonal();
            for(int x = 0; x < epdDiagonal.length; x++) {
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

    /**
     * not implemented in this simulator
     */
    @Override
    public void simulate() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
