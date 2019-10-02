/*
 * Copyright (C) 2019 m.lachetta
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

/**
 *
 * @author m.lachetta
 */
public class SingleAngleColorCalculator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // creating meta data object
        DmdSimulator.MetaData meta = new DmdSimulator.MetaData();
        meta.outDir = args[0] + "/";
        meta.gpuActive = Boolean.valueOf(args[1]);
        
        int lambdaStart = Integer.parseInt(args[2]);
        int lambdaEnd = Integer.parseInt(args[3]);
        meta.lambdas = new int[lambdaEnd - lambdaStart + 1];
        for (int i = 0; i <= lambdaEnd-lambdaStart; i++) meta.lambdas[i] = lambdaStart + i;

        meta.nrX = 20;
        meta.nrY = 20;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 11.7;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -80;
        meta.phiOutEnd = 80;
        meta.thetaOutStart = -80;
        meta.thetaOutEnd = 80;
        meta.outStepSize = 0.1;

        meta.phiInStart = -45;
        meta.phiInEnd = 45;
        //meta.thetaInStart = -60;
        //meta.thetaInEnd = 60;
        meta.inStepSize = 0.2;

        //meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\SLM_0,40_1,75_33_wl532_ang0_pha0.bmp");
        meta.bmp = new Image(meta.nrX, meta.nrY);
        //meta.bmp.setAll(1);
        
        for(int lambda : meta.lambdas) System.out.println(lambda);
        int inSteps = (int) ((meta.phiInEnd - meta.phiInStart) / meta.inStepSize);
        Image epdByLambda = new Image(inSteps, meta.lambdas.length);
        epdByLambda.setTitle(lambdaStart + "_" + lambdaEnd + "_epd");
        epdByLambda.show();
        for(int lambda : meta.lambdas) {
            System.out.println(lambda);
            double[] epdDiagonal = DmdSimulator.simulateColorFastDiagonal(meta, lambda, false);
            for(int i = 0; i < epdDiagonal.length; i++) {
                epdByLambda.set(i, lambda-lambdaStart, (float) epdDiagonal[i]);
            }
            epdByLambda.repaint();
        }
        epdByLambda.saveAsTiff(meta.outDir + epdByLambda.getTitle() + ".tif", meta);
        epdByLambda.close();
    }
    
}
