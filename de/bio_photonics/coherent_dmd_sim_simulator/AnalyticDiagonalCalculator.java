/*
 * Copyright (C) 2020 m.lachetta
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
public class AnalyticDiagonalCalculator {
    
    private static double calcDiffractionOrder(double inAngle, double tiltAngle, double waveLength, double latticeConstant) {
        double inAng = inAngle * Math.PI / 180;
        double tiltAng = tiltAngle * Math.PI / 180;
        double lambda = waveLength * 1e-9;
        double latConst = latticeConstant / Math.sqrt(2) * 1e-6;
        return (Math.sin(inAng) + Math.sin(-inAng - 2*tiltAng)) * latConst / lambda;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        DmdSimulator.MetaData meta = new DmdSimulator.MetaData();
        meta.outDir = args[0] + "/";
        meta.gpuActive = Boolean.valueOf(args[1]);
        
        int lambdaStart = Integer.parseInt(args[2]);
        int lambdaEnd = Integer.parseInt(args[3]);
        meta.lambdas = new int[lambdaEnd - lambdaStart + 1];
        for (int i = 0; i <= lambdaEnd-lambdaStart; i++) meta.lambdas[i] = lambdaStart + i;

        //meta.nrX = 20;
        //meta.nrY = 20;

        meta.latticeConstant = 7.56;
        //meta.fillFactor = 0.92;
        meta.tiltAngle = -12.0;

        //meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        //meta.phiOutStart = -80;
        //meta.phiOutEnd = 80;
        //meta.thetaOutStart = -80;
        //meta.thetaOutEnd = 80;
        //meta.outStepSize = 0.05;

        meta.phiInStart = -90;
        meta.phiInEnd = 90;
        //meta.thetaInStart = -60;
        //meta.thetaInEnd = 60;
        meta.inStepSize = 0.05;

        //meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\SLM_0,40_1,75_33_wl532_ang0_pha0.bmp");
        //meta.bmp = new Image(meta.nrX, meta.nrY);
        int width = (int) ((meta.phiInEnd - meta.phiInStart) / meta.inStepSize);
        int height = meta.lambdas.length;
        Image diagonalEpd = new Image(width, height);
        Image outAngle = new Image(width, height);
        diagonalEpd.setTitle(lambdaStart + "_" + lambdaEnd + "_epd_" + (int)(meta.tiltAngle*10) + "_analytic");
        outAngle.setTitle(lambdaStart + "_" + lambdaEnd + "_out_" + (int)(meta.tiltAngle*10) + "_analytic");
        for (int y = 0; y < height; y++) {
            int waveLength = meta.lambdas[y];
            for (int x = 0; x < width; x++) {
                double phi = meta.phiInStart + meta.inStepSize * x;
                double inAngle = phi;//Math.atan(Math.sqrt(2)*Math.tan(phi*Math.PI/180))*180/Math.PI;
                //System.out.println(inAngle + " " + meta.tiltAngle + " " + waveLength + " " + meta.latticeConstant);
                double n = calcDiffractionOrder(inAngle, meta.tiltAngle, waveLength, meta.latticeConstant);
                float value = (float) Math.pow(Math.sin(n*Math.PI), 2.0);
                //System.out.println(x + " " + y);
                diagonalEpd.set(x, y, value);
                
                outAngle.set(x, y, (float) (-inAngle - 2*meta.tiltAngle));
            }
        }
        diagonalEpd.saveAsTiff(meta.outDir + diagonalEpd.getTitle() + ".tif", meta);
        outAngle.saveAsTiff(meta.outDir + outAngle.getTitle() + ".tif", meta);
        
    }
    
}
