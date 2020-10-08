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

/**
 * class which implements the blaze condition approach along the diagonal
 * of the dmd
 * @author Mario
 */
public class AnalyticDiagonalCalculator {
    
    /**
     * calculates the corresponing diffraction order for the input parameters
     * @param inAngle in degrees
     * @param tiltAngle in degrees
     * @param waveLength in nm
     * @param latticeConstant in µm
     * @return corresponding diffraction order
     */
    public static double calcDiffractionOrder(double inAngle, double tiltAngle, double waveLength, double latticeConstant) {
        double inAng = inAngle * Math.PI / 180;
        double tiltAng = tiltAngle * Math.PI / 180;
        double outAng = -inAng + 2*tiltAng;
        double lambda = waveLength * 1e-9;
        double latConst = latticeConstant / Math.sqrt(2) * 1e-6;
        return (Math.sin(inAng) + Math.sin(outAng)) * latConst / lambda;
    }

    /**
     * starts the blaze condition approcch, values in this method need
     * to be adjusted for the desired system conditions
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        DmdSimulationCore.MetaData meta = new DmdSimulationCore.MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = false;
        
        int lambdaStart = 400;
        int lambdaEnd = 700;
        meta.lambdas = new int[lambdaEnd - lambdaStart + 1];
        for (int i = 0; i <= lambdaEnd-lambdaStart; i++) meta.lambdas[i] = lambdaStart + i;

        meta.latticeConstant = 7.56;
        meta.tiltAngle = -12.0;

        meta.phiInStart = -90;
        meta.phiInEnd = 90;
        meta.inStepSize = 0.2;
        
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
                double alpha = phi;//Math.atan(Math.sqrt(2)*Math.tan(phi*Math.PI/180))*180/Math.PI;
                //System.out.println(inAngle + " " + meta.tiltAngle + " " + waveLength + " " + meta.latticeConstant);
                double n = calcDiffractionOrder(alpha, meta.tiltAngle, waveLength, meta.latticeConstant);
                float value = (float) Math.sqrt(Math.pow(Math.sin(n*Math.PI), 2.0));
                diagonalEpd.set(x, y, value);
                
                outAngle.set(x, y, (float) (-alpha + 2*meta.tiltAngle));
            }
        }
        diagonalEpd.saveAsTiff(meta.outDir + diagonalEpd.getTitle() + ".tif", meta);
        outAngle.saveAsTiff(meta.outDir + outAngle.getTitle() + ".tif", meta);
        
    }
    
}
