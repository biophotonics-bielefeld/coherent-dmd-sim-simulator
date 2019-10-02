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
 * Class for simulating multiple colors
 * @author m.lachetta
 */
public class MultiColCalculator {
    final DmdSimulator.MetaData meta;
    final Dmd dmd;
    
    MultiColCalculator(DmdSimulator.MetaData meta) {
        this.meta = meta;
        double mirrorSize = Math.sqrt(meta.latticeConstant*meta.latticeConstant*meta.fillFactor);
        double gap = meta.latticeConstant-mirrorSize;
        dmd = new Dmd(meta.nrX, meta.nrY, mirrorSize, gap);
    }
    
    public static void main(String[] args) {
        // creating meta data object
        DmdSimulator.MetaData meta = new DmdSimulator.MetaData();
        try {
            if (args.length >= 2 && args.length < 19) {
                meta.outDir = args[0] + "/";
                meta.gpuActive = Boolean.valueOf(args[1]);
                //meta.beamDiameter = 8000;
                meta.lambdas = new int[args.length - 2];
                for (int i = 0; i < args.length-2; i++) meta.lambdas[i] = Integer.parseInt(args[i+2]);

                meta.nrX = 20;
                meta.nrY = 20;

                meta.latticeConstant = 7.56;
                meta.fillFactor = 0.92;
                meta.tiltAngle = 12;
                
                meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

                meta.phiOutStart = -80;
                meta.phiOutEnd = 80;
                meta.thetaOutStart = -80;
                meta.thetaOutEnd = 80;
                meta.outStepSize = 0.1;

                meta.phiInStart = -45;
                meta.phiInEnd = 45;
                meta.thetaInStart = -45;
                meta.thetaInEnd = 45;
                meta.inStepSize = 0.2;

                //meta.bmp = Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\SLM_0,40_1,75_33_wl532_ang0_pha0.bmp");
                meta.bmp = new Image(meta.nrX, meta.nrY);
                //meta.bmp.setAll(1);
            } else if (args.length >= 19) {
                meta.outDir = args[0] + "/";
                meta.gpuActive = Boolean.valueOf(args[1]);
                meta.lambdas = new int[args.length - 18];
                for (int i = 0; i < args.length-18; i++) meta.lambdas[i] = Integer.parseInt(args[i+18]);

                meta.nrX = Integer.parseInt(args[2]);
                meta.nrY = Integer.parseInt(args[3]);

                meta.latticeConstant = Double.parseDouble(args[4]);
                meta.fillFactor = Double.parseDouble(args[5]);
                meta.tiltAngle = Double.parseDouble(args[6]);

                meta.phiOutStart = Double.parseDouble(args[7]);
                meta.phiOutEnd = Double.parseDouble(args[8]);
                meta.thetaOutStart = Double.parseDouble(args[9]);
                meta.thetaOutEnd = Double.parseDouble(args[10]);
                meta.outStepSize = Double.parseDouble(args[11]);

                meta.phiInStart = Double.parseDouble(args[12]);
                meta.phiInEnd = Double.parseDouble(args[13]);
                meta.thetaInStart = Double.parseDouble(args[14]);
                meta.thetaInEnd = Double.parseDouble(args[15]);
                meta.inStepSize = Double.parseDouble(args[16]);

                try {
                    int flipState = Integer.parseInt(args[17]);
                    meta.bmp = new Image(meta.nrX, meta.nrY);
                    meta.bmp.setAll(flipState);
                } catch (NumberFormatException ex) {
                    meta.bmp = Image.readBitmap(args[17]);
                }
            } else {
                String error = "Usage: [0:directory for saving] [1:true or false for gpuActive]"
                        + "[2:lambda0]  [3:lambda1] ..."
                        + "\n or \n"
                        + "[0:directory for saving] [1:true or false for gpuActive]"
                        + "[1:beamDiameter in µm] [2:nrX mirrors] [3:nrYmirrors] "
                        + "[4:single mirrorSize in µm] [5:dmd fill factor] "
                        + "[6:tiltAngle] [7:phiOutStart] [8:phiOutEnd] "
                        + "[9:thetaOutStart] [10:thetaOutEnd] [11:outStepSize] "
                        + "[12:phiInStart] [13:phiInEnd] [14:thetaInStart] "
                        + "[15:thetaInEnd] [16:inStepSize]"
                        + "[17:flipState as bitMap path] [18:lambda0]  [19:lambda1] ...\n"
                        + "all angles in degree";
                System.err.println(error);
                throw new RuntimeException(error);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        
        for(int lambda : meta.lambdas) System.out.println(lambda);
        for(int lambda : meta.lambdas) DmdSimulator.simulateColorFast(meta, lambda, false);
        
    }
}
