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
package de.bio_photonics.fiji_plugin;

import de.bio_photonics.coherent_dmd_sim_simulator.BlazeConditionApproachCalculator;
import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * implementation of the blaze condition approach as fiji plugin
 * @author Mario
 */
public class BlazeConditionApproach implements PlugIn {
    
    DmdSimulationCore.MetaData meta;

    @Override
    public void run(String string) {
        
        GenericDialog gd1 = new GenericDialog("Blaze Condition Approach");
        
        gd1.addMessage("General Options");
        //gd.addCheckbox("GPU Support", false);
        gd1.addStringField("Storage Directory", "D:\\dmd-simulator-images\\fiji-plugin-test", 100);
        gd1.addMessage("DMD Options");
        gd1.addNumericField("Lattice Constant", 7.56, 3, 5, "µm");
        gd1.addNumericField("All Mirrors Tilt Angle", -12.0, 2, 5, "°");
        
        gd1.addMessage("Incidence Parameters");
        gd1.addNumericField("Wavelength Start", 400, 0, 5, "nm");
        gd1.addNumericField("Wavelength End", 700, 0, 5, "nm");
        gd1.addNumericField("Alpha In Start", -90, 1, 5, "°");
        gd1.addNumericField("Alpha In End", 90, 1, 5, "°");
        gd1.addNumericField("In Step Size", 0.2, 1, 5, "°");
        
        gd1.showDialog();
        if (gd1.wasCanceled()) return;
        
        meta = new DmdSimulationCore.MetaData();
        
        meta.gpuActive = false;
        meta.outDir = gd1.getNextString() + "/";
        meta.latticeConstant = gd1.getNextNumber();
        meta.tiltAngle = gd1.getNextNumber();
        int lambdaStart = (int) gd1.getNextNumber();
        int lambdaEnd = (int) gd1.getNextNumber();
        meta.lambdas = new int[lambdaEnd - lambdaStart + 1];
        for (int i = 0; i <= lambdaEnd-lambdaStart; i++) meta.lambdas[i] = lambdaStart + i;
        meta.phiInStart = gd1.getNextNumber();
        meta.phiInEnd = gd1.getNextNumber();
        meta.inStepSize = gd1.getNextNumber();
        
        
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
                double n = BlazeConditionApproachCalculator.calcDiffractionOrder(alpha, meta.tiltAngle, waveLength, meta.latticeConstant);
                float value = (float) Math.sqrt(Math.pow(Math.sin(n*Math.PI), 2.0));
                diagonalEpd.set(x, y, value);
                
                outAngle.set(x, y, (float) (-alpha + 2*meta.tiltAngle));
            }
        }
        diagonalEpd.saveAsTiff(meta.outDir + diagonalEpd.getTitle() + ".tif", meta);
        outAngle.saveAsTiff(meta.outDir + outAngle.getTitle() + ".tif", meta);
        IJ.log("Saved at: " + meta.outDir);
    }
    
    public static void main(String[] args) {
        BlazeConditionApproach a  = new BlazeConditionApproach();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	a.run("");
    }
}
