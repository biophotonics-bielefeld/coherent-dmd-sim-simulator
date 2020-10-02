/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.fiji_plugin;

import de.bio_photonics.coherent_dmd_sim_simulator.AnalyticDiagonalCalculator;
import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 *
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
                //System.out.println(inAngle + " " + meta.tiltAngle + " " + waveLength + " " + meta.latticeConstant);
                double n = AnalyticDiagonalCalculator.calcDiffractionOrder(alpha, meta.tiltAngle, waveLength, meta.latticeConstant);
                float value = (float) Math.sqrt(Math.pow(Math.sin(n*Math.PI), 2.0));
                //System.out.println(x + " " + y);
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
