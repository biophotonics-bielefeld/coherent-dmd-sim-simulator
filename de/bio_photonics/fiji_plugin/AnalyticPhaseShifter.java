/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.fiji_plugin;

import de.bio_photonics.coherent_dmd_sim_simulator.AccurateSimulator;
import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore.MetaData;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.io.File;

/**
 *
 * @author Mario
 */
public class AnalyticPhaseShifter implements PlugIn {
    
    MetaData meta;

    @Override
    public void run(String string) {
        
        GenericDialog gd = new GenericDialog("Analytic Phase Shifter");
        
        gd.addMessage("Gelenerel Options");
        //gd.addCheckbox("GPU Support", false);
        gd.addStringField("Storing Directory", "D:\\dmd-simulator-images\\fiji-plugin-test", 100);
        gd.addMessage("DMD Options");
        gd.addNumericField("Nr mirrors X", 50, 0);
        gd.addNumericField("Nr mirrors Y", 50, 0);
        gd.addNumericField("Lattice Constant", 7.56, 3, 5, "µm");
        gd.addNumericField("Fill Factor", 0.92, 3);
        gd.addNumericField("Mirrors Tilt Angle", 12.0, 2, 5, "°");
        gd.addStringField("Tilt State BMP", "D:\\dmd-simulator-images\\interesting patterns\\circles-50.bmp", 100);
        
        gd.addMessage("Incident Options");
        gd.addNumericField("Wavelength", 532, 0, 5, "nm");
        gd.addNumericField("Beam Size Multiplier", 0.5, 2);
        gd.addNumericField("Phi In", -21, 1, 5, "°");
        gd.addNumericField("Theta In", 21, 1, 5, "°");
        
        gd.addMessage("Diffracted Options");
        gd.addNumericField("Phi Out Start", -15, 1, 5, "°");
        gd.addNumericField("Phi Out End", 15, 1, 5, "°");
        gd.addNumericField("Theta Out Start", -15, 1, 5, "°");
        gd.addNumericField("Theta Out End", 15, 1, 5, "°");
        gd.addNumericField("Out Step Size", 0.02, 3, 5, "°");
        
        gd.showDialog();
        
        meta = new MetaData();
        
        meta.gpuActive = false; //gd.getNextBoolean();
        meta.outDir = gd.getNextString() + "/";
        File f = new File(meta.outDir);
        if (!(f.exists() && f.isDirectory())) {
            String message = meta.outDir + " is not a directory";
            IJ.error(message);
            throw new RuntimeException(message);
        }
        
        
        meta.nrX = (int) gd.getNextNumber();
        meta.nrY = (int) gd.getNextNumber();
        meta.latticeConstant = gd.getNextNumber();
        meta.fillFactor = gd.getNextNumber();
        meta.tiltAngle = gd.getNextNumber();
        String bmpString = gd.getNextString();
        f = new File(bmpString);
        if (!(f.exists() && f.isFile())) {
            String message = bmpString + " is not a file";
            IJ.error(message);
            throw new RuntimeException(message);
        }
        meta.bmp = Image.readBitmap(bmpString);
        
        meta.lambdas = new int[]{(int) gd.getNextNumber()};
        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant * gd.getNextNumber());
        meta.phiInStart = gd.getNextNumber();
        meta.phiInEnd = meta.phiInStart + 1;
        meta.thetaInStart = gd.getNextNumber();
        meta.thetaInEnd = meta.thetaInStart + 1;
        meta.inStepSize = 1;
        
        meta.phiOutStart = gd.getNextNumber();
        meta.phiOutEnd = gd.getNextNumber();
        meta.thetaOutStart = gd.getNextNumber();
        meta.thetaOutEnd = gd.getNextNumber();
        meta.outStepSize = gd.getNextNumber();
        
        IJ.log("Preparing simulation...");
        AccurateSimulator as = new AccurateSimulator(meta, true);
        IJ.log("Starting simulation...");
        as.simulate();
        IJ.log("Finnished simulation.");
    }
    
    public static void main(String[] args) {
        AnalyticPhaseShifter a  = new AnalyticPhaseShifter();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	a.run("");
    }
    
}
