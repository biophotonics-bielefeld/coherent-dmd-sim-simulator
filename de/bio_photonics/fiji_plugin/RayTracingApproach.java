/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.fiji_plugin;

import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import de.bio_photonics.coherent_dmd_sim_simulator.RayTracingSimulator;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.io.File;

/**
 *
 * @author Mario
 */
public class RayTracingApproach implements PlugIn {
    
    DmdSimulationCore.MetaData meta;

    @Override
    public void run(String string) {
        
        GenericDialog gd1 = new GenericDialog("Ray Tracing Approach");
        
        gd1.addMessage("General Options");
        //gd.addCheckbox("GPU Support", false);
        gd1.addStringField("Storage Directory", "D:\\dmd-simulator-images\\fiji-plugin-test", 100);
        gd1.addMessage("DMD Options");
        gd1.addNumericField("No. mirrors X", 50, 0);
        gd1.addNumericField("No. mirrors Y", 50, 0);
        gd1.addNumericField("Lattice Constant", 7.56, 3, 5, "µm");
        gd1.addNumericField("Fill Factor", 0.92, 3);
        gd1.addNumericField("Mirror Tilt Angle", 12.0, 2, 5, "°");
        gd1.addStringField("Tilt State Image (*.bmp)", "D:\\dmd-simulator-images\\interesting patterns\\circles-50.bmp", 100);
        
        gd1.addMessage("Incidence Parameters");
        gd1.addNumericField("Wavelength", 532, 0, 5, "nm");
        gd1.addNumericField("Beam Size Multiplier", 0.5, 2);
        gd1.addNumericField("Phi In", -21, 1, 5, "°");
        gd1.addNumericField("Theta In", 21, 1, 5, "°");
        gd1.addNumericField("Amount Of Rays", 100000, 0);
        
        gd1.addMessage("Diffraction Parameters");
        gd1.addNumericField("Phi Out Start", -15, 1, 5, "°");
        gd1.addNumericField("Phi Out End", 15, 1, 5, "°");
        gd1.addNumericField("Theta Out Start", -15, 1, 5, "°");
        gd1.addNumericField("Theta Out End", 15, 1, 5, "°");
        gd1.addNumericField("Out Step Size", 0.02, 3, 5, "°");
        
        gd1.showDialog();
        if (gd1.wasCanceled()) return;
        
        meta = new DmdSimulationCore.MetaData();
        
        meta.gpuActive = false; //gd.getNextBoolean();
        meta.outDir = gd1.getNextString() + "/";
        File f = new File(meta.outDir);
        if (!(f.exists() && f.isDirectory())) {
            String message = meta.outDir + " is not a directory";
            IJ.error(message);
            throw new RuntimeException(message);
        }
        
        
        meta.nrX = (int) gd1.getNextNumber();
        meta.nrY = (int) gd1.getNextNumber();
        meta.latticeConstant = gd1.getNextNumber();
        meta.fillFactor = gd1.getNextNumber();
        meta.tiltAngle = gd1.getNextNumber();
        String bmpString = gd1.getNextString();
        f = new File(bmpString);
        if (!(f.exists() && f.isFile())) {
            String message = bmpString + " is not a file";
            IJ.error(message);
            throw new RuntimeException(message);
        }
        meta.bmp = Image.readBitmap(bmpString);
        
        meta.lambdas = new int[]{(int) gd1.getNextNumber()};
        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant * gd1.getNextNumber());
        meta.phiInStart = gd1.getNextNumber();
        meta.phiInEnd = meta.phiInStart + 1;
        meta.thetaInStart = gd1.getNextNumber();
        meta.thetaInEnd = meta.thetaInStart + 1;
        meta.inStepSize = 1;
        int nrRays = (int) gd1.getNextNumber();
        
        meta.phiOutStart = gd1.getNextNumber();
        meta.phiOutEnd = gd1.getNextNumber();
        meta.thetaOutStart = gd1.getNextNumber();
        meta.thetaOutEnd = gd1.getNextNumber();
        meta.outStepSize = gd1.getNextNumber();
        
        IJ.log("Preparing simulation...");
        RayTracingSimulator rts = new RayTracingSimulator(meta, nrRays);
        IJ.log("Starting simulation...");
        rts.simulate();
        IJ.log("Finnished simulation.");
    }
    
    public static void main(String[] args) {
        RayTracingApproach a  = new RayTracingApproach();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	a.run("");
    }
    
}
