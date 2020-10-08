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

import de.bio_photonics.coherent_dmd_sim_simulator.CoarseSimulator;
import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.io.File;

/**
 * implementation of the grating approach as fiji plugin
 * @author Mario
 */
public class GratingApproach implements PlugIn {
    
    DmdSimulationCore.MetaData meta;

    @Override
    public void run(String string) {
        
        GenericDialog gd1 = new GenericDialog("Grating Approach");
        
        gd1.addMessage("General Options");
        //gd.addCheckbox("GPU Support", false);
        gd1.addStringField("Storage Directory", "D:\\dmd-simulator-images\\fiji-plugin-test", 100);
        gd1.addMessage("DMD Options");
        gd1.addNumericField("No. mirrors X", 50, 0);
        gd1.addNumericField("No. mirrors Y", 50, 0);
        gd1.addNumericField("Lattice Constant", 7.56, 3, 5, "µm");
        gd1.addNumericField("Fill Factor", 0.92, 3);
        gd1.addNumericField("All Mirrors Tilt Angle", -12.0, 2, 5, "°");
        
        gd1.addMessage("Incidence Parameters");
        gd1.addNumericField("Wavelength", 532, 0, 5, "nm");
        gd1.addNumericField("Beam Size Multiplier", 0.5, 2);
        gd1.addNumericField("Phi In Start", -60, 1, 5, "°");
        gd1.addNumericField("Phi In End", 60, 1, 5, "°");
        gd1.addNumericField("Theta In Start", -60, 1, 5, "°");
        gd1.addNumericField("Theta In End", 60, 1, 5, "°");
        gd1.addNumericField("In Step Size", 0.2, 1, 5, "°");
        
        gd1.addMessage("Diffraction Parameters");
        gd1.addNumericField("Phi Out Start", -80, 1, 5, "°");
        gd1.addNumericField("Phi Out End", 80, 1, 5, "°");
        gd1.addNumericField("Theta Out Start", -80, 1, 5, "°");
        gd1.addNumericField("Theta Out End", 80, 1, 5, "°");
        gd1.addNumericField("Out Step Size", 0.1, 3, 5, "°");
        
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
        meta.bmp = new Image(meta.nrX, meta.nrY);
        
        meta.lambdas = new int[]{(int) gd1.getNextNumber()};
        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant * gd1.getNextNumber());
        meta.phiInStart = gd1.getNextNumber();
        meta.phiInEnd = gd1.getNextNumber();
        meta.thetaInStart = gd1.getNextNumber();
        meta.thetaInEnd = gd1.getNextNumber();
        meta.inStepSize = gd1.getNextNumber();
        
        if (meta.phiInStart + meta.inStepSize > meta.phiInEnd) meta.phiInEnd = meta.phiInStart + meta.inStepSize;
        if (meta.thetaInStart + meta.inStepSize > meta.thetaInEnd) meta.thetaInEnd = meta.thetaInStart + meta.inStepSize;
        
        meta.phiOutStart = gd1.getNextNumber();
        meta.phiOutEnd = gd1.getNextNumber();
        meta.thetaOutStart = gd1.getNextNumber();
        meta.thetaOutEnd = gd1.getNextNumber();
        meta.outStepSize = gd1.getNextNumber();
        
        double widthOut = (meta.phiOutEnd - meta.phiOutStart) / meta.outStepSize;
        double heightOut = (meta.thetaOutEnd - meta.thetaOutStart) / meta.outStepSize;
        double widthIn = (meta.phiInEnd - meta.phiInStart) / meta.inStepSize;
        double heightIn = (meta.thetaInEnd - meta.thetaInStart) / meta.inStepSize;
        int bytesPerPixel = 4;
        double reqRam = (int) (widthIn * heightIn * widthOut * heightOut * bytesPerPixel * 0.0000011) * 0.001;
        
        GenericDialog gd2 = new GenericDialog("Save diffraction images or analyzing envelope to bightes diffraction order only?");
        gd2.addMessage("Saving all diffraction images requires about " + reqRam + " Gb of RAM");
        gd2.enableYesNoCancel("Yes, save diffraction images", "No,  analyzing envelope to bightes diffraction order only");
        gd2.showDialog();
        boolean saveRaws;
        if (gd2.wasCanceled()) return;
        else if (gd2.wasOKed()) saveRaws = true;
        else saveRaws = false;
        
        IJ.log("Preparing simulation...");
        CoarseSimulator cs = new CoarseSimulator(meta, true, saveRaws);
        IJ.log("Starting simulation...");
        cs.simulate();
        IJ.log("Finnished simulation.");
    }
    
    public static void main(String[] args) {
        GratingApproach a  = new GratingApproach();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	a.run("");
    }
    
}
