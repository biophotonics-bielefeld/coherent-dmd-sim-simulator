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

import de.bio_photonics.coherent_dmd_sim_simulator.AnalyticPhaseShiftingApproachSimulator;
import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore.MetaData;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.io.File;

/**
 * implementation of the analytic phase shifting approach as fiji plugin
 * @author Mario
 */
public class AnalyticPhaseShifter implements PlugIn {
    
    MetaData meta;

    @Override
    public void run(String string) {
        
        GenericDialog gd = new GenericDialog("Analytic Phase Shifting Approach");
        
        gd.addMessage("General Options");
        //gd.addCheckbox("GPU Support", false);
        gd.addStringField("Storage Directory", "D:\\dmd-simulator-images\\fiji-plugin-test", 100);
        gd.addMessage("DMD Options");
        gd.addNumericField("No. mirrors X", 50, 0);
        gd.addNumericField("No. mirrors Y", 50, 0);
        gd.addNumericField("Lattice Constant", 7.56, 3, 5, "µm");
        gd.addNumericField("Fill Factor", 0.92, 3);
        gd.addNumericField("Mirror Tilt Angle", 12.0, 2, 5, "°");
        gd.addStringField("Tilt State Image (*.bmp)", "D:\\dmd-simulator-images\\interesting patterns\\circles-50.bmp", 100);
        
        gd.addMessage("Incidence Parameters");
        gd.addNumericField("Wavelength", 532, 0, 5, "nm");
        gd.addNumericField("Beam Size Multiplier", 0.5, 2);
        gd.addNumericField("Phi In", -21, 1, 5, "°");
        gd.addNumericField("Theta In", 21, 1, 5, "°");
        
        gd.addMessage("Diffraction Parameters");
        gd.addNumericField("Phi Out Start", -15, 1, 5, "°");
        gd.addNumericField("Phi Out End", 15, 1, 5, "°");
        gd.addNumericField("Theta Out Start", -15, 1, 5, "°");
        gd.addNumericField("Theta Out End", 15, 1, 5, "°");
        gd.addNumericField("Out Step Size", 0.02, 3, 5, "°");
        
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
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
        
        meta.lambdas = new double[]{ gd.getNextNumber()};
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
        AnalyticPhaseShiftingApproachSimulator as = new AnalyticPhaseShiftingApproachSimulator(meta, true);
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
