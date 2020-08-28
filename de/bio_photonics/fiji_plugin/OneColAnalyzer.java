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
package de.bio_photonics.fiji_plugin;

import de.bio_photonics.coherent_dmd_sim_simulator.Utilities;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.util.Map;

/**
 *
 * @author m.lachetta
 */
public class OneColAnalyzer implements PlugIn {

    @Override
    public void run(String arg) {
        GenericDialog gd = new GenericDialog("OneColAnalyzer");
        gd.addMessage("Choose images");
        String[] openWindows = ij.WindowManager.getImageTitles();
        if (openWindows.length < 3) {
            IJ.error("Need min 3 open images to run this PlugIn");
            return;
        }
        gd.addChoice("Channel 1 - intensity max X", openWindows, openWindows[1]);
        gd.addChoice("Channel 1 - intensity max Y", openWindows, openWindows[2]);
        gd.addChoice("Channel 1 - envelope peak distance", openWindows, openWindows[0]);
        gd.addNumericField("max envelope peak disstance", 0.1, 4, 8, "°");
        
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
        ImagePlus phiOut1 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus thetaOut1 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus epd1 = ij.WindowManager.getImage(gd.getNextChoice());
        double oneColorThreshold = gd.getNextNumber();
        
        Map<String, String> phiOut1Meta = Utilities.infoStringToMap(phiOut1.getInfoProperty());
        Map<String, String> thetaOut1Meta = Utilities.infoStringToMap(thetaOut1.getInfoProperty());
        Map<String, String> epd1Meta = Utilities.infoStringToMap(epd1.getInfoProperty());
        
        double inStepSize;
        double phiInStart;
        double thetaInStart;
        double outStepSize;
        double phiOutStart;
        double phiOutEnd;
        double thetaOutStart;
        double thetaOutEnd;
        int lambda;
        try {
            phiInStart = Double.parseDouble(phiOut1Meta.get("phiInStart")); // not checked
            thetaInStart = Double.parseDouble(phiOut1Meta.get("thetaInStart")); // not checked
            outStepSize = Double.parseDouble(phiOut1Meta.get("outStepSize"));
            double outStepSize2 = Double.parseDouble(thetaOut1Meta.get("outStepSize"));
            double outStepSize3 = Double.parseDouble(epd1Meta.get("outStepSize"));
            if (outStepSize != outStepSize2 || outStepSize != outStepSize3) {
                IJ.error("outStepSize missmatch");
                return;
            }
            inStepSize = Double.parseDouble(phiOut1Meta.get("inStepSize"));
            double inStepSize2 = Double.parseDouble(thetaOut1Meta.get("inStepSize"));
            double inStepSize3 = Double.parseDouble(epd1Meta.get("inStepSize"));
            if (inStepSize != inStepSize2 || inStepSize != inStepSize3) {
                IJ.error("inStepSize missmatch");
                return;
            }
            phiOutStart = Double.parseDouble(phiOut1Meta.get("phiOutStart")); // not checked
            phiOutEnd = Double.parseDouble(phiOut1Meta.get("phiOutEnd")); // not checked
            thetaOutStart = Double.parseDouble(phiOut1Meta.get("thetaOutStart")); // not checked
            thetaOutEnd = Double.parseDouble(phiOut1Meta.get("thetaOutEnd")); // not checked
            lambda = Integer.parseInt(phiOut1Meta.get("lambdas").split(" ")[0]); // not checked
        } catch (NumberFormatException ex) {
            IJ.error("stepSize exception: " + ex);
            return;
        }
        oneColorThreshold /= outStepSize;
        int width = phiOut1.getWidth();
        int height = phiOut1.getHeight();
        if (thetaOut1.getWidth() != width || thetaOut1.getHeight() != height) {
            IJ.error("Dimension missmatch Channel 1 - X with Channel 1 - Y");
            return;
        }
        
        int outWidth = (int) ((phiOutEnd - phiOutStart) / outStepSize);
        int outHeight = (int) ((thetaOutEnd - thetaOutStart) / outStepSize);
        Image outImg = new Image(outWidth, outHeight);
        outImg.setTitle(lambda + "_out_angles");
        outImg.show();
        float[][] epd1Array = epd1.getProcessor().getFloatArray();
        float[][] phiOutArray = phiOut1.getProcessor().getFloatArray();
        float[][] thetaOutArray = thetaOut1.getProcessor().getFloatArray();
        for (int th = 0; th < height; th++) for (int ph = 0; ph < width; ph++) {
            if (epd1Array[ph][th] <= oneColorThreshold) {
                float phiOut = phiOutArray[ph][th];
                float thetaOut = thetaOutArray[ph][th];
                int phOut = (int) ((phiOut - phiOutStart) / outStepSize);
                int thOut = (int) ((thetaOut - thetaOutStart) / outStepSize);
                float temp = outImg.get(phOut, thOut);
                outImg.set(phOut, thOut, temp+1);
                outImg.repaint();
                //IJ.log("Channel 1: " + ph + " " + th);
            }
        }
    }
    
    /** main method */
    public static void main( String [] args ) {

	
	OneColAnalyzer a  = new OneColAnalyzer();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	a.run("");
    }

}
