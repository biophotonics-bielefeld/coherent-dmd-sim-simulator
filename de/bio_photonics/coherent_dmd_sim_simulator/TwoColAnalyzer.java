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

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author m.lachetta
 */
public class TwoColAnalyzer implements PlugIn {

    @Override
    public void run(String arg) {
        GenericDialog gd = new GenericDialog("Analying ray tracer results");
        gd.addMessage("Choose images");
        String[] openWindows = ij.WindowManager.getImageTitles();
        if (openWindows.length < 6) {
            IJ.error("Need min 6 open images to run this PlugIn");
            return;
        }
        gd.addChoice("Channel 1 - intensity max X", openWindows, openWindows[1]);
        gd.addChoice("Channel 1 - intensity max Y", openWindows, openWindows[2]);
        gd.addChoice("Channel 1 - envelope peak distance", openWindows, openWindows[0]);
        gd.addChoice("Channel 2 - intensity max X", openWindows, openWindows[4]);
        gd.addChoice("Channel 2 - intensity max Y", openWindows, openWindows[5]);
        gd.addChoice("Channel 2 - envelope peak distance", openWindows, openWindows[3]);
        gd.addNumericField("max envelope peak disstance", 0.2, 4, 8, "°");
        gd.addNumericField("max two color disstance", 0.001, 4, 8, "°");
        
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
        ImagePlus phiOut1 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus thetaOut1 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus epd1 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus phiOut2 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus thetaOut2 = ij.WindowManager.getImage(gd.getNextChoice());
        ImagePlus epd2 = ij.WindowManager.getImage(gd.getNextChoice());
        double oneColorThreshold = gd.getNextNumber();
        double towColorThreshold = gd.getNextNumber();
        
        Map<String, String> phiOut1Meta = infoStringToMap(phiOut1.getInfoProperty());
        Map<String, String> thetaOut1Meta = infoStringToMap(thetaOut1.getInfoProperty());
        Map<String, String> epd1Meta = infoStringToMap(epd1.getInfoProperty());
        Map<String, String> phiOut2Meta = infoStringToMap(phiOut2.getInfoProperty());
        Map<String, String> thetaOut2Meta = infoStringToMap(thetaOut2.getInfoProperty());
        Map<String, String> epd2Meta = infoStringToMap(epd2.getInfoProperty());
        
        double inStepSize;
        double phiInStart;
        double thetaInStart;
        double outStepSize;
        try {
            phiInStart = Double.parseDouble(phiOut1Meta.get("phiInStart")); // not checked
            thetaInStart = Double.parseDouble(phiOut1Meta.get("thetaInStart")); // not checked
            outStepSize = Double.parseDouble(phiOut1Meta.get("outStepSize"));
            double outStepSize2 = Double.parseDouble(thetaOut1Meta.get("outStepSize"));
            double outStepSize3 = Double.parseDouble(epd1Meta.get("outStepSize"));
            double outStepSize4 = Double.parseDouble(phiOut2Meta.get("outStepSize"));
            double outStepSize5 = Double.parseDouble(thetaOut2Meta.get("outStepSize"));
            double outStepSize6 = Double.parseDouble(epd2Meta.get("outStepSize"));
            if (outStepSize != outStepSize2 || outStepSize != outStepSize3
                    || outStepSize != outStepSize4 || outStepSize != outStepSize5
                    || outStepSize != outStepSize6) {
                IJ.error("outStepSize missmatch");
                return;
            }
            inStepSize = Double.parseDouble(phiOut1Meta.get("inStepSize"));
            double inStepSize2 = Double.parseDouble(thetaOut1Meta.get("inStepSize"));
            double inStepSize3 = Double.parseDouble(epd1Meta.get("inStepSize"));
            double inStepSize4 = Double.parseDouble(phiOut2Meta.get("inStepSize"));
            double inStepSize5 = Double.parseDouble(thetaOut2Meta.get("inStepSize"));
            double inStepSize6 = Double.parseDouble(epd2Meta.get("inStepSize"));
            if (inStepSize != inStepSize2 || inStepSize != inStepSize3
                    || inStepSize != inStepSize4 || inStepSize != inStepSize5
                    || inStepSize != inStepSize6) {
                IJ.error("inStepSize missmatch");
                return;
            }
        } catch (NumberFormatException ex) {
            IJ.error("stepSize exception: " + ex);
            return;
        }
        
        int width = phiOut1.getWidth();
        int height = phiOut1.getHeight();
        if (thetaOut1.getWidth() != width || thetaOut1.getHeight() != height) {
            IJ.error("Dimension missmatch Channel 1 - X with Channel 1 - Y");
            return;
        }
        if (phiOut2.getWidth() != width || phiOut2.getHeight() != height) {
            IJ.error("Dimension missmatch Channel 1 - X with Channel 2 - X");
            return;
        }
        if (thetaOut2.getWidth() != width || thetaOut2.getHeight() != height) {
            IJ.error("Dimension missmatch Channel 1 - X with Channel 2 - Y");
            return;
        }
        
        float[][] epd1Array = epd1.getProcessor().getFloatArray();
        float[][] epd2Array = epd2.getProcessor().getFloatArray();
        List<int[]> oneColCanidates1 = new ArrayList<>();
        List<int[]> oneColCanidates2 = new ArrayList<>();
        for (int th = 0; th < height; th++) for (int ph = 0; ph < width; ph++) {
            if (epd1Array[ph][th] <= oneColorThreshold) {
                oneColCanidates1.add(new int[] {ph, th});
                //IJ.log("Channel 1: " + ph + " " + th);
            }
            if (epd2Array[ph][th] <= oneColorThreshold) {
                oneColCanidates2.add(new int[] {ph, th});
                //IJ.log("Channel 2: " + ph + " " + th);
            }
        }
        IJ.log(oneColCanidates1.size() + " and " + oneColCanidates2.size() + " oneColorCanidates");
        float[][] phiOut1Array = phiOut1.getProcessor().getFloatArray();
        float[][] thetaOut1Array = thetaOut1.getProcessor().getFloatArray();
        float[][] phiOut2Array = phiOut2.getProcessor().getFloatArray();
        float[][] thetaOut2Array = thetaOut2.getProcessor().getFloatArray();
        
        //for (int[] c1 : oneColCanidates1) IJ.log("c1: " + c1[0] + " " + c1[1]);
        //for (int[] c2 : oneColCanidates2) IJ.log("c2: " + c2[0] + " " + c2[1]);
        IJ.log("List of two color canidates \n"
                + "phiIn1 \t thetaIn1 \t phiIn2 \t thetaIn2 \t phiOut1Array[ph1][th1] \t thetaOut1Array[ph1][th1] \t phiOut2Array[ph2][th2] \t thetaOut2Array[ph2][th2] \t twoColDistance \t epd1Array[ph1][th1] \t epd2Array[ph2][th2] \t distanceFromZero1 \t distanceFromZero2");
        
        List<int[]> twoColCanidates = new ArrayList<>();
        for (int[] c1 : oneColCanidates1) {
            int ph1 = c1[0], th1 = c1[1];
//            if (ph1 == 54 && th1 == 270) {
//                IJ.log("C1 case found!");
//            }
            for (int[] c2 : oneColCanidates2) {
                int ph2 = c2[0], th2 = c2[1];
                float pDistance = phiOut2Array[ph2][th2]-phiOut1Array[ph1][th1];
                float tDistance = thetaOut2Array[ph2][th2]-thetaOut1Array[ph1][th1];
                double distanceFromZero1 = Math.sqrt(phiOut1Array[ph1][th1]*phiOut1Array[ph1][th1] + thetaOut1Array[ph1][th1]*thetaOut1Array[ph1][th1]);
                double distanceFromZero2 = Math.sqrt(phiOut2Array[ph2][th2]*phiOut2Array[ph2][th2] + thetaOut2Array[ph2][th2]*thetaOut2Array[ph2][th2]);
                double twoColDistance = Math.sqrt(pDistance*pDistance + tDistance*tDistance);
                //IJ.log(ph1 + " " + th1 + " " + ph2 + " " + th2 + " " + distance);
                if (twoColDistance < towColorThreshold) {
                    double phiIn1 = phiInStart + ph1 * inStepSize;
                    double thetaIn1 = thetaInStart + th1 * inStepSize; 
                    double phiIn2 = phiInStart + ph2 * inStepSize;
                    double thetaIn2 = thetaInStart + th2 * inStepSize; 
                    //IJ.log("" + epd2Array[ph2][th2]);
                    IJ.log(phiIn1 + " \t " + thetaIn1 + " \t " + phiIn2 + " \t " + thetaIn2 + " \t " + phiOut1Array[ph1][th1] + " \t " + thetaOut1Array[ph1][th1] + " \t " + phiOut2Array[ph2][th2] + " \t " + thetaOut2Array[ph2][th2] + " \t " + twoColDistance + " \t " + epd1Array[ph1][th1] + " \t " + epd2Array[ph2][th2] + " \t " + distanceFromZero1 + " \t " + distanceFromZero2);
                    //IJ.log("Found: " + ph1 + " " + th1 + " " + ph2 + " " + th2 + " " + distance);
                    twoColCanidates.add(new int[] {ph1, th1, ph2, th2});
                }
//                if (ph1 == 54 && th1 == 270) {
//                    IJ.log("c2: " + c2[0] + " " + c2[1]);
//                }
//                if (ph2 == 213 && th2 == 110) {
//                    IJ.log("C2 case found!");
//                }
//                if (ph1 == 54 && th1 == 270 && ph2 == 213 && th2 == 110) {
//                    IJ.log("C1 & C2 case found!");
//                }
            }
        }
        IJ.log("Number of twoColCanidates: " + twoColCanidates.size());
//        IJ.log(phiOut1Array[54][270] + " " + thetaOut1Array[54][270]);
//        IJ.log(phiOut2Array[213][110] + " " + thetaOut2Array[213][110]);

    }
    
    private static Map infoStringToMap(String info) {
        Map<String, String> map = new HashMap<>();
        String[] split = info.split("\n");
        for (int i = 0; i < split.length; i++) {
            String[] kv = split[i].split(": ");
            if (kv.length == 2) map.put(kv[0], kv[1]);
        }
        return map;
    }
    
    /** main method */
    public static void main( String [] args ) {

	
	TwoColAnalyzer a  = new TwoColAnalyzer();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	a.run("");

	/*
	double gratMin = Double.parseDouble( args[0] );
	double gratMax = Double.parseDouble( args[1] );
	int nrPhases   = Integer.parseInt( args[2] );
	int nrDirs     = Integer.parseInt( args[3] );

	gs.calculate(gratMin, gratMax, nrPhases, 
	    nrDirs, 3./180*Math.PI, 20, 0.02, false, 400); */
    }

}
