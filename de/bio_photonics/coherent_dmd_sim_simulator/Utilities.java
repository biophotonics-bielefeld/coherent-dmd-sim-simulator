/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.coherent_dmd_sim_simulator;

import de.bio_photonics.coherent_dmd_sim_simulator.DmdSimulationCore.MetaData;
import de.bio_photonics.coherent_dmd_sim_simulator.Image;
import ij.IJ;
import ij.gui.GenericDialog;
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Mario
 */
public class Utilities {
    
//    static MetaData askGeneralMetaData() {
//        
//        GenericDialog gd = new GenericDialog("Analytic Phase Shifter");
//        
//        gd.addMessage("Gelenerel Options");
//        //gd.addCheckbox("GPU Support", false);
//        gd.addStringField("Storing Directory", "D:\\dmd-simulator-images\\fiji-plugin-test", 100);
//        
//        
//        gd.addMessage("DMD Options");
//        gd.addNumericField("Nr mirrors X", 50, 0);
//        gd.addNumericField("Nr mirrors Y", 50, 0);
//        gd.addNumericField("Lattice Constant", 7.56, 3, 5, "µm");
//        gd.addNumericField("Fill Factor", 0.92, 3);
//        gd.addNumericField("Mirrors Tilt Angle", 12.0, 2, 5, "°");
//        gd.addStringField("Tilt State BMP", "D:\\dmd-simulator-images\\interesting patterns\\circles-50.bmp", 100);
//        
//        gd.showDialog();
//        
//        MetaData meta = new MetaData();
//        
//        
//        meta.gpuActive = false; //gd.getNextBoolean();
//        meta.outDir = gd.getNextString() + "/";
//        File f = new File(meta.outDir);
//        if (!(f.exists() && f.isDirectory())) {
//            String message = meta.outDir + " is not a directory";
//            IJ.error(message);
//            throw new RuntimeException(message);
//        }
//        
//        
//        meta.nrX = (int) gd.getNextNumber();
//        meta.nrY = (int) gd.getNextNumber();
//        meta.latticeConstant = gd.getNextNumber();
//        meta.fillFactor = gd.getNextNumber();
//        meta.tiltAngle = gd.getNextNumber();
//        String bmpString = gd.getNextString();
//        f = new File(bmpString);
//        if (!(f.exists() && f.isFile())) {
//            String message = bmpString + " is not a file";
//            IJ.error(message);
//            throw new RuntimeException(message);
//        }
//        meta.bmp = Image.readBitmap(bmpString);
//        
//        return meta;
//    }
    
    public static Map infoStringToMap(String info) {
        Map<String, String> map = new HashMap<>();
        String[] split = info.split("\n");
        for (int i = 0; i < split.length; i++) {
            String[] kv = split[i].split(": ");
            if (kv.length == 2) map.put(kv[0], kv[1]);
        }
        PrintStream out = System.out;
        return map;
    }
    
}
