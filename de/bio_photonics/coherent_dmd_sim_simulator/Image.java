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
package de.bio_photonics.coherent_dmd_sim_simulator;

import ij.ImagePlus;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.FloatProcessor;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * simple class to handle 2 dimensinal float images
 * @author Mario
 */
public class Image {
    
    final private int width, height;
    private float[] array;
    private ImagePlus ip;
    
    /**
     * Creates a new black (0) image
     * @param width
     * @param height 
     */
    public Image(int width, int height) {
        this.width = width;
        this.height = height;
        array = new float[width * height];
        FloatProcessor fp = new FloatProcessor(width, height, array);
        ip = new ImagePlus("dmd-ray-tracer image", fp);
    }
    
    /**
     * Creates a new image with the given float array
     * @param width
     * @param height
     * @param array 
     */
    public Image(int width, int height, float[] array) {
        if (width * height != array.length) throw new RuntimeException("Missmacht in dimension");
        this.width = width;
        this.height = height;
        this.array = array;
        FloatProcessor fp = new FloatProcessor(width, height, array);
        ip = new ImagePlus("dmd-ray-tracer image", fp);
    }
    
    /**
     * replaces this image array
     * @param img the new image array
     */
    public void set(Image img) {
        if (width * height != img.array.length) throw new RuntimeException("Missmacht in dimension");
        array = img.array;
    }
    
    /**
     * pixel wise image multiplication
     * @param other 
     */
    public void multiply(Image other) {
        if (width != other.width || height != other.height) throw new RuntimeException("Missmacht in dimension");
        for (int i = 0; i<array.length; i++) {
            array[i] *= other.array[i];
        }
    }
    
    /**
     * multiplies every pixel with the value
     * @param value 
     */
    public void multiply(float value) {
        for (int i = 0; i<array.length; i++) {
            array[i] *= value;
        }
    }
    
    /**
     * adds the value to every pixel
     * @param value 
     */
    public void add(float value) {
        for (int i = 0; i<array.length; i++) {
            array[i] += value;
        }
    }
    
    
    /**
     * multiplies two images pixelwise
     * @param i1
     * @param i2
     * @return 
     */
    public static Image multiply(Image i1, Image i2) {
        if (i1.width != i2.width || i1.height != i2.height) throw new RuntimeException("Missmacht in dimension");
        Image img = new Image(i1.width, i1.height);
        for (int i = 0; i<img.array.length; i++) {
            img.array[i] = i1.array[i] * i2.array[i];
        }
        return img;
    }
    
    /**
     * converts the float array of this image to a double array
     * @return 
     */
    public double[][] asDoubleArray() {
        double[][] dArray = new double[height][width];
        for(int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                dArray[y][x] = get(x, y);
            }
        }
        return dArray;
    }
    
    /**
     * logarithmizes this image
     */
    public void log() {
        for (int i = 0; i<array.length; i++) {
            array[i] = (float) Math.log10(array[i]);
        }
    }
    
    /**
     * Sets the title of the connected ImagePlus instance
     * @param title 
     */
    public void setTitle(String title) {
        ip.setTitle(title);
    }
    
    /**
     * 
     * @return title of the connected ImagePlus instance
     */
    public String getTitle() {
        return ip.getTitle();
    }
    
    /**
     * finds the maximum in this image
     * @return {x,y,maxValue}
     */
    public int[] findMax() {
        float max = -Float.MAX_VALUE;
        int iMax = -1;
        for (int i = 0; i<array.length; i++) {
            if (array[i] > max) {
                max = array[i];
                iMax = i;
            }
        }
        int x = iMax % width;
        int y = iMax / width;
        int[] ret = {x,y,(int) max};
        return ret;
    }
    
    /**
     * sets a value in this image
     * @param x
     * @param y
     * @param value 
     */
    public void set(int x, int y, float value) {
        if (x >= width || y >= height) throw new IndexOutOfBoundsException();
        array[width * y + x] = value;
    }
    
    /**
     * sets all pixel in this image
     * @param value 
     */
    public void setAll(float value) {
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                set(x, y, value);
            }
        }
    }
    
    /**
     * 
     * @param x
     * @param y
     * @return value of the x y pixel
     */
    public float get(int x, int y) {
        return array[width * y + x];
    }
    
    /**
     * 
     * @return width of this image
     */
    public int getWidth() {
        return width;
    }
    
    /**
     * 
     * @return height of this image
     */
    public int getHeight() {
        return height;
    }
    
    /**
     * 
     * @return max value of this image
     */
    public float getMaxValue() {
        float max = -Float.MAX_VALUE;
        for (float f : array) {
            if (f > max) max = f;
        }
        return max;
    }
    
    /**
     * updates the floatProcessor of the connected ImagePlus instance
     */
    public void update() {
        FloatProcessor fp = new FloatProcessor(width, height, array);
        ip.setProcessor(fp);
    }
    
    /**
     * 
     * @return the floatProcessor of the connected ImagePlus instance
     */
    public FloatProcessor getFloatProcessor() {
        return new FloatProcessor(width, height, array);
    }
    
    /**
     * shows this image
     */
    public void show() {
        update();
        ip.show();
    }
    
    /**
     * repaint this image
     */
    public void repaint() {
        update();
        ip.repaintWindow();
    }
    
    /**
     * closes this image
     */
    public void close() {
        ip.close();
    }
    
    /**
     * saves this image
     * @param outfile path to save this image
     */
    public void saveAsTiff(String outfile) {
        update();
        new FileSaver(ip).saveAsTiff(outfile);
    }
    
    /**
     * saves this image
     * @param outfile path to save this image
     * @param meta meta object
     */
    public void saveAsTiff(String outfile, DmdSimulationCore.MetaData meta) {
        update();
        ip.setProperty("Info", meta.toString());
        new FileSaver(ip).saveAsTiff(outfile);
    }
    
    /**
     * reads a bitmap, if(pixel==0) 0 else 1
     * @param file file path to the bitmap
     * @return the image of the created bitmap
     */
    public static Image readBitmap(String file) {
        ImagePlus ip = new Opener().openImage(file);
        int[][] intArray = ip.getProcessor().getIntArray();
        int width = intArray.length;
        int height = intArray[0].length;
        Image bmp = new Image(width, height);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (intArray[x][y] == 0) bmp.set(x, y, 0);
                else bmp.set(x, y, 1);
            }
        }
        return bmp;
    }
    
    /**
     * generates an image of concentric black/white circles
     * @param size width & height of this image
     * @return 
     */
    public static Image generateCircles(int size) {
        Image img = new Image(size, size);
        for(int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                int X = x - size/2;
                int Y = y - size/2;
                double r = Math.sqrt(X*X+Y*Y);
                int R = (int) r;
                if (R % 4 < 2) img.set(x, y, 0);
                else img.set(x, y, 1);
            }
        }
        return img;
    }
    
    private static String[] readFile(String path) {
        try {
            File file = new File(path);
            Scanner reader = new Scanner(file);
            ArrayList<String> lines = new ArrayList<>();
            while (reader.hasNextLine()) {
                lines.add(reader.nextLine());
            }
            reader.close();
            String[] lineArray = new String[lines.size()];
            lineArray = lines.toArray(lineArray);
            return lineArray;
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
    
    private static float[][] linesToFloats(String[] lines) {
        int N = lines.length;
        int M = lines[0].split(" ").length;
        float[][] floatArray = new float[M][N];
        for(int n = 0; n < N; n++) {
            String[] numbers = lines[n].split(" ");
            for (int m = 0; m < M; m++) {
                Float value = Float.valueOf(numbers[m]);
                floatArray[m][n] = value;
            }
        }
        return floatArray;
    }
    
    /**
     * converts outputs of the c++ implementation of the ray tracing approach
     * (https://github.com/hsandmeyer/dmd_traycing) into images
     * @param path
     * @return 
     */
    public static Image readImageFromTxtFile(String path) {
        String[] lines = readFile(path);
        float[][] floatArray = linesToFloats(lines);
        float[] values = floatArray[2];
        int width = (int) Math.sqrt(values.length);
        Image img = new Image(width, width);
        for(int x = 0; x < width; x++) {
            for (int y = 0; y < width; y++) {
                img.set(x, y, values[width * x + y]);
            }
        }
        return img;
    }
    
    /**
     * for testing
     * @param args 
     */
    public static void main(String args[]) {
        
        Image img = generateCircles(50);
        img.saveAsTiff("D:\\dmd-simulator-images\\interesting patterns\\circles-50.tif");
        /*
        String folder = "D:\\Gits\\dmd_traycing\\build\\";
        String fileWoEnding = "ray_tracing";
        Image img = readImageFromTxtFile(folder + fileWoEnding + ".txt");
        img.saveAsTiff(folder + fileWoEnding + ".tif");
        fileWoEnding = "one_mirror_ref_image";
        img = readImageFromTxtFile(folder + fileWoEnding + ".txt");
        img.saveAsTiff(folder + fileWoEnding + ".tif");
        fileWoEnding = "analytic";
        img = readImageFromTxtFile(folder + fileWoEnding + ".txt");
        img.saveAsTiff(folder + fileWoEnding + ".tif");
         */
    }
}
