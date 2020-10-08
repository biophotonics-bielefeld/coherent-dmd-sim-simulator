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

/**
 * Class which represents a dmd with a diagonal tilt axis
 * the mathematic model is described in the mathematica file dmd.nb 
 * @author m.lachetta
 */
public class Dmd {
    
    int nrX, nrY;
    double mirrorWidth, mirrorHeight, gapX, gapY, dmdWidth, dmdHeight;
    
    /**
     * Constructor for rectangular mirrors in the dmd
     * @param nrX amount of mirrors in x
     * @param nrY amount of mirrors in y
     * @param mirrorWidth width of a single mirror in µm
     * @param mirrorHeight height of a single mirror in µm
     * @param gapX gap in x between mirrors in µm
     * @param gapY gap in y between mirrors in µm
     */
    public Dmd(int nrX, int nrY, double mirrorWidth,
            double mirrorHeight, double gapX, double gapY) {
        this.nrX = nrX;
        this.nrY = nrY;
        this.mirrorWidth = mirrorWidth;
        this.mirrorHeight = mirrorHeight;
        this.gapX = gapX;
        this.gapY = gapY;
        this.dmdWidth = nrX * (mirrorWidth + gapX) - gapX;
        this.dmdHeight = nrY * (mirrorHeight + gapY) - gapY;
    }
    
    /**
     * Constructor for quadratic mirrors on the dmd
     * @param nrX amount of mirrors in x
     * @param nrY amount of mirrors in y
     * @param mirrorSize width & height of a single mirror in µm
     * @param gap gap in x & y between mirrors in µm
     */
    public Dmd(int nrX, int nrY, double mirrorSize, double gap) {
        this(nrX, nrY, mirrorSize, mirrorSize, gap, gap);
    }
    
    /**
     * Calculates the x coordinate dependent on mx, my, s, t in dmd.nb
     * @param mx mirror in x
     * @param my mirror in y
     * @param tiltAngleInDegree tilt of the mirror in degree
     * @param s s position on the mirror
     * @param t t position on the mirror
     * @return x coordinate in µm
     */
    public double getX(int mx, int my, double tiltAngleInDegree, double s, double t) {
        rangeCheck(mx, my, tiltAngleInDegree, s, t);
        double a = 2 * Math.PI * tiltAngleInDegree / 360;
        double cos = Math.cos(a);
        return (mirrorWidth + gapX) * mx + 0.5 * t *(1 - cos) + s * (0.5 * (1 - cos) + cos) - dmdWidth / 2;
    }
    
    /**
     * Calculates the y coordinate dependent on mx, my, s, t in dmd.nb
     * @param mx mirror in x
     * @param my mirror in y
     * @param tiltAngleInDegree tilt of the mirror in degree
     * @param s s position on the mirror
     * @param t t position on the mirror
     * @return x coordinate in µm
     */
    public double getY(int mx, int my, double tiltAngleInDegree, double s, double t) {
        rangeCheck(mx, my, tiltAngleInDegree, s, t);
        double a = 2 * Math.PI * tiltAngleInDegree / 360;
        double cos = Math.cos(a);
        return (mirrorHeight + gapY) * my + 0.5 * s * (1 - cos) + t * (0.5 * (1 - cos) + cos) - dmdHeight / 2;
    }
    
    /**
     * Calculates the z coordinate dependent on mx, my, s, t in dmd.nb
     * @param mx mirror in x
     * @param my mirror in y
     * @param tiltAngleInDegree tilt of the mirror in degree
     * @param s s position on the mirror
     * @param t t position on the mirror
     * @return z coordinate in µm
     */
    public double getZ(int mx, int my, double tiltAngleInDegree, double s, double t) {
        rangeCheck(mx, my, tiltAngleInDegree, s, t);
        double a = 2 * Math.PI * tiltAngleInDegree / 360;
        double sin = Math.sin(a);
        return (- s * sin + t * sin) / Math.sqrt(2);
    }
    
    /**
     * Calculates the x,y,z coordinates dependent on mx, my, s, t in dmd.nb
     * @param mx mirror in x
     * @param my mirror in y
     * @param tiltAngleInDegree tilt of the mirror in degree
     * @param s s position on the mirror
     * @param t t position on the mirror
     * @return x,y,z coordinates as vector
     */
    public Vector getCoordinates(int mx, int my, double tiltAngleInDegree, double s, double t) {
        double x = getX(mx, my, tiltAngleInDegree, s, t);
        double y = getY(mx, my, tiltAngleInDegree, s, t);
        double z = getZ(mx, my, tiltAngleInDegree, s, t);
        return new Vector(x,y,z);
    }
    
    /**
     * Checks ranges of mx, my, tiltAngle, s, t
     * @param mx mirror in x
     * @param my mirror in y
     * @param tiltAngleInDegree tilt of the mirror in degree
     * @param s s position on the mirror
     * @param t t position on the mirror
     */
    private void rangeCheck(int mx, int my, double tiltAngleInDegree, double s, double t) {
        if (mx < 0 || mx >= nrX)
            throw new RuntimeException("mx has to be between 0 and " + nrX);
        if (my < 0 || my >= nrY)
            throw new RuntimeException("my has to be between 0 and " + nrY);
        if (Math.abs(tiltAngleInDegree) > 90)
            throw new RuntimeException("tiltAngleInDegree has to be between -90 and 90");
        if (s < 0 || s > mirrorWidth)
            throw new RuntimeException("s has to be between 0 and " + mirrorWidth);
        if (t < 0 || t > mirrorHeight)
            throw new RuntimeException("t has to be between 0 and " + mirrorHeight);
    }
    
    /**
     * Creates an image with the topology of the dmd
     * @param nrMirrorsX amount of mirrors in x
     * @param nrMirrorsY amount of mirrors in y
     * @param pixelsPerMicron pixels per micron in x & y
     * @param tilt tilts all mirrors into this angle
     * @return dmd topology image
     */
    public Image getSurfaceView(int nrMirrorsX, int nrMirrorsY, int pixelsPerMicron, double tilt) {
        double[][] tiltAngleInDegree = new double[nrMirrorsY][nrMirrorsX];
        for (int my = 0; my < nrMirrorsY; my++) {
            for (int mx = 0; mx < nrMirrorsX; mx++) {
                tiltAngleInDegree[my][mx] = tilt;
            }
        }
        return getSurfaceView(nrMirrorsX, nrMirrorsY, pixelsPerMicron, tiltAngleInDegree);
    }
    
    /**
     * Creates an image with the topology of the dmd
     * @param nrMirrorsX amount of mirrors in x
     * @param nrMirrorsY amount of mirrors in y
     * @param pixelsPerMicron pixels per micron in x & y
     * @param tiltAngleInDegree tilts the mirror defined by this array [y][x]
     * @return dmd topology image
     */
    public Image getSurfaceView(int nrMirrorsX, int nrMirrorsY, int pixelsPerMicron, double[][] tiltAngleInDegree) {
        int stPerMirror = pixelsPerMicron * 10; // calculates how many s & t per mirror
        double offsetX = dmdWidth / 2; // shift to only get positive x & y values
        double offsetY = dmdHeight / 2; // shift to only get positive x & y values
        double offsetZ = 0;
        int width = (int) (nrMirrorsX * (mirrorWidth + gapX) * pixelsPerMicron);
        int height = (int) (nrMirrorsY * (mirrorHeight + gapY) * pixelsPerMicron);
        Image surface = new Image(width, height);
        surface.setAll(Float.NaN);
        for (int my = 0; my < nrMirrorsY; my++) {
            for (int mx = 0; mx < nrMirrorsX; mx++) {
                for (int ss = 0; ss < stPerMirror; ss++) {
                    double s = mirrorWidth / stPerMirror * ss;
                    for (int tt = 0; tt < stPerMirror; tt++) {
                        double t = mirrorHeight / stPerMirror * tt;
                        double x = getX(mx, my, tiltAngleInDegree[my][mx], s, t);
                        double y = getY(mx, my, tiltAngleInDegree[my][mx], s, t);
                        double z = getZ(mx, my, tiltAngleInDegree[my][mx], s, t);
                        
                        int setX = (int) ((x + offsetX) * pixelsPerMicron);
                        int setY = (int) ((y + offsetY) * pixelsPerMicron);
                        float setZ = (float) (z + offsetZ);
                        try {
                            surface.set(setX, setY, setZ);
                        } catch (IndexOutOfBoundsException ex) {
                            // if calculated x|y value is out of the image boundaries, the z value will be ignorred
                        }
                    }
                }
            }
        }
        return surface;
    }
    
    
    /**
     * for testing shows a topography
     * @param args 
     */
    public static void main(String[] args) {
        int nrX = 1920;
        int nrY = 1080;
        double latticeConstant = 7.56;
        double fillFactor = 0.92;
        double mirrorSize = Math.sqrt(latticeConstant*latticeConstant*fillFactor);
        double gap = latticeConstant-mirrorSize;
        
        Dmd dmd = new Dmd(nrX, nrY, mirrorSize, gap);
        
        
        Image pattern = new Image(100, 100);//Image.readBitmap("C:\\Users\\m.lachetta\\Downloads\\DLP6500_1,35_1,75_33_wl532_ang2_pha0.bmp");
        pattern.multiply(24);
        pattern.add(-12);
        pattern.saveAsTiff("D:\\dmd-simulator-images\\pattern.tif");
        Image surfaceView = dmd.getSurfaceView(1, 1, 1000, pattern.asDoubleArray());
        surfaceView.saveAsTiff("D:\\dmd-simulator-images\\dmd_10.tif");
        
        //surfaceView = dmd.getSurfaceView(40, 40, 25, pattern.asDoubleArray());
        //surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_40.tif");
        
//        surfaceView = dmd.getSurfaceView(2, 2, 500, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_2.tif");
//        
//        surfaceView = dmd.getSurfaceView(5, 5, 200, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_5.tif");
//        
//        surfaceView = dmd.getSurfaceView(10, 10, 100, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_10.tif");
//        
//        surfaceView = dmd.getSurfaceView(20, 20, 50, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_20.tif");
//        
//        surfaceView = dmd.getSurfaceView(100, 100, 10, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_100.tif");
//        
//        surfaceView = dmd.getSurfaceView(200, 200, 5, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_200.tif");
//        
//        surfaceView = dmd.getSurfaceView(500, 500, 2, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_500.tif");
//        
//        surfaceView = dmd.getSurfaceView(1920, 1080, 1, 12);
//        surfaceView.saveAsTiff("G:\\dmd-ray-tracer-imgs\\dmd_full.tif");
    }
}
