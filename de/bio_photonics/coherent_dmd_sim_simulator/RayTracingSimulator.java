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

import ij.IJ;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * class which implements the ray tracing approach
 * @author Mario
 */
public class RayTracingSimulator extends AbstractSimulator {
    
    final int nrRays;
    
    public RayTracingSimulator(DmdSimulationCore.MetaData meta, int nrRays) {
        super(meta);
        this.nrRays = nrRays;
    }
    
    /**
     * class which represents a single ray
     */
    static private class Ray {
        static Random rNumbGen = new Random();
        final Vector supportVec, directionVec;
        
        /**
         * creats a ray based on a equation of line
         * @param sVec support vector
         * @param dVec direction vector
         */
        Ray (Vector sVec, Vector dVec) {
            this.supportVec = sVec;
            this.directionVec = dVec;
        }
        
        /**
         * creates a gaussian distributed random ray in the plane clamped by the
         * two plane axis vectors
         * @param direction direction vector of the ray
         * @param planeAxis1 first axis of the spawn plane
         * @param planeAxis2 second axis of the spawn plane
         * @param fwhm gaussian full width half maximum
         * @return created ray
         */
        static Ray createGaussianRandomRay(Vector direction, Vector planeAxis1, Vector planeAxis2, double fwhm) {
            double sigma = fwhm/2.3548;
            double x = rNumbGen.nextGaussian() * sigma;
            double y = rNumbGen.nextGaussian() * sigma;
            Vector support = Vector.add(Vector.times(x, planeAxis1), Vector.times(y, planeAxis2));
            return new Ray(support, direction);
        }
    }
    
    /**
     * calculates based on a equation system the intersection point of a ray
     * with a mirror of the dmd
     * @param r ray
     * @param mx mirror index in x
     * @param my mirror index in y
     * @return double array of intersection coordinates {s, t, u, tilt angle in degree} 
     */
    double[] calcIntersection(Ray r, int mx, int my) {
        double ax = r.directionVec.getX();
        double ay = r.directionVec.getY();
        double az = r.directionVec.getZ();
        double hx = r.supportVec.getX();
        double hy = r.supportVec.getY();
        double hz = r.supportVec.getZ();
        double m = meta.latticeConstant;
        double gammaD = dsc.tiltStates[my][mx] ? tiltAngle : -tiltAngle;
        double gammaR = gammaD /180.0 * Math.PI;
        
        double denominator = (2*az*Math.cos(gammaR)+Math.sqrt(2)*ax*Math.sin(gammaR)-Math.sqrt(2)*ay*Math.sin(gammaR));
        double s = -(((-ay*hz-az*(-hy+m*my))*(-(1/2)*az*(1-Math.cos(gammaR))+(ax*Math.sin(gammaR))/Math.sqrt(2))-(-ax*hz-az*(-hx+m*mx))*(-az*(1/2*(1-Math.cos(gammaR))+Math.cos(gammaR))+(ay*Math.sin(gammaR))/Math.sqrt(2)))
                /(-az*az*Math.cos(gammaR)-(ax*az*Math.sin(gammaR))/Math.sqrt(2)+(ay*az*Math.sin(gammaR))/Math.sqrt(2)));
        double t = -((az*hx-az*hy-ax*hz+ay*hz-az*m*mx+az*m*my-az*hx*Math.cos(gammaR)-az*hy*Math.cos(gammaR)+ax*hz*Math.cos(gammaR)+ay*hz*Math.cos(gammaR)+az*m*mx*Math.cos(gammaR)+az*m*my*Math.cos(gammaR)+Math.sqrt(2)*ay*hx*Math.sin(gammaR)-Math.sqrt(2)*ax*hy*Math.sin(gammaR)-Math.sqrt(2)*ay*m*mx*Math.sin(gammaR)+Math.sqrt(2)*ax*m*my*Math.sin(gammaR))
                /denominator);
        double u = -((2*hz*Math.cos(gammaR)+Math.sqrt(2)*hx*Math.sin(gammaR)-Math.sqrt(2)*hy*Math.sin(gammaR)-Math.sqrt(2)*m*mx*Math.sin(gammaR)+Math.sqrt(2)*m*my*Math.sin(gammaR))
                /denominator);
        
        return new double[]{s,t,u, gammaD};
    }
    
    /**
     * finds intersections with the dmd and checks if the found intersection
     * is physicaly valid
     * @param r ray
     * @return null if no valid intersection was found or the vector of the
     * intersection koordinate with the dmd
     */
    Vector findValidIntersection(Ray r) {
        double smallestU  = Double.MAX_VALUE;
        int finalMx = -1, finalMy = -1;
        double finalS = -1, finalT = -1, finalGamma = 0;
        for (int my = 0; my < nrY; my++) {
            for (int mx = 0; mx < nrY; mx++) {
                double[] intersection = calcIntersection(r, mx, my);
                double s = intersection[0];
                double t = intersection[1];
                double u = intersection[2];
                if (s > mirrorSize || t > mirrorSize || s < 0 || t < 0) continue;
                if (u < smallestU) {
                    smallestU = u;
                    finalMx = mx;
                    finalMy = my;
                    finalS = s;
                    finalT = t;
                    finalGamma = intersection[3];
                }
            }
        }
        if (finalMx >= 0 && finalMy >= 0) {
            Vector validIntersection = dsc.dmd.getCoordinates(finalMx, finalMy, finalGamma, finalS, finalT);
            return validIntersection;
        } else return null;
        
    }

    @Override
    public void simulate() {
        
        Vector inBeam = new Vector(phiInStart / 180. * Math.PI, thetaInStart / 180. * Math.PI);
        inBeam.times(-1);
        Vector[] intersections = new Vector[nrRays];
        Vector[] orthoNormalBasis = Vector.getOrthoNormalBasis(inBeam, new Vector(1,0,0), new Vector(0,1,0));
        Vector planeAxis1 = orthoNormalBasis[1];
        Vector planeAxis2 = orthoNormalBasis[2];
        IntStream intersectionStream = IntStream.range(0, nrRays);
        intersectionStream.forEach(ray -> {
            Ray r = Ray.createGaussianRandomRay(inBeam, planeAxis1, planeAxis2, beamDiameter);
            r.supportVec.setX(r.supportVec.getX() + dsc.dmd.dmdWidth/2);
            r.supportVec.setY(r.supportVec.getY() + dsc.dmd.dmdHeight/2);
            intersections[ray] = findValidIntersection(r);
        });
        
        // init images for simulations
        Image intensity = new Image(width, height);
        intensity.setTitle(String.valueOf(lambda) + "_intensity");
        intensity.show();
        
        Complex[][] field = new Complex[dsc.tMax][dsc.pMax];
        for (int th = 0; th < dsc.tMax; th++) {
            for (int ph = 0; ph < dsc.pMax; ph++) {
                field[th][ph] = new Complex(0, 0);
            }
        }
        
        for (int th = 0; th < dsc.tMax; th++) {
            double theta = dsc.thetaMinR + th * dsc.outStepSizeR;
            final int finalTh = th;
            IntStream pStream = IntStream.range(0, dsc.pMax).parallel();
            pStream.forEach(ph -> {
                double phi = dsc.phiMinR + ph * dsc.outStepSizeR;
                for (int inters = 0; inters < nrRays; inters++) {
                    if (intersections[inters] == null) continue;
                    double argIn = 2*Math.PI/lambdaUm*intersections[inters].dotProduct(inBeam);
                    Vector out = new Vector(phi, theta);

                    double argOut = 2*Math.PI/lambdaUm*intersections[inters].dotProduct(out);
                    double arg = argIn - argOut;
                    double re = Math.cos(arg);
                    double im = Math.sin(arg);

                    field[finalTh][ph].add(new Complex(re, im));
                }
            });
            if (th % 10 == 0) {
                IJ.log(th + "/" + dsc.tMax);
                intensity.set(DmdSimulationCore.buildIntensityImage(field));
                intensity.repaint();
            }
            
            
        }
        intensity.set(DmdSimulationCore.buildIntensityImage(field));
        intensity.repaint();
        intensity.saveAsTiff(outDir + lambda + "_ray_intensity.tif", meta);
        intensity.close();
    }
    
    /**
     * starts the ray tracing approach, values in this method need
     * to be adjusted for the desired system conditions
     * @param args 
     */
    public static void main(String[] args) {
        // creating meta data object
        DmdSimulationCore.MetaData meta = new DmdSimulationCore.MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = false;
        
        meta.lambdas = new int[]{532};

        meta.nrX = 50;
        meta.nrY = 50;

        meta.latticeConstant = 7.56;
        meta.fillFactor = 0.92;
        meta.tiltAngle = 12.0;

        meta.beamDiameter = (int) (Math.min(meta.nrX, meta.nrY) * meta.latticeConstant / 2.0);

        meta.phiOutStart = -1.8;
        meta.phiOutEnd = 8.2;
        meta.thetaOutStart = -8.2;
        meta.thetaOutEnd = 1.8;
        meta.outStepSize = 0.02;

        meta.phiInStart = -21;
        meta.thetaInStart = 21;
        meta.inStepSize = 1.0;
        
        meta.bmp = Image.readBitmap("D:\\dmd-simulator-images\\interesting patterns\\dots-tri-50.bmp");
        //meta.bmp = new Image(meta.nrX, meta.nrY);
        //meta.bmp.setAll(1);
        
        int nrOfRays = 10000;
        RayTracingSimulator rts = new RayTracingSimulator(meta, nrOfRays);
        long timeStart = System.currentTimeMillis();
        rts.simulate();
        System.out.println("Time in seconds: " + ((System.currentTimeMillis() - timeStart) * 0.001));
    }
    
}
