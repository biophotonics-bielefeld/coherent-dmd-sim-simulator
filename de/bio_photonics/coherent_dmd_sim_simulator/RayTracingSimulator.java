/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.bio_photonics.coherent_dmd_sim_simulator;

import ij.IJ;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.stream.IntStream;

/**
 *
 * @author Mario
 */
public class RayTracingSimulator extends AbstractSimulator {
    
    
    final int nrRays;
    
    public RayTracingSimulator(DmdSimulationCore.MetaData meta, int nrRays) {
        super(meta);
        this.nrRays = nrRays;
    }
    
    static private class Ray {
        static Random rNumbGen = new Random();
        final Vector supportVec, directionVec;
        
        Ray (Vector sVec, Vector dVec) {
            this.supportVec = sVec;
            this.directionVec = dVec;
        }
        
        static Ray createGaussianRandomRay(Vector direction, double fwhm) {
            double sigma = fwhm/2.3548;
            double x = rNumbGen.nextGaussian() * sigma;
            double y = rNumbGen.nextGaussian() * sigma;
            //System.out.println(x + "\t" + y);
            Vector support = new Vector(x, y, 0);
            support = support.projectOnPlane(direction);
            return new Ray(support, direction);
        }
        
        Vector get(double param) {
            double x = supportVec.getX() + directionVec.getX() * param;
            double y = supportVec.getY() + directionVec.getY() * param;
            double z = supportVec.getZ() + directionVec.getZ() * param;
            return new Vector(x, y, z);
        }
    }
    
    /**
     * 
     * @param r
     * @param mx
     * @param my
     * @return double array of s t u intersection coordinates
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
                //System.out.println(s + " " + t + " " + u);
                if (s > mirrorSize || t > mirrorSize || s < 0 || t < 0) continue;
                if (u < smallestU) {
                    smallestU = u;
                    finalMx = mx;
                    finalMy = my;
                    finalS = s;
                    finalT = t;
                    finalGamma = intersection[3];
                    //System.out.println(smallestU + " " + finalMx + " " + finalMy + " " + finalS + " " + finalT);
                }
            }
        }
        if (finalMx >= 0 && finalMy >= 0) {
            Vector validIntersection = dsc.dmd.getCoordinates(finalMx, finalMy, finalGamma, finalS, finalT);
            //System.out.println(validIntersection.getX() + "\t" + validIntersection.getY() + "\t" + validIntersection.getZ());
            return validIntersection;
        } else return null;
        
    }

    @Override
    public void simulate() {
        
        Vector inBeam = new Vector(phiInStart / 180. * Math.PI, thetaInStart / 180. * Math.PI);
        inBeam.times(-1);
        //Ray[] rays = new Ray[nrRays];
        Vector[] intersections = new Vector[nrRays];
        
        IntStream intersectionStream = IntStream.range(0, nrRays);
        intersectionStream.forEach(ray -> {
            Ray r = Ray.createGaussianRandomRay(inBeam, beamDiameter);
            r.supportVec.setX(r.supportVec.getX() + dsc.dmd.dmdWidth/2);
            r.supportVec.setY(r.supportVec.getY() + dsc.dmd.dmdHeight/2);
            intersections[ray] = findValidIntersection(r);
        });
        try {
            FileWriter fw = new FileWriter("D:\\dmd-simulator-images\\diff_points.txt");
            for(Vector inters : intersections) {
                if (inters == null) continue;
                //System.out.println(inters.getX() + "\t" + inters.getY() + "\t" + inters.getZ());
                fw.write(inters.getX() + "\t" + inters.getY() + "\t" + inters.getZ() + "\n");
            }
            fw.close();
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        /*
        for (int ray = 0; ray < nrRays; ray++) {
            Ray r = Ray.createGaussianRandomRay(inBeam, beamDiameter);
            r.supportVec.setX(r.supportVec.getX() + dsc.dmd.dmdWidth/2);
            r.supportVec.setY(r.supportVec.getY() + dsc.dmd.dmdHeight/2);
            intersections[ray] = findValidIntersection(r);
        }
        */
        // init images for simulations
        Image intensity = new Image(width, height);
        //Image phase = new Image(width, height);
        intensity.setTitle(String.valueOf(lambda) + "_intensity");
        //phase.setTitle(String.valueOf(lambda) + "_phase");
        intensity.show();
        //phase.show();
        
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
                    double argIn = 2*Math.PI/lambdaUm*intersections[inters].times(inBeam);
                    Vector out = new Vector(phi, theta);

                    double argOut = 2*Math.PI/lambdaUm*intersections[inters].times(out);
                    double arg = argIn - argOut;
                    //System.out.println(inBeam + "\t" + out);
                    double re = Math.cos(arg);
                    double im = Math.sin(arg);

                    field[finalTh][ph].add(new Complex(re, im));
                }
            });
            if (th % 10 == 0) {
                IJ.log(th + "/" + dsc.tMax);
                intensity.set(DmdSimulationCore.buildIntensityImage(field));
                intensity.repaint();
                //phase.set(DmdSimulationCore.buildPhaseImage(field));
                //phase.repaint();
            }
            
            
        }
        intensity.set(DmdSimulationCore.buildIntensityImage(field));
        intensity.repaint();
        intensity.saveAsTiff(outDir + lambda + "_ray_intensity.tif", meta);
        //phase.saveAsTiff(outDir + lambda + "_ray_phase.tif", meta);
        
        // cloeses all images
        intensity.close();
        //phase.close();
    }
    
    public static void main(String[] args) {
        // creating meta data object
        DmdSimulationCore.MetaData meta = new DmdSimulationCore.MetaData();
        meta.outDir = "D:\\dmd-simulator-images\\";
        meta.gpuActive = false;
        
        int lambdaStart = 532;
        int lambdaEnd = 700;
        int lambdaStepSize = 100;
        int nrLambdas = (lambdaEnd - lambdaStart) / lambdaStepSize + 1;
        meta.lambdas = new int[(lambdaEnd - lambdaStart) / lambdaStepSize + 1];
        for (int i = 0; i < nrLambdas; i++) {
            meta.lambdas[i] = lambdaStart + i*lambdaStepSize;
            System.out.println(i + " " + meta.lambdas[i]);
        }

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
        meta.phiInEnd = -20;
        meta.thetaInStart = 21;
        meta.thetaInEnd = 22;
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
