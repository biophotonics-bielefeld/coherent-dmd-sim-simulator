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
 * Simple 3 dimensional vector class containing doubles
 * @author m.lachetta
 */
public class Vector {
    
    private double x, y, z;
    
    /**
     * Creates a new Vector
     * @param x
     * @param y
     * @param z 
     */
    public Vector(double x, double y, double z) {
        set(x,y,z);
    }
    
    public Vector createCopy() {
        return new Vector(x, y, z);
    }
    
    /**
     * Creates a normalized Vector dependent on the incidence/reflection angle
     * x = z * tan(phi);
     * y = z * tan(theta);
     * z = sqrt(tan(phi)^2 + tan(theta)^2 + 1);
     * @param phi x angle in rad
     * @param theta y angle in rad
     */
    public Vector(double phi, double theta) {
        double tp = Math.tan(phi);
        double tt = Math.tan(theta);
        z = Math.sqrt(1 / (tp*tp + tt*tt + 1));
        x = z*tp;
        y = z*tt;
    }
    
    @Override
    public String toString() {
        return x + "\t" + y + "\t" + z;
    }
    
    public void set(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    
    public void setX(double x) {
        this.x = x;
    }
    
    public void setY(double y) {
        this.y = y;
    }
    
    public void setZ(double z) {
        this.z = z;
    }
    
    public double getX() {
        return x;
    }
    
    public double getY() {
        return y;
    }
    
    public double getZ() {
        return z;
    }
    
    /**
     * 
     * @return this vector as [3]array
     */
    public double[] getArray() {
        double[] array = {x,y,z};
        return array;
    }
    
    /**
     * add a vector to this
     * @param v the other vector
     */
    public void add(Vector v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }
    
    public static Vector add(Vector v1, Vector v2) {
        return new Vector(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
    }
    
    public static Vector minus(Vector v1, Vector v2) {
        return new Vector(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
    }
    
    /**
     * calculates the dot product
     * @param v the other vector
     * @return dot product
     */
    public double dotProduct(Vector v) {
        return x * v.x + y * v.y + z * v.z;
    }
    
    /**
     * multiplies this vector with a scalar
     * @param a scalar
     */
    public void times(double a) {
        x *= a;
        y *= a;
        z *= a;
    }
    
    public static Vector times(double a, Vector v) {
        return new Vector(a*v.x, a*v.y, a*v.z);
    }
    
    /**
     * divides this vector by a scalar
     * @param a 
     */
    public void divide (double a) {
        if (a == 0) throw new ArithmeticException("Division by zero");
        times((1.0 / a));
    }
    
    public double getAbs() {
        return x*x+y*y+z*z;
    }
    
    /**
     * calculates length with the 2-norm
     * @return length of this vector
     */
    public double getNorm() {
        return Math.sqrt(getAbs());
    }
    
    /**
     * sets the length/norm of this vector to 1
     */
    public void normalize() {
        divide(getNorm());
    }
    
    public Vector projectOnPlane(Vector planeNormal) {
        double t = - (this.dotProduct(planeNormal)) / planeNormal.getAbs();
        return add(this, times(t, planeNormal));
    }
    
    
    
    public static Vector[] getOrthoNormalBasis(Vector w1, Vector w2, Vector w3) {
        Vector v1 = w1.createCopy();
        v1.normalize();
        Vector v2 = minus(w2, times(v1.dotProduct(w2), v1));
        v2.normalize();
        Vector v3 = minus(minus(w3, times(v1.dotProduct(w3), v1)), times(v2.dotProduct(w3), v2));
        v3.normalize();
        return new Vector[]{v1, v2, v3};
    }
    
    /**
     * main method for testing
     * @param args 
     */
    public static void main(String[] args) {
        Vector v = new Vector(Math.PI*0.25, Math.PI*0.25*2);
        System.out.println(v);
        Vector[] basis = getOrthoNormalBasis(v, new Vector(1,0,0), new Vector(0,1,0));
        Vector v1 = basis[0];
        Vector v2 = basis[1];
        Vector v3 = basis[2];
        System.out.println(v1);
        System.out.println(v2);
        System.out.println(v3);
        System.out.println("");
        System.out.println(v1.dotProduct(v2));
        System.out.println(v1.dotProduct(v3));
        System.out.println(v2.dotProduct(v3));
    }
    
}
