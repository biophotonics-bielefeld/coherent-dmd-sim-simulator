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
        return "{" + x + "," + y + "," + z + "}";
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
    
    /**
     * calculates the dot product
     * @param v the other vector
     * @return dot product
     */
    public double times(Vector v) {
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
    
    /**
     * divides this vector by a scalar
     * @param a 
     */
    public void divide (double a) {
        if (a == 0) throw new ArithmeticException("Division by zero");
        times((1.0 / a));
    }
    
    /**
     * calculates length with the 2-norm
     * @return length of this vector
     */
    public double getNorm() {
        return Math.sqrt(x*x+y*y+z*z);
    }
    
    /**
     * sets the length/norm of this vector to 1
     */
    public void normalize() {
        divide(getNorm());
    }
    
    /**
     * main method for testing
     * @param args 
     */
    public static void main(String[] args) {
        Vector v = new Vector(Math.PI*0.25, Math.PI*0.25);
        System.out.println(v);
    }
    
}
