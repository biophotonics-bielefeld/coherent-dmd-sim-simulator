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
 * Simple class for complex numbers
 * @author m.lachetta
 */
public class Complex {
    double re, im;
    
    /**
     * 
     * @param real real part
     * @param imag imaginary part
     */
    public Complex(double real, double imag) {
        re = real;
        im = imag;
    }
    
    @Override
    public String toString() {
        return "{" + re + "," + im + "}";
    }
    
    /**
     * adds a complex to this complex
     * @param other 
     */
    public void add(Complex other) {
        re += other.re;
        im += other.im;
    }
    
    /**
     * multiplies this complex with the other 
     * @param other 
     */
    public void multi(Complex other) {
        double a = re;
        double b = im;
        double c = other.re;
        double d = other.im;
        re = a*c-b*d;
        im = a*d+b*c;
    }
    
    /**
     * divides this complex by the other
     * @param other 
     */
    public void divide(Complex other) {
        double a = re;
        double b = im;
        double c = other.re;
        double d = other.im;
        re = (a*c+b*d)/(c*c+d*d);
        im = (b*c-a*d)/(c*c+d*d);
    }
    
    /**
     * 
     * @return the absolute value of this complex
     */
    public double abs() {
        return Math.sqrt(re*re + im*im);
    }
    
    /**
     * 
     * @return square of the absolute value
     */
    public double square() {
        return Math.pow(abs(), 2.);
    }
    
    /**
     * 
     * @return argument of this complex number
     */
    public double arg() {
        double r = abs();
        double arg = Math.acos(re / r);
        if (im >= 0) return arg;
        else return -arg;
    }
    
    /**
     * sets the real part of this complex
     * @param re 
     */
    public void setRe(double re) {
        this.re = re;
    }
    
    /**
     * sets the imaginary part of this complex
     * @param im 
     */
    public void setIm(double im) {
        this.im = im;
    }
    
    /**
     * 
     * @return the real part of this complex
     */
    public double getRe() {
        return re;
    }
    
    /**
     * 
     * @return the real part of this complex
     */
    public double getIm() {
        return im;
    }
}
