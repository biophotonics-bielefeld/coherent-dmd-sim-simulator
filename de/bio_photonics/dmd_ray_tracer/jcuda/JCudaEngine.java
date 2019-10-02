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
package de.bio_photonics.dmd_ray_tracer.jcuda;

import de.bio_photonics.coherent_dmd_sim_simulator.Complex;
import de.bio_photonics.coherent_dmd_sim_simulator.Vector;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.BufferOverflowException;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import jcuda.Sizeof;
import jcuda.Pointer;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import static jcuda.driver.JCudaDriver.cuCtxCreate;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import static jcuda.driver.JCudaDriver.cuMemAlloc;
import static jcuda.driver.JCudaDriver.cuMemFree;
import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;
import static jcuda.driver.JCudaDriver.cuMemcpyHtoD;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleLoad;
import static jcuda.driver.JCudaDriver.cuModuleGetGlobal;

/**
 * this class acts as an interface between the simulations and JCUDA http://www.jcuda.org/
 * @author m.lachetta
 */
public class JCudaEngine {
    
    private final CUdevice device;
    private final CUcontext context;
    private final Map<String, Module> modules;
    
    /**
     * Initialize the driver and create a context for the first device.
     */
    public JCudaEngine() {
        cuInit(0);
        device = new CUdevice();
        cuDeviceGet(device, 0);
        context = new CUcontext();
        cuCtxCreate(context, 0, device);
        modules = new ConcurrentHashMap<>();
    }
    
    /**
     * class to handle CUmodules
     */
    private class Module {
        
        private final CUmodule cuModule;
        private final Map<String, CUfunction> cuFunctions;
        
        /**
         * 
         * @param moduleName free choosable name for this module
         * @param cuFilePath path to the related cu file
         * @throws IOException 
         */
        private Module(String moduleName, String cuFilePath) throws IOException {
            String ptxFileName = preparePtxFile(cuFilePath);
            cuModule = new CUmodule();
            cuModuleLoad(cuModule, ptxFileName);
            cuFunctions = new ConcurrentHashMap<>();
            
        }
    
        /**
         * copies a float array to the constant memory of the gpu
         * @param constantName name of the __constant__ float array
         * @param hostData data to copy
         */
        private void writeConstant(String constantName, float[] hostData) {
             // Obtain the pointer to the constant memory, and print some info
            CUdeviceptr constantMemoryPointer = new CUdeviceptr();
            long constantMemorySizeArray[] = { 0 };
            cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray, 
                cuModule, constantName);
            int constantMemorySize = (int)constantMemorySizeArray[0];
            int maxElements = constantMemorySize / Sizeof.FLOAT;
            if (hostData.length > maxElements) throw new BufferOverflowException();
            cuMemcpyHtoD(constantMemoryPointer, Pointer.to(hostData), constantMemorySize);
        }
        
        /**
         * copies an int array to the constant memory of the gpu
         * @param constantName name of the __constant__ int array
         * @param hostData data to copy
         */
        private void writeConstant(String constantName, int[] hostData) {
             // Obtain the pointer to the constant memory, and print some info
            CUdeviceptr constantMemoryPointer = new CUdeviceptr();
            long constantMemorySizeArray[] = { 0 };
            cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray, 
                cuModule, constantName);
            int constantMemorySize = (int)constantMemorySizeArray[0];
            int maxElements = constantMemorySize / Sizeof.INT;
            if (hostData.length > maxElements) throw new BufferOverflowException();
            cuMemcpyHtoD(constantMemoryPointer, Pointer.to(hostData), constantMemorySize);
        }
    }
    
    /**
     * loads a CUDA module
     * @param moduleName free choosable name for this module
     * @param cuFilePath path to the related cu file
     * @throws IOException 
     */
    public void loadModule(String moduleName, String cuFilePath) throws IOException {
        Module module = new Module(moduleName, cuFilePath);
        modules.put(moduleName, module);
    }
    
    /**
     * copies a double array to the constant memory of the gpu as floats
     * @param moduleName name of the related CUDA module
     * @param constantName name of the __constant__ float array
     * @param dData data to copy
     */
    public void writeConstant(String moduleName, String constantName, double[] dData) {
        float[] fData = new float[dData.length];
        for (int i = 0; i < dData.length; i++) fData[i] = (float) dData[i];
        writeConstant(moduleName, constantName, fData);
    }
    
    /**
     * copies an int array to the constant memory of the gpu
     * @param moduleName name of the related CUDA module
     * @param constantName name of the __constant__ int array
     * @param data data to copy
     */
    public void writeConstant(String moduleName, String constantName, int[] data) {
        modules.get(moduleName).writeConstant(constantName, data);
    }
    
    /**
     * copies an float array to the constant memory of the gpu
     * @param moduleName name of the related CUDA module
     * @param constantName name of the __constant__ float array
     * @param fData data to copy
     */
    public void writeConstant(String moduleName, String constantName, float[] fData) {
        modules.get(moduleName).writeConstant(constantName, fData);
    }
    
    /**
     * Class to manage arrays
     */
    public abstract class Array {
        protected final CUdeviceptr cuPointer;
        public abstract void hostToDevice();
        public abstract void deviceToHost();
        
        protected Array() {
            cuPointer = new CUdeviceptr();
        }
        
        @Override
        protected void finalize() {
            cuMemFree(cuPointer);
        }
    }
    
    /**
     * Class to manage int arrays on CPU & GPU 
     */
    public class IntArray extends Array{
        
        private int[] array;
        
        private IntArray(int[] array) {
            super();
            this.array = array;
            cuMemAlloc(cuPointer, array.length * Sizeof.INT);
        }
        
        public int[] getArray() {return array;}
        
        /**
         * copies data from the CPU to the GPU
         */
        public void hostToDevice() {
            cuMemcpyHtoD(cuPointer, Pointer.to(array), array.length * Sizeof.INT);
        }
        
        /**
         * copies data from the GPU to the CPU
         */
        public void deviceToHost() {
            cuMemcpyDtoH(Pointer.to(array), cuPointer, array.length * Sizeof.INT);
        }
    }
    
    /**
     * Class to manage float arrays on CPU & GPU 
     */
    public class FloatArray extends Array{
        
        private float[] array;
        
        private FloatArray(float[] array) {
            super();
            this.array = array;
            cuMemAlloc(cuPointer, array.length * Sizeof.FLOAT);
        }
        
        public float[] getArray() {return array;}
        
        /**
         * copies data from the CPU to the GPU
         */
        public void hostToDevice() {
            cuMemcpyHtoD(cuPointer, Pointer.to(array), array.length * Sizeof.FLOAT);
        }
        
        /**
         * copies data from the GPU to the CPU
         */
        public void deviceToHost() {
            cuMemcpyDtoH(Pointer.to(array), cuPointer, array.length * Sizeof.FLOAT);
        }
        
        /**
         * converts this 2*w*h array into a 2D complex array
         * @param width
         * @param height
         * @return 2D complex array
         */
        public Complex[][] toComplex2d(int width, int height) {
            Complex[][] cArray = new Complex[height][width];
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    float re = array[y * 2 * width + 2 * x];
                    float im = array[y * 2 * width + 2 * x + 1];
                    Complex c = new Complex(re, im);
                    cArray[y][x] = c;
                }
            }
            return cArray;
        }
    }
    
    /**
     * Allocates memory for a float array on the GPU connected to the CPU array
     * @param array native data
     * @return gpu/cpu array
     */
    public FloatArray createFloatArray (float[] array) {
        return new FloatArray(array);
    }
    
    /**
     * Allocates memory for a float array on the GPU connected to the CPU array
     * @param dArray native data
     * @return gpu/cpu array
     */
    public FloatArray createFloatArray(double[][] dArray) {
        int width = dArray[0].length;
        int height = dArray.length;
        float[] array = new float[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                array[y * width + x] = (float) dArray[y][x];
            }
        }
        return new FloatArray(array);
    }

    /**
     * Allocates memory for a float array on the GPU connected to the CPU array
     * @param cArray native data
     * @return gpu/cpu array
     */
    public FloatArray createFloatArray(Complex[][] cArray) {
        int width = cArray[0].length;
        int height = cArray.length;
        float[] array = new float[width * height * 2];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                array[y * 2 * width + 2 * x] = (float) cArray[y][x].getRe();
                array[y * 2 * width + 2 * x + 1] = (float) cArray[y][x].getIm();
            }
        }
        return new FloatArray(array);
    }
    
    /**
     * Allocates memory for a float array on the GPU connected to the CPU array
     * @param v native data
     * @return gpu/cpu array
     */
    public FloatArray createFloatArray(Vector v) {
        float[] array = new float[3];
        array[0] = (float) v.getX();
        array[1] = (float) v.getY();
        array[2] = (float) v.getZ();
        return new FloatArray(array);
    }
    
    /**
     * Allocates memory for a float array on the GPU connected to the CPU array
     * @param vArray native data
     * @return gpu/cpu array
     */
    public FloatArray createFloatArray(Vector[][] vArray) {
        int width = vArray[0].length;
        int height = vArray.length;
        float[] array = new float[width * height * 3];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                array[y * 3 * width + 3 * x] = (float) vArray[y][x].getX();
                array[y * 3 * width + 3 * x + 1] = (float) vArray[y][x].getY();
                array[y * 3 * width + 3 * x + 2] = (float) vArray[y][x].getZ();
            }
        }
        return new FloatArray(array);
    }
    
    /**
     * Allocates memory for an int array on the GPU connected to the CPU array
     * @param array native data
     * @return gpu/cpu array
     */
    public IntArray createIntArray(int[] array) {
        return new IntArray(array);
    }
    
    /**
     * Allocates memory for an int array on the GPU connected to the CPU array
     * @param bArray native data
     * @return gpu/cpu array
     */
    public IntArray createIntArray(boolean[][] bArray) {
        int width = bArray[0].length;
        int height = bArray.length;
        int[] array = new int[width * height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                array[y * width + x] = bArray[y][x] ? 1 : 0;
            }
        }
        return new IntArray(array);
    }
    
    /**
     * original from http://www.jcuda.org/
     * Starts a CUDA kernel
     * @param moduleName name of the related CUDA module
     * @param functionName name of the kernel function
     * @param blockSizeX CUDA blockSizeX
     * @param numElements
     * @param mStart
     * @param mEnd
     * @param arrays
     */
    public void launchKernel(String moduleName, String functionName, int blockSizeX,
            int numElements, int mStart, int mEnd, Array... arrays) {
        // Set up the kernel parameters: A pointer to an array
        // of pointers which point to the actual values.
        Pointer[] pointers = new Pointer[arrays.length + 3];
        pointers[0] = Pointer.to(new int[]{numElements});
        pointers[1] = Pointer.to(new int[]{mStart});
        pointers[2] = Pointer.to(new int[]{mEnd});
        for(int i = 0; i < arrays.length; i ++) {
            CUdeviceptr cuPointer = arrays[i].cuPointer;
            pointers[i+3] = Pointer.to(cuPointer);
        }
        Pointer kernelParameters = Pointer.to(pointers);

        // Call the kernel function.
        int gridSizeX = (int)Math.ceil((double)numElements / blockSizeX);
        CUfunction function = modules.get(moduleName).cuFunctions.get(functionName);
        cuLaunchKernel(function,
            gridSizeX,  1, 1,      // Grid dimension
            blockSizeX, 1, 1,      // Block dimension
            0, null,               // Shared memory size and stream
            kernelParameters, null // Kernel- and extra parameters
        );
        cuCtxSynchronize();
    }
    
    public void launchKernel(String moduleName, String functionName, int blockSizeX,
            int numElements, double latticeConstant, double inBeamX, double inBeamY, Array... arrays) {
        // Set up the kernel parameters: A pointer to an array
        // of pointers which point to the actual values.
        Pointer[] pointers = new Pointer[arrays.length + 4];
        pointers[0] = Pointer.to(new int[]{numElements});
        pointers[1] = Pointer.to(new double[]{latticeConstant});
        pointers[2] = Pointer.to(new double[]{inBeamX});
        pointers[3] = Pointer.to(new double[]{inBeamY});
        for(int i = 0; i < arrays.length; i ++) {
            CUdeviceptr cuPointer = arrays[i].cuPointer;
            pointers[i+4] = Pointer.to(cuPointer);
        }
        Pointer kernelParameters = Pointer.to(pointers);

        // Call the kernel function.
        int gridSizeX = (int)Math.ceil((double)numElements / blockSizeX);
        CUfunction function = modules.get(moduleName).cuFunctions.get(functionName);
        cuLaunchKernel(function,
            gridSizeX,  1, 1,      // Grid dimension
            blockSizeX, 1, 1,      // Block dimension
            0, null,               // Shared memory size and stream
            kernelParameters, null // Kernel- and extra parameters
        );
        cuCtxSynchronize();
    }
    
    public void launchKernel(String moduleName, String functionName, int blockSizeX,
            int numElements, double latticeConstant, double inBeamX, double inBeamY,
            double inBeamZ, double alpha, Array... arrays) {
        // Set up the kernel parameters: A pointer to an array
        // of pointers which point to the actual values.
        Pointer[] pointers = new Pointer[arrays.length + 6];
        pointers[0] = Pointer.to(new int[]{numElements});
        pointers[1] = Pointer.to(new double[]{latticeConstant});
        pointers[2] = Pointer.to(new double[]{inBeamX});
        pointers[3] = Pointer.to(new double[]{inBeamY});
        pointers[4] = Pointer.to(new double[]{inBeamZ});
        pointers[5] = Pointer.to(new double[]{alpha});
        for(int i = 0; i < arrays.length; i ++) {
            CUdeviceptr cuPointer = arrays[i].cuPointer;
            pointers[i+6] = Pointer.to(cuPointer);
        }
        Pointer kernelParameters = Pointer.to(pointers);

        // Call the kernel function.
        int gridSizeX = (int)Math.ceil((double)numElements / blockSizeX);
        CUfunction function = modules.get(moduleName).cuFunctions.get(functionName);
        cuLaunchKernel(function,
            gridSizeX,  1, 1,      // Grid dimension
            blockSizeX, 1, 1,      // Block dimension
            0, null,               // Shared memory size and stream
            kernelParameters, null // Kernel- and extra parameters
        );
        cuCtxSynchronize();
    }
    
    /**
     * copied from http://www.jcuda.org/
     * Enables or disables exceptions. By default, the methods of this class
     * only return the CUresult error code from the underlying CUDA function.
     * If exceptions are enabled, a CudaException with a detailed error message
     * will be thrown if a method is about to return a result code that is not
     * CUresult.CUDA_SUCCESS
     * @param b Whether exceptions are enabled
     */
    public static void setExceptionsEnabled(boolean b) {
        JCudaDriver.setExceptionsEnabled(b);
    }
    
    /**
     * loads a function from a module
     * @param moduleName name of the related CUDA module
     * @param functionName name of the kernel function
     */
    public void loadFunktion(String moduleName, String functionName) {
        Module module = modules.get(moduleName);
        CUfunction function = new CUfunction();
        cuModuleGetFunction(function, module.cuModule, functionName);
        module.cuFunctions.put(functionName, function);
    }
    
    /**
     * copied from http://www.jcuda.org/
     * The extension of the given file name is replaced with "ptx".
     * If the file with the resulting name does not exist, it is
     * compiled from the given file using NVCC. The name of the
     * PTX file is returned.
     *
     * @param cuFileName The name of the .CU file
     * @return The name of the PTX file
     * @throws IOException If an I/O error occurs
     */
    private static String preparePtxFile(String cuFileName) throws IOException
    {
        int endIndex = cuFileName.lastIndexOf('.');
        if (endIndex == -1)
        {
            endIndex = cuFileName.length()-1;
        }
        String ptxFileName = cuFileName.substring(0, endIndex+1)+"ptx";
        File ptxFile = new File(ptxFileName);
        if (ptxFile.exists())
        {
            return ptxFileName;
        }

        File cuFile = new File(cuFileName);
        if (!cuFile.exists())
        {
            throw new IOException("Input file not found: "+cuFileName);
        }
        String modelString = "-m"+System.getProperty("sun.arch.data.model");
        String command =
            "nvcc -lm " + modelString + " -ptx "+
            cuFile.getPath()+" -o "+ptxFileName;

        System.out.println("Executing\n"+command);
        Process process = Runtime.getRuntime().exec(command);

        String errorMessage =
            new String(toByteArray(process.getErrorStream()));
        String outputMessage =
            new String(toByteArray(process.getInputStream()));
        int exitValue = 0;
        try
        {
            exitValue = process.waitFor();
        }
        catch (InterruptedException e)
        {
            Thread.currentThread().interrupt();
            throw new IOException(
                "Interrupted while waiting for nvcc output", e);
        }

        if (exitValue != 0)
        {
            System.out.println("nvcc process exitValue "+exitValue);
            System.out.println("errorMessage:\n"+errorMessage);
            System.out.println("outputMessage:\n"+outputMessage);
            throw new IOException(
                "Could not create .ptx file: "+errorMessage);
        }

        System.out.println("Finished creating PTX file");
        return ptxFileName;
    }
    
    /**
     * copied from http://www.jcuda.org/
     * Fully reads the given InputStream and returns it as a byte array
     *
     * @param inputStream The input stream to read
     * @return The byte array containing the data from the input stream
     * @throws IOException If an I/O error occurs
     */
    private static byte[] toByteArray(InputStream inputStream)
        throws IOException
    {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        byte buffer[] = new byte[8192];
        while (true)
        {
            int read = inputStream.read(buffer);
            if (read == -1)
            {
                break;
            }
            baos.write(buffer, 0, read);
        }
        return baos.toByteArray();
    }
}
