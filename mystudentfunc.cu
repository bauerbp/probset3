/* Udacity Homework 3
   HDR Tone-mapping

  Background HDR
  ==============

  A High Dynamic Range (HDR) image contains a wider variation of intensity
  and color than is allowed by the RGB format with 1 byte per channel that we
  have used in the previous assignment.  

  To store this extra information we use single precision floating point for
  each channel.  This allows for an extremely wide range of intensity values.

  In the image for this assignment, the inside of church with light coming in
  through stained glass windows, the raw input floating point values for the
  channels range from 0 to 275.  But the mean is .41 and 98% of the values are
  less than 3!  This means that certain areas (the windows) are extremely bright
  compared to everywhere else.  If we linearly map this [0-275] range into the
  [0-255] range that we have been using then most values will be mapped to zero!
  The only thing we will be able to see are the very brightest areas - the
  windows - everything else will appear pitch black.

  The problem is that although we have cameras capable of recording the wide
  range of intensity that exists in the real world our monitors are not capable
  of displaying them.  Our eyes are also quite capable of observing a much wider
  range of intensities than our image formats / monitors are capable of
  displaying.

  Tone-mapping is a process that transforms the intensities in the image so that
  the brightest values aren't nearly so far away from the mean.  That way when
  we transform the values into [0-255] we can actually see the entire image.
  There are many ways to perform this process and it is as much an art as a
  science - there is no single "right" answer.  In this homework we will
  implement one possible technique.

  Background Chrominance-Luminance
  ================================

  The RGB space that we have been using to represent images can be thought of as
  one possible set of axes spanning a three dimensional space of color.  We
  sometimes choose other axes to represent this space because they make certain
  operations more convenient.

  Another possible way of representing a color image is to separate the color
  information (chromaticity) from the brightness information.  There are
  multiple different methods for doing this - a common one during the analog
  television days was known as Chrominance-Luminance or YUV.

  We choose to represent the image in this way so that we can remap only the
  intensity channel and then recombine the new intensity values with the color
  information to form the final image.

  Old TV signals used to be transmitted in this way so that black & white
  televisions could display the luminance channel while color televisions would
  display all three of the channels.
  

  Tone-mapping
  ============

  In this assignment we are going to transform the luminance channel (actually
  the log of the luminance, but this is unimportant for the parts of the
  algorithm that you will be implementing) by compressing its range to [0, 1].
  To do this we need the cumulative distribution of the luminance values.

  Example
  -------

  input : [2 4 3 3 1 7 4 5 7 0 9 4 3 2]
  min / max / range: 0 / 9 / 9

  histo with 3 bins: [4 7 3]

  cdf : [4 11 14]


  Your task is to calculate this cumulative distribution by following these
  steps.

*/




#include "reference_calc.cpp"
#include "utils.h"
#include <stdio.h>

/*__global__ void minMax(float* sortingArray,
                  float* min_logLum,
                  float* max_logLum,
                  int rows, int cols)
{
    float min = sortingArray[0];
    float max = sortingArray[0];
    int theArraySizehalf = rows*cols/2;
    int2 index2d= make_int2( ( blockIdx.x * blockDim.x ) + threadIdx.x, ( blockIdx.y * blockDim.y ) + threadIdx.y );
    int index1d = index2d.x + index2d.y * cols;
    int counter=1;
    for(int i=1; i < theArraySizehalf; i=i*2){
        if(index%(2^counter)==0 && index < (rows*cols)){
            if(sortingArray[index1d]>sortingArray[index1d+i]){
               
                
            
    
    
    
    return;
}*/

__global__ void global_min_kernel(float* sortingArray)
{
    int myId = threadIdx.x + blockDim.x * blockIdx.x;
    int tid  = threadIdx.x;

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            if(sortingArray[myId] > sortingArray[myId + s])
            {
                sortingArray[myId] = sortingArray[myId + s];
            }
        }
        __syncthreads();        // make sure all adds at one stage are done!
    }

    // only thread 0 writes result for this block back to global mem
    if (tid == 0)
    {
        sortingArray[blockIdx.x] = sortingArray[myId];
    }
}

__global__ void global_max_kernel(float* sortingArray)
{
    int myId = threadIdx.x + blockDim.x * blockIdx.x;
    int tid  = threadIdx.x;

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            if(sortingArray[myId] < sortingArray[myId + s])
            {
                sortingArray[myId] = sortingArray[myId + s];
            }
        }
        __syncthreads();        // make sure all adds at one stage are done!
    }

    // only thread 0 writes result for this block back to global mem
    if (tid == 0)
    {
        sortingArray[blockIdx.x] = sortingArray[myId];
    }
}


__global__ void simple_histo(int *d_bins, const float *d_in, int* runningArray, const int numBins, int lumRange, int lumMin, int totalSize)
{
    
    int myId = threadIdx.x + blockDim.x * blockIdx.x;
    int myItem = d_in[myId];
    int myBin = (myItem - lumMin) / lumRange * numBins;
    runningArray[myId] = myBin;
    /*
    if(myId==0){
        for(int i=0; i<totalSize; i++){
            d_bins[runningArray[i]] ++;
        }
    }*/
    atomicAdd(&(d_bins[myBin]), 1);
}

void your_histogram_and_prefixsum(const float* const d_logLuminance,
                                  unsigned int* const d_cdf,
                                  float &min_logLum,
                                  float &max_logLum,
                                  const size_t numRows,
                                  const size_t numCols,
                                  const size_t numBins)
{
  //TODO
  /*Here are the steps you need to implement
    1) find the minimum and maximum value in the input logLuminance channel
       store in min_logLum and max_logLum
    2) subtract them to find the range
    3) generate a histogram of all the values in the logLuminance channel using
       the formula: bin = (lum[i] - lumMin) / lumRange * numBins
    4) Perform an exclusive scan (prefix sum) on the histogram to get
       the cumulative distribution of luminance values (this should go in the
       incoming d_cdf pointer which already has been allocated for you)       */
       
    printf("woot1");
       
    float* sortingArray;
    checkCudaErrors(cudaMalloc(&sortingArray, numRows*numCols));
    checkCudaErrors(cudaMemcpy(sortingArray,  d_logLuminance,   numRows*numCols, cudaMemcpyHostToDevice));
    
    const dim3 blockSize(32, 32);
    const dim3 gridSize(numCols*numRows/1024);
    
    global_min_kernel<<<gridSize, blockSize>>>(sortingArray);
    min_logLum = sortingArray[0];
    global_max_kernel<<<gridSize, blockSize>>>(sortingArray);
    max_logLum = sortingArray[0];
    int range = max_logLum-min_logLum;
    printf("woot2");
    int* runningArray;
    
    cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
    
    checkCudaErrors(cudaMalloc(&runningArray, numRows*numCols));
    int totalSize=numCols*numRows;
    for(int i=0; i<totalSize; i++){
        runningArray[i]=0;
    }
    checkCudaErrors(cudaMemcpy(runningArray,  runningArray,   numRows*numCols, cudaMemcpyHostToDevice));
    int* d_bins;
    checkCudaErrors(cudaMalloc(&d_bins, numBins));
    for(int i=0; i<numBins; i++){
        d_bins[i]=0;
    }
    
    printf("woot3");
    checkCudaErrors(cudaMemcpy(d_bins,  d_bins,   numBins, cudaMemcpyHostToDevice));
    simple_histo<<<gridSize, blockSize>>>(d_bins, d_logLuminance, runningArray, numBins, range, min_logLum, totalSize);
    
    cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
    
    d_cdf[0]=0;
    for(int i=1; i<numBins; i++){
        d_cdf[i]=d_cdf[i-1]+d_bins[i-1];
    }
    printf("woot4");
}