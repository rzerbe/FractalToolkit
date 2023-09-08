#include "Cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#if _DEBUG
#   define cudaCall(cuda_func, ...) { cudaError_t status = cuda_func(__VA_ARGS__); cudaAssert((status), __FILE__, #cuda_func, __LINE__); }
#else
#   define cudaCall(cuda_func, ...) { cudaError_t status = cuda_func(__VA_ARGS__); }
#endif

inline void cudaAssert(cudaError_t status, const char* file, const char* func, int line)
{
    if (status != cudaSuccess)
    {
        std::stringstream ss;
        ss << "Error: " << cudaGetErrorString(status) << std::endl;
        ss << "Func: " << func << std::endl;
        ss << "File: " << file << std::endl;
        ss << "Line: " << line << std::endl;

        throw std::runtime_error(ss.str());
    }
}

__device__ double* escape_count_ssaa;
__device__ double* magnitude_ssaa;
__device__ double* escape_count;

double* escape_count_ssaa_d;
double* magnitude_ssaa_d;
double* escape_count_d;

__global__ void downsample(
    double* escape_count_ssaa,
    double* magnitude_ssaa,
    double* escape_count,
    const std::uint32_t max_iterations,
    const int width,
    const int height,
    const int ssaa,
    const int smooth_shading)
{
    const int y_pixel = threadIdx.y + blockIdx.y * blockDim.y;
    const int x_pixel = threadIdx.x + blockIdx.x * blockDim.x;

    int acccessor = y_pixel * width + x_pixel;

    // Pixel
    int index;
    int smallestIndex;
    double min;
    double norm;

    min = max_iterations;
    smallestIndex = ((y_pixel * ssaa + 0) * width * ssaa) + (x_pixel * ssaa + 0);
    // Sub Pixel
    for (int y = 0; y < ssaa; ++y)
    {
        for (int x = 0; x < ssaa; ++x)
        {
            index = ((y_pixel * ssaa + y) * width * ssaa) + (x_pixel * ssaa + x);
            if (escape_count_ssaa[index] < min)
            {
                smallestIndex = index;
                min = escape_count_ssaa[index];
            }
        }
    }

    if (smooth_shading == 1)
    {
        norm = min + 1 - log2f(log2f(sqrt(magnitude_ssaa[smallestIndex]))) / log2f(2);
        if (norm > max_iterations || norm < -1000 || min == max_iterations - 1)
        {
            norm = max_iterations;
        }

        escape_count[acccessor] = norm;
    }
    else
    {
        escape_count[acccessor] = min;
    }
}

__global__ void mandelbrot_kernel(
    double* escape_count,
    double* magnitude,
    const std::uint32_t max_iterations,
    const int width, 
    const int height, 
    const double x_scale,
    const double y_scale,
    const double x_min, 
    const double y_min)
{
    const int y_pixel = threadIdx.y + blockIdx.y * blockDim.y;
    const int x_pixel = threadIdx.x + blockIdx.x * blockDim.x;

    int acccessor = y_pixel * width + x_pixel;

    if (y_pixel >= height || x_pixel >= width)
    {
        return;
    }

    double y = y_min + (y_pixel * x_scale);
    double x = x_min + (x_pixel * y_scale);

    double yy = y * y;
    double zx = (x - 0.25) * (x - 0.25) + yy;

    if (zx*(zx + x - 0.25) - 0.25 * yy < 0 || (x + 1.0) * (x + 1.0) + yy - 0.0625 < 0)
    {
        escape_count[acccessor] = max_iterations;
        magnitude[acccessor] = INFINITY;
        return;
    }

    std::uint32_t curr_iterations = 0;
    double zy, zx2, zy2, mag;
    zx = zy = zx2 = zy2 = mag = 0.0;

    do {
        zy = 2.0 * zx * zy + y;
        zx = zx2 - zy2 + x;
        zx2 = zx * zx;
        zy2 = zy * zy;
        mag = zx2 + zy2;
    } while (curr_iterations++ < max_iterations && mag < 4.0);

    escape_count[acccessor] = curr_iterations;
    magnitude[acccessor] = mag;
}

namespace mandelbrot_cuda
{
    bool InitCUDA(Properties* properties)
    {
        const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
        const int width = properties->Get_General_Width();
        const int height = properties->Get_General_Height();

        const int width_ssaa = width * ssaa;
        const int height_ssaa = height * ssaa;

        int count = 0;
        int i = 0;

        cudaCall(cudaGetDeviceCount, &count);
        if (count == 0) {
            fprintf(stderr, "There is no device.\n");
            return false;
        }

        for (i = 0; i < count; i++) {
            cudaDeviceProp prop;
            if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
                if (prop.major >= 1) {
                    break;
                }
            }
        }
        if (i == count) {
            fprintf(stderr, "There is no device supporting CUDA.\n");
            return false;
        }
        cudaCall(cudaSetDevice, i);
        cudaCall(cudaFree, 0);

        const size_t buffer_size_ssaa = width_ssaa * height_ssaa * sizeof(double);
        const size_t buffer_size = width * height * sizeof(double);

        // Allocate memory
        cudaMalloc((void**)&escape_count_ssaa_d, buffer_size_ssaa);
        cudaMalloc((void**)&magnitude_ssaa_d, buffer_size_ssaa);
        cudaMallocHost((void**)&escape_count_d, buffer_size);

        // Zero memory and assign to global device symbol
        //cudaMemset(escape_count_ssaa_d, 0, buffer_size_ssaa);
        //cudaMemset(magnitude_ssaa_d, 0, buffer_size_ssaa);
        //cudaMemset(escape_count_d, 0, buffer_size);
        cudaCall(cudaMemcpyToSymbol, escape_count_ssaa, &escape_count_ssaa_d, sizeof(double*));
        cudaCall(cudaMemcpyToSymbol, magnitude_ssaa, &magnitude_ssaa_d, sizeof(double*));
        cudaCall(cudaMemcpyToSymbol, escape_count, &escape_count_d, sizeof(double*));

        printf("CUDA initialized.\n");
        return true;
    }

    bool DestroyCUDA() {
        cudaCall(cudaFree, escape_count_ssaa_d);
        cudaCall(cudaFree, magnitude_ssaa_d);
        cudaCall(cudaFreeHost, escape_count_d);

        printf("CUDA destroyed.\n");
        return true;
    }

    template<class T, typename... A>
    float launch_kernel(T& kernel, dim3 work, A&&... args)
    {
        int device;
        cudaDeviceProp props;
        cudaGetDevice(&device);
        cudaGetDeviceProperties(&props, device);

        int threadBlocks;
        if (props.major == 2)
        {
            threadBlocks = 8;
        }
        else if (props.major == 3)
        {
            threadBlocks = 16;
        }
        else
        {
            threadBlocks = 32;
        }

        int blockSize;
        std::uint32_t minGridSize;
        cudaOccupancyMaxPotentialBlockSize((int*)&minGridSize, &blockSize, kernel, 0, 0);

        int maxActiveBlocks = 0;
        do
        {
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxActiveBlocks, kernel, blockSize, 0);

            if (blockSize < props.warpSize || maxActiveBlocks >= threadBlocks)
            {
                break;
            }

            blockSize -= props.warpSize;
        } while (true);

        int blockSizeDimX, blockSizeDimY;
        blockSizeDimX = blockSizeDimY = (int)pow(2, ceil(log(sqrt(blockSize)) / log(2)));

        while (blockSizeDimX * blockSizeDimY > blockSize)
        {
            blockSizeDimY--;
        }

        dim3 block(blockSizeDimX, blockSizeDimY);
        dim3 grid((work.x + block.x - 1) / block.x, (work.y + block.y - 1) / block.y);
        grid.x = grid.x > minGridSize ? grid.x : minGridSize;
        grid.y = grid.y > minGridSize ? grid.y : minGridSize;


        float occupancy = (maxActiveBlocks * blockSize / props.warpSize) / (float)(props.maxThreadsPerMultiProcessor / props.warpSize);

        std::cout << "Grid of size " << grid.x * grid.y << std::endl;
        std::cout << "Launched blocks of size " << blockSize << std::endl;
        std::cout << "Theoretical occupancy " << occupancy * 100.0f << "%" << std::endl;


        cudaEvent_t start;
        cudaCall(cudaEventCreate, &start);

        cudaEvent_t stop;
        cudaCall(cudaEventCreate, &stop);

        cudaCall(cudaEventRecord, start, 0);

        kernel << < grid, block >> > (std::forward<A>(args)...);

        cudaCall(cudaGetLastError);
        cudaCall(cudaEventRecord, stop, 0);
        cudaCall(cudaEventSynchronize, stop);

        float elapsed_time;
        cudaCall(cudaEventElapsedTime, &elapsed_time, start, stop);

        cudaCall(cudaEventDestroy, start);
        cudaCall(cudaEventDestroy, stop);

        cudaProfilerStop();

        return elapsed_time;
    }

    std::uint64_t generate_mandelbrot(double* escape_count_ret, Properties* properties)
    {
        const int width = properties->Get_General_Width();
        const int height = properties->Get_General_Height();
        const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
        const int max_iterations = properties->Get_Mandelbrot_MaxIterations();
        double x_min = properties->Get_Mandelbrot_Point_Cr() - properties->Get_Mandelbrot_Point_Rr();
        double x_max = properties->Get_Mandelbrot_Point_Cr() + properties->Get_Mandelbrot_Point_Rr();
        double y_min = properties->Get_Mandelbrot_Point_Ci() - properties->Get_Mandelbrot_Point_Rr() * ((double)height / width);
        double y_max = properties->Get_Mandelbrot_Point_Ci() + properties->Get_Mandelbrot_Point_Rr() * ((double)height / width);
        const int smooth_shading = properties->Get_Mandelbrot_SmoothShading();

        const int width_ssaa = width * ssaa;
        const int height_ssaa = height * ssaa;

        const double x_scale = (x_max - x_min) / width_ssaa;
        const double y_scale = (y_max - y_min) / height_ssaa;

        const int size_d = width * height * sizeof(double);

        float elapsed_time = launch_kernel(mandelbrot_kernel, 
            dim3(width_ssaa, height_ssaa),
            escape_count_ssaa_d,
            magnitude_ssaa_d,
            max_iterations, 
            width_ssaa,
            height_ssaa,
            x_scale,
            y_scale,
            x_min,
            y_min);

        elapsed_time += launch_kernel(downsample,
            dim3(width, height),
            escape_count_ssaa_d,
            magnitude_ssaa_d,
            escape_count_d,
            max_iterations,
            width,
            height,
            ssaa,
            smooth_shading);

        cudaCall(cudaMemcpy, escape_count_ret, escape_count_d, size_d, cudaMemcpyDeviceToHost);

        return static_cast<std::uint64_t>(elapsed_time);
    }
}