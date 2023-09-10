#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <memory>
#include <sstream>
#include <iostream>
#include <numeric>

#include "cuda_occupancy.h"
#include "cuda_runtime.h"
#include "cuda_profiler_api.h"

#include "../../Program/Image.h"
#include "../../Program/Properties.h"

namespace mandelbrot_cuda
{
    bool InitCUDA(Properties* properties);
    bool DestroyCUDA();
    std::uint64_t generate_mandelbrot(double* escape_count, Properties* properties);
}