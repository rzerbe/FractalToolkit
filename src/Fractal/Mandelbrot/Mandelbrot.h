#pragma once
#include <vector>
#include <map>
#include <chrono>
#include <thread>
#include <immintrin.h>
#include <intrin.h>
#include <omp.h>

#include "MandelbrotCuda.h"
#include "../../Program/Image.h"
#include "../../Program/Render.h"
#include "../../Program/Properties.h"
#include "../../Program/InputControl.h"
class Mandelbrot
{
public:
	Mandelbrot();
	Mandelbrot(Properties* properties, Render* render);

	~Mandelbrot();

	void fractal_render(int output);
	std::uint64_t cuda_fractal();

	Image* image;

	std::function<void(InputControl*, int, int, int)> onMouseClickxs;

private:
	Properties* properties;
	Render* render;

	Image* palette;

	double* escape_count;
	double* escape_count_ssaa;
	double* magnitude;

	void avx_fractal_32(int i_init, int j_init, int i_fin, int j_fin);
	void avx_fractal_64(int i_init, int j_init, int i_fin, int j_fin);
	void c_fractal_64(int i_init, int j_init, int i_fin, int j_fin);

	void write_buffer();
	void write_region(int x_init, int y_init, int x_fin, int y_fin);
	void downsample();
	void downsample_region(int i_init, int j_init, int i_fin, int j_fin);
};

