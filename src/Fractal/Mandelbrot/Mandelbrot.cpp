#include "Mandelbrot.h"

Mandelbrot::Mandelbrot()
{
}

Mandelbrot::Mandelbrot(Properties* properties, Render* render)
	: properties(properties), render(render)
{
	const std::string palette = properties->Get_Mandelbrot_Palette();
	const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();

	const int pixels = width * height;
	const int sub_pixels = ssaa * ssaa;

	this->palette = new Image(palette);
	this->escape_count = new double[pixels];
	this->escape_count_ssaa = new double[pixels * sub_pixels];
	this->magnitude = new double[pixels * sub_pixels];

	this->image = new Image(width, height);
	std::vector<Image::pixel_t*> image_rows(height);
	image_rows[0] = &image->data[0];
	for (std::int32_t i = 1; i < height; i++)
	{
		image_rows[i] = image_rows[i - 1] + width;
	}

	mandelbrot_cuda::InitCUDA(properties);


	// Binds
	std::function<void(InputControl*, int, int, int)> MouseButtonCallback = [=](auto self, int button, int action, int mods)
	{
		const int width = properties->Get_General_Width();
		const int height = properties->Get_General_Height();

		double x_mouse, y_mouse;
		glfwGetCursorPos(render->window, &x_mouse, &y_mouse);

		const double p_cr = properties->Get_Mandelbrot_Point_Cr();
		const double p_ci = properties->Get_Mandelbrot_Point_Ci();
		const double p_rr = properties->Get_Mandelbrot_Point_Rr();

		const double cr_min = p_cr - p_rr;
		const double cr_max = p_cr + p_rr;
		const double ci_min = p_ci - p_rr * ((double)height / width);
		const double ci_max = p_ci + p_rr * ((double)height / width);

		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)	// Left Mouse Button
		{
			double cr_center = (cr_max + cr_min) / 2;
			double ci_center = (ci_max + ci_min) / 2;
			double cr_window = (cr_max - cr_min);
			double ci_window = (ci_max - ci_min);
			std::cout << "MouseUp on " << x_mouse << " " << y_mouse << std::endl;

			properties->Set_Mandelbrot_Point_Cr(p_cr + cr_window * (x_mouse - (double)width / 2) / width);
			properties->Set_Mandelbrot_Point_Ci(p_ci + ci_window * (y_mouse - (double)height / 2) / height);

			//memset(calculation.updatePixel, true, config.resolutionX * config.resolutionY);
			//calculation.screenStill = 0;
			fractal_render(1);
			printf("Window is Cr:[%9e, %9e] Ci:[%9e, %9e] Rr:[%9e] \n",
				cr_min, cr_max,
				ci_min, ci_max,
				p_rr
			);
		}
	};

	std::function<void(InputControl*, double, double)> ScrollCallback = [=](auto self, double x, double y)
	{
		double scroll_factor = 1.05;
		if (y > 0) // Scroll forward
		{
			scroll_factor = 1 / scroll_factor;
		}
		else if (y < 0) // Scroll back
		{
			scroll_factor;
		}

		const int width = properties->Get_General_Width();
		const int height = properties->Get_General_Height();

		const double p_cr = properties->Get_Mandelbrot_Point_Cr();
		const double p_ci = properties->Get_Mandelbrot_Point_Ci();
		double p_rr = properties->Get_Mandelbrot_Point_Rr();

		const double cr_min = p_cr - p_rr;
		const double cr_max = p_cr + p_rr;
		const double ci_min = p_ci - p_rr * ((double)height / width);
		const double ci_max = p_ci + p_rr * ((double)height / width);

		p_rr *= scroll_factor;
		properties->Set_Mandelbrot_Point_Rr(p_rr);

		//memset(calculation.updatePixel, true, config.resolutionX * config.resolutionY);
		//calculation.screenStill = 0;
		fractal_render(1);
		printf("Window is Cr:[%9e, %9e] Ci:[%9e, %9e] Rr:[%9e] \n",
			cr_min, cr_max,
			ci_min, ci_max,
			p_rr
		);
	};


	render->inputControl->MouseButtonCallback = MouseButtonCallback;
	render->inputControl->ScrollCallback = ScrollCallback;
	render->screenBufferCPU = (unsigned char*) &image->data[0];
}

Mandelbrot::~Mandelbrot()
{
	delete &palette;
	delete[] escape_count;
	delete[] escape_count_ssaa;
	delete[] magnitude;
}

void Mandelbrot::fractal_render(int output)
{
	auto begin = std::chrono::high_resolution_clock::now();

	const int threads = properties->Get_General_Threads();
	const int instruction_set = properties->Get_General_InstructionSet();
	const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();

	switch (instruction_set)
	{
	case 0:
		c_fractal_64(0, 0, width * ssaa, height * ssaa);
		downsample();
		break;
	case 1:
		avx_fractal_64(0, 0, width * ssaa, height * ssaa);
		downsample();
		break;
	case 2:
		avx_fractal_32(0, 0, width * ssaa, height * ssaa);
		downsample();
		break;
	case 3:
		std::uint64_t elapsed_time = cuda_fractal();
		std::cout << "Generate time: " << elapsed_time << " ms." << std::endl;
		break;
	}

	// Post Processing
	write_buffer();

	render->screenBufferCPU = (unsigned char*) &image->data[0];
	render->redraw = (bool*)true;

	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("Fractal Generate Time %d\n", elasped);
}

std::uint64_t Mandelbrot::cuda_fractal() 
{
	auto begin = std::chrono::high_resolution_clock::now();

	uint64_t r = mandelbrot_cuda::generate_mandelbrot(escape_count, properties);

	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("CUDA Thread finished in %d milliseconds \n", elasped);

	return r;
}

void Mandelbrot::avx_fractal_32(int i_init, int j_init, int i_fin, int j_fin)
{
	auto begin = std::chrono::high_resolution_clock::now();

	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();
	const int threads = properties->Get_General_Threads();
	const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
	const int max_iterations = properties->Get_Mandelbrot_MaxIterations();

	const double x_min = properties->Get_Mandelbrot_Point_Cr() - properties->Get_Mandelbrot_Point_Rr();
	const double x_max = properties->Get_Mandelbrot_Point_Cr() + properties->Get_Mandelbrot_Point_Rr();
	const double y_min = properties->Get_Mandelbrot_Point_Ci() - properties->Get_Mandelbrot_Point_Rr() * ((double)height / width);
	const double y_max = properties->Get_Mandelbrot_Point_Ci() + properties->Get_Mandelbrot_Point_Rr() * ((double)height / width);

	const int width_ssaa = width * ssaa;
	const int height_ssaa = height * ssaa;

	__m256 xmin = _mm256_set1_ps(x_min);
	__m256 ymin = _mm256_set1_ps(y_min);
	__m256 xscale = _mm256_set1_ps((x_max - x_min) / width_ssaa);
	__m256 yscale = _mm256_set1_ps((y_max - y_min) / height_ssaa);
	__m256 threshold = _mm256_set1_ps(4);
	__m256 one = _mm256_set1_ps(1);
	__m256 zero25 = _mm256_set1_ps(0.25);
	__m256 zero0625 = _mm256_set1_ps(0.0625);
	__m256 zero = _mm256_set1_ps(0);
	__m256 all_one = _mm256_set_ps(1, 1, 1, 1, 1, 1, 1, 1);

	#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
	for (int y = 0; y < height_ssaa; y++)
	{
		for (int x = 0; x < width_ssaa; x += 8)
		{
			__m256 mx = _mm256_set_ps(
				x + 7, x + 6, x + 5, x + 4,
				x + 3, x + 2, x + 1, x + 0);
			__m256 my = _mm256_set1_ps(y);
			__m256 cr = _mm256_add_ps(_mm256_mul_ps(mx, xscale), xmin);
			__m256 ci = _mm256_add_ps(_mm256_mul_ps(my, yscale), ymin);


			//check whether point lies in main cardioid or check whether point lies in period-2 bulb
			__m256 ci_sqr = _mm256_mul_ps(ci, ci);
			__m256 cr_minus_zero25 = _mm256_sub_ps(cr, zero25);
			__m256 cr_minus_zero25_sqr = _mm256_mul_ps(cr_minus_zero25, cr_minus_zero25);

			__m256 cr_plus_1 = _mm256_add_ps(cr, one);
			__m256 cr_plus_1_sqr = _mm256_mul_ps(cr_plus_1, cr_plus_1);

			// (cr - 0.25) * (cr - 0.25) + ci * ci;
			__m256 q = _mm256_add_ps(cr_minus_zero25_sqr, ci_sqr);
			// (q * (q + cr - 0.25) - 0.25 * ci * ci < 0 
			__m256 cardioid = _mm256_sub_ps((_mm256_mul_ps(q, (_mm256_add_ps(q, cr_minus_zero25)))), _mm256_mul_ps(zero25, ci_sqr));
			// (cr + 1) * (cr + 1) + ci * ci - 0.0625 < 0
			__m256 period_2 = _mm256_sub_ps(_mm256_add_ps(cr_plus_1_sqr, ci_sqr), zero0625);

			__m256 bitmask = _mm256_or_ps(_mm256_cmp_ps(cardioid, zero, _CMP_LT_OS), _mm256_cmp_ps(period_2, zero, _CMP_LT_OS));

			double* ecs = escape_count_ssaa + y * width_ssaa + x;
			double* mag = magnitude + y * width_ssaa + x;

			if (isnan(bitmask.m256_f32[0]) && isnan(bitmask.m256_f32[1]) && isnan(bitmask.m256_f32[2]) && isnan(bitmask.m256_f32[3]) &&
				isnan(bitmask.m256_f32[4]) && isnan(bitmask.m256_f32[5]) && isnan(bitmask.m256_f32[6]) && isnan(bitmask.m256_f32[7])) {
				ecs[0] = max_iterations;
				ecs[1] = max_iterations;
				ecs[2] = max_iterations;
				ecs[3] = max_iterations;
				ecs[4] = max_iterations;
				ecs[5] = max_iterations;
				ecs[6] = max_iterations;
				ecs[7] = max_iterations;

				mag[0] = INFINITY;
				mag[1] = INFINITY;
				mag[2] = INFINITY;
				mag[3] = INFINITY;
				mag[4] = INFINITY;
				mag[5] = INFINITY;
				mag[6] = INFINITY;
				mag[7] = INFINITY;
				continue;
			}

			__m256 zr = zero;
			__m256 zi = zero;
			int k = 1;
			__m256 mk = _mm256_set1_ps(k);
			__m256 mag2 = zero;

			while (++k < max_iterations)
			{
				/* Compute z1 from z0 */
				__m256 zr2 = _mm256_mul_ps(zr, zr);
				__m256 zi2 = _mm256_mul_ps(zi, zi);
				__m256 zrzi = _mm256_mul_ps(zr, zi);
				/* zr1 = zr0 * zr0 - zi0 * zi0 + cr */
				/* zi1 = zr0 * zi0 + zr0 * zi0 + ci */
				zr = _mm256_add_ps(_mm256_sub_ps(zr2, zi2), cr);
				zi = _mm256_add_ps(_mm256_add_ps(zrzi, zrzi), ci);

				/* Increment k */
				zr2 = _mm256_mul_ps(zr, zr);
				zi2 = _mm256_mul_ps(zi, zi);
				mag2 = _mm256_add_ps(zr2, zi2);
				__m256 mask = _mm256_cmp_ps(mag2, threshold, _CMP_LT_OS);
				mk = _mm256_add_ps(_mm256_and_ps(mask, one), mk);

				/* Early bailout? */
				if (_mm256_testz_ps(mask, _mm256_set1_ps(-1)))
				{
					break;
				}
			}

			for (int i = 0; i < 8; i++)
			{
				if (y * width_ssaa + x + i > height_ssaa * width_ssaa)
				{
					break;
				}
				ecs[i] = mk.m256_f32[i];
				mag[i] = mag2.m256_f32[i];
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("AVX_32 finished in %d milliseconds \n", elasped);
}

void Mandelbrot::avx_fractal_64(int i_init, int j_init, int i_fin, int j_fin)
{
	auto begin = std::chrono::high_resolution_clock::now();

	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();
	const int threads = properties->Get_General_Threads();
	const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
	const int max_iterations = properties->Get_Mandelbrot_MaxIterations();

	const double x_min = properties->Get_Mandelbrot_Point_Cr() - properties->Get_Mandelbrot_Point_Rr();
	const double x_max = properties->Get_Mandelbrot_Point_Cr() + properties->Get_Mandelbrot_Point_Rr();
	const double y_min = properties->Get_Mandelbrot_Point_Ci() - properties->Get_Mandelbrot_Point_Rr() * ((double)height / width);
	const double y_max = properties->Get_Mandelbrot_Point_Ci() + properties->Get_Mandelbrot_Point_Rr() * ((double)height / width);

	const int width_ssaa = width * ssaa;
	const int height_ssaa = height * ssaa;

	__m256d xmin = _mm256_set1_pd(x_min);
	__m256d ymin = _mm256_set1_pd(y_min);
	__m256d xscale = _mm256_set1_pd((x_max - x_min) / width_ssaa);
	__m256d yscale = _mm256_set1_pd((y_max - y_min) / height_ssaa);
	__m256d threshold = _mm256_set1_pd(4);
	__m256d one = _mm256_set1_pd(1);
	__m256d zero25 = _mm256_set1_pd(0.25);
	__m256d zero0625 = _mm256_set1_pd(0.0625);
	__m256d zero = _mm256_set1_pd(0);
	__m256d all_one = _mm256_set_pd(1, 1, 1, 1);

	#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
	for (int y = 0; y < height_ssaa; y++) 
	{
		for (int x = 0; x < width_ssaa; x += 4) 
		{
			__m256d mx = _mm256_set_pd(x + 3, x + 2, x + 1, x + 0);
			__m256d my = _mm256_set1_pd(y);
			__m256d cr = _mm256_add_pd(_mm256_mul_pd(mx, xscale), xmin);
			__m256d ci = _mm256_add_pd(_mm256_mul_pd(my, yscale), ymin);


			//check whether point lies in main cardioid or check whether point lies in period-2 bulb
			__m256d ci_sqr = _mm256_mul_pd(ci, ci);
			__m256d cr_minus_zero25 = _mm256_sub_pd(cr, zero25);
			__m256d cr_minus_zero25_sqr = _mm256_mul_pd(cr_minus_zero25, cr_minus_zero25);

			__m256d cr_plus_1 = _mm256_add_pd(cr, one);
			__m256d cr_plus_1_sqr = _mm256_mul_pd(cr_plus_1, cr_plus_1);

			// (cr - 0.25) * (cr - 0.25) + ci * ci;
			__m256d q = _mm256_add_pd(cr_minus_zero25_sqr, ci_sqr);
			// (q * (q + cr - 0.25) - 0.25 * ci * ci < 0 
			__m256d cardioid = _mm256_sub_pd((_mm256_mul_pd(q, (_mm256_add_pd(q, cr_minus_zero25)))), _mm256_mul_pd(zero25, ci_sqr));
			// (cr + 1) * (cr + 1) + ci * ci - 0.0625 < 0
			__m256d period_2 = _mm256_sub_pd(_mm256_add_pd(cr_plus_1_sqr, ci_sqr), zero0625);

			__m256d bitmask = _mm256_or_pd(_mm256_cmp_pd(cardioid, zero, _CMP_LT_OS), _mm256_cmp_pd(period_2, zero, _CMP_LT_OS));

			double* ecs = escape_count_ssaa + y * width_ssaa + x;
			double* mag = magnitude + y * width_ssaa + x;

			if (isnan(bitmask.m256d_f64[0]) && isnan(bitmask.m256d_f64[1]) && isnan(bitmask.m256d_f64[2]) && isnan(bitmask.m256d_f64[3])) {
				ecs[0] = max_iterations;
				ecs[1] = max_iterations;
				ecs[2] = max_iterations;
				ecs[3] = max_iterations;

				mag[0] = INFINITY;
				mag[1] = INFINITY;
				mag[2] = INFINITY;
				mag[3] = INFINITY;
				continue;
			}

			__m256d zr = zero;
			__m256d zi = zero;
			int k = 1;
			__m256d mk = _mm256_set1_pd(k);
			__m256d mag2 = zero;

			while (++k < max_iterations)
			{
				/* Compute z1 from z0 */
				__m256d zr2 = _mm256_mul_pd(zr, zr);
				__m256d zi2 = _mm256_mul_pd(zi, zi);
				__m256d zrzi = _mm256_mul_pd(zr, zi);
				/* zr1 = zr0 * zr0 - zi0 * zi0 + cr */
				/* zi1 = zr0 * zi0 + zr0 * zi0 + ci */
				zr = _mm256_add_pd(_mm256_sub_pd(zr2, zi2), cr);
				zi = _mm256_add_pd(_mm256_add_pd(zrzi, zrzi), ci);

				/* Increment k */
				zr2 = _mm256_mul_pd(zr, zr);
				zi2 = _mm256_mul_pd(zi, zi);
				__m256d mag2 = _mm256_add_pd(zr2, zi2);
				__m256d mask = _mm256_cmp_pd(mag2, threshold, _CMP_LT_OS);
				mk = _mm256_add_pd(_mm256_and_pd(mask, one), mk);

				/* Early bailout? */
				if (_mm256_testz_pd(mask, _mm256_set1_pd(-1)))
				{
					break;
				}
			}

			for (int i = 0; i < 4; i++) 
			{
				ecs[i] = mk.m256d_f64[i];
				mag[i] = mag2.m256d_f64[i];
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("AVX_64 finished in %d milliseconds \n", elasped);
}

void Mandelbrot::c_fractal_64(int i_init, int j_init, int i_fin, int j_fin)
{
	auto begin = std::chrono::high_resolution_clock::now();

	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();
	const int threads = properties->Get_General_Threads();
	const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
	const int max_iterations = properties->Get_Mandelbrot_MaxIterations();

	const double x_min = properties->Get_Mandelbrot_Point_Cr() - properties->Get_Mandelbrot_Point_Rr();
	const double x_max = properties->Get_Mandelbrot_Point_Cr() + properties->Get_Mandelbrot_Point_Rr();
	const double y_min = properties->Get_Mandelbrot_Point_Ci() - properties->Get_Mandelbrot_Point_Rr() * ((double)height/width);
	const double y_max = properties->Get_Mandelbrot_Point_Ci() + properties->Get_Mandelbrot_Point_Rr() * ((double)height/width);

	const int width_ssaa = width * ssaa;
	const int height_ssaa = height * ssaa;

	const double x_scale = (x_max - x_min) / width_ssaa;
	const double y_scale = (y_max - y_min) / height_ssaa;

	#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
	for (int j = j_init; j < j_fin; ++j) // y axis
	{
		double ci = y_min + j * x_scale;

		for (int i = i_init; i < i_fin; ++i) // x axis
		{
			double cr = x_min + i * y_scale;

			//check whether point lies in main cardioid or check whether point lies in period-2 bulb
			double q = (cr - 0.25) * (cr - 0.25) + ci * ci;
			if (q * (q + cr - 0.25) - 0.25 * ci * ci < 0 || (cr + 1) * (cr + 1) + ci * ci - 0.0625 < 0)
			{
				escape_count_ssaa[j * width_ssaa + i] = max_iterations;
				magnitude[j * width_ssaa + i] = INFINITY;
				continue;
			}

			/*if (calculation.screenOptimization == true)
			{
				if (calculation.updatePixel[j * width + i] == false)
				{
					shortcut = true;
				}
			}*/

			double zr, zi, zMagSqr, zr_prev;
			int k;

			zr = zi = zMagSqr = zr_prev = 0;

			for (k = 1; k < max_iterations; ++k)
			{
				zr_prev = zr;
				zr = zr * zr - zi * zi + cr;
				zi = 2 * zr_prev * zi + ci;

				zMagSqr = zr * zr + zi * zi;
				if (zMagSqr > 4.0)
				{
					break;
				}
			}
			escape_count_ssaa[j * width_ssaa + i] = k;
			magnitude[j * width_ssaa + i] = zMagSqr;
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("Standard finished in %d milliseconds \n", elasped);
}

void Mandelbrot::write_buffer()
{
	auto begin = std::chrono::high_resolution_clock::now();

	const int threads = properties->Get_General_Threads();
	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();

	//image.data.resize(width * height);

	write_region(0, 0, width, height);

	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

	printf("Write Buffer Time %d milliseconds\n",
		elasped
	);
}

void Mandelbrot::write_region(int x_init, int y_init, int x_fin, int y_fin)
{
	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();
	const int threads = properties->Get_General_Threads();
	const int max_iterations = properties->Get_Mandelbrot_MaxIterations();
	const int shading_mode = properties->Get_Mandelbrot_ShadingMode();

	const int colord = palette->data.size(); // update this
	const int colorn = palette->data.size();

	#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
	for (int y = y_init; y < y_fin; y++)
	{
		for (int x = x_init; x < x_fin; x++)
		{
			int pixel = y * width + x;
			if (escape_count[pixel] == max_iterations || isnan(escape_count[pixel]))
			{
				image->data[pixel].r = 0;
				image->data[pixel].g = 0;
				image->data[pixel].b = 0;
				image->data[pixel].a = 255;
			}
			else
			{
				//Linear Interpolation
				if (shading_mode == 0)
				{
					double k = escape_count[pixel];

					double kscaled = (k - colord * (int)(k / colord)) / colord * colorn;

					Image::pixel_t color_start = palette->data[(int)kscaled];
					Image::pixel_t color_end = palette->data[((int)kscaled + 1) % colorn];

					image->data[pixel].r = (int)((color_end.r - color_start.r) * (kscaled - (int)(kscaled)) + color_start.r);
					image->data[pixel].g = (int)((color_end.g - color_start.g) * (kscaled - (int)(kscaled)) + color_start.g);
					image->data[pixel].b = (int)((color_end.b - color_start.b) * (kscaled - (int)(kscaled)) + color_start.b);

					//image->data[pixel].r = color_start.r;
					//image->data[pixel].g = color_start.g;
					//image->data[pixel].b = color_start.b;


					image->data[pixel].a = 255;

				}
				//Monochromatic
				else if (shading_mode == 1)
				{
					image->data[pixel].r = max_iterations - escape_count[y * width + x];
					image->data[pixel].g = max_iterations - escape_count[y * width + x];
					image->data[pixel].b = max_iterations - escape_count[y * width + x];
					image->data[pixel].a = 255;
				}
			}
		}
	}
}

void Mandelbrot::downsample()
{
	auto begin = std::chrono::high_resolution_clock::now();

	int width = properties->Get_General_Width();
	int height = properties->Get_General_Height();

	downsample_region( 0, 0, width, height);

	auto end = std::chrono::high_resolution_clock::now();
	auto elasped = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

	printf("Downsampling Time %d milliseconds\n",
		elasped
	);
}

void Mandelbrot::downsample_region(int i_init, int j_init, int i_fin, int j_fin)
{

	const int width = properties->Get_General_Width();
	const int height = properties->Get_General_Height();
	const int threads = properties->Get_General_Threads();
	const int ssaa = properties->Get_Mandelbrot_SuperSamplingAntiAliasing();
	const int max_iterations = properties->Get_Mandelbrot_MaxIterations();
	const int shading_mode = properties->Get_Mandelbrot_ShadingMode();
	const int smooth_shading = properties->Get_Mandelbrot_SmoothShading();

	// Pixel
	int index;
	int smallestIndex;
	double min;
	double norm;
	#pragma omp parallel for schedule(dynamic, 1) private(index, smallestIndex, min, norm) num_threads(threads)
	for (int j = j_init; j < j_fin; ++j)
	{
		for (int i = i_init; i < i_fin; ++i)
		{

			min = max_iterations;
			smallestIndex = ((j * ssaa + 0) * width * ssaa) + (i * ssaa + 0);
			// Sub Pixel
			for (int y = 0; y < ssaa; ++y)
			{
				for (int x = 0; x < ssaa; ++x)
				{

					index = ((j * ssaa + y) * width * ssaa) + (i * ssaa + x);
					if (escape_count_ssaa[index] < min)
					{
						smallestIndex = index;
						min = escape_count_ssaa[index];
					}
				}
			}

			if (smooth_shading == 1)
			{
				//norm = min + 1 - log(log(sqrt(magnitude[smallestIndex]))) / log(2);

				double logzmag = abs(log(magnitude[smallestIndex])) / 2;
				double potential = log(logzmag / log(2)) / log(2);
				norm = min + 1 - potential;

				if (norm > max_iterations || norm < -1000 || min == max_iterations - 1)
				{
					norm = max_iterations;
				}

				escape_count[j * width + i] = norm;
			}
			else
			{
				escape_count[j * width + i] = min;
			}
		}
	}
}