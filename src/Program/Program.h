#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "Cuda.h"
#include "Image.h"
#include "../Fractal/Mandelbrot/Mandelbrot.h"
#include "Properties.h"
#include "InputControl.h"
#include "Render.h"

// CPU
#include <chrono>
#include <thread>
#include <immintrin.h>
#include <intrin.h>


class Program
{
public:

    Program(Properties properties);

    ~Program();

    void display_help();

    int run();

private:
    Properties properties;
};