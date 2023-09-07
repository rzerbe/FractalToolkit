#pragma once
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

class Properties
{
public:
	Properties();
	~Properties();

	bool LoadProperties(std::string filename);
	void LoadDefaults();
	nlohmann::json data;

	int Get_General_Width();
	int Get_General_Height();
	int Get_General_Threads();
	int Get_General_InstructionSet();

	double Get_Mandelbrot_Point_Ci();
	double Get_Mandelbrot_Point_Cr();
	double Get_Mandelbrot_Point_Rr();

	int Get_Mandelbrot_MaxIterations();
	int Get_Mandelbrot_ShadingMode();
	int Get_Mandelbrot_SmoothShading();
	int Get_Mandelbrot_SuperSamplingAntiAliasing();
	std::string Get_Mandelbrot_Filename();
	std::string Get_Mandelbrot_Palette();

	void Set_Mandelbrot_Window_X_Min(double val);
	void Set_Mandelbrot_Window_X_Max(double val);
	void Set_Mandelbrot_Window_Y_Min(double val);
	void Set_Mandelbrot_Window_Y_Max(double val);

	void Set_Mandelbrot_Point_Ci(double val);
	void Set_Mandelbrot_Point_Cr(double val);
	void Set_Mandelbrot_Point_Rr(double val);
};

