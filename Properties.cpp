#include "Properties.h"

Properties::Properties() {
}

Properties::~Properties() {
}

bool Properties::LoadProperties(std::string filename) {
	try 
	{
		std::ifstream ifstream(filename);
		this->data = nlohmann::json::parse(ifstream);
	}
	catch (std::exception& ex)
	{
		return false;
	}

	return true;
}

void Properties::LoadDefaults() {
	data["General"]["Width"] = 1024;
	data["General"]["Height"] = 1024;
	data["General"]["Threads"] = 1;
	data["General"]["Instruction_Set"] = 0;

	data["Mandelbrot"]["Window"]["X"]["Min"] = -1.00;
	data["Mandelbrot"]["Window"]["X"]["Max"] = 1.00;
	data["Mandelbrot"]["Window"]["Y"]["Min"] = -1.00;
	data["Mandelbrot"]["Window"]["Y"]["Max"] = 1.00;

	data["Mandelbrot"]["Max_Iterations"] = 256;
	data["Mandelbrot"]["Shading_Mode"] = 0;
	data["Mandelbrot"]["Smooth_Shading"] = 0;

	data["Mandelbrot"]["Super_Sampling_Anti_Aliasing"] = 1;
	data["Mandelbrot"]["Filename"] = "mandelbrot.bmp";
	data["Mandelbrot"]["Palette"] = "palette.json";
}

int Properties::Get_General_Width()
{
	return data["General"]["Width"];
}
int Properties::Get_General_Height() 
{
	return data["General"]["Height"];
}
int Properties::Get_General_Threads() 
{
	return data["General"]["Threads"];
}
int Properties::Get_General_InstructionSet() 
{
	return data["General"]["Instruction_Set"];
}

double Properties::Get_Mandelbrot_Point_Ci()
{
	return data["Mandelbrot"]["Point"]["Ci"];
}
double Properties::Get_Mandelbrot_Point_Cr()
{
	return data["Mandelbrot"]["Point"]["Cr"];
}
double Properties::Get_Mandelbrot_Point_Rr()
{
	return data["Mandelbrot"]["Point"]["Rr"];
}

int Properties::Get_Mandelbrot_MaxIterations()
{
	return data["Mandelbrot"]["Max_Iterations"];
}
int Properties::Get_Mandelbrot_ShadingMode()
{
	return data["Mandelbrot"]["Shading_Mode"];
}
int Properties::Get_Mandelbrot_SmoothShading()
{
	return data["Mandelbrot"]["Smooth_Shading"];
}
int Properties::Get_Mandelbrot_SuperSamplingAntiAliasing()
{
	return data["Mandelbrot"]["Super_Sampling_Anti_Aliasing"];
}
std::string Properties::Get_Mandelbrot_Filename()
{
	return data["Mandelbrot"]["Filename"];
}
std::string Properties::Get_Mandelbrot_Palette()
{
	return data["Mandelbrot"]["Palette"];
}

void Properties::Set_Mandelbrot_Window_X_Min(double val)
{
	data["Mandelbrot"]["Window"]["X"]["Min"] = val;
}
void Properties::Set_Mandelbrot_Window_X_Max(double val)
{
	data["Mandelbrot"]["Window"]["X"]["Max"] = val;
}
void Properties::Set_Mandelbrot_Window_Y_Min(double val)
{
	data["Mandelbrot"]["Window"]["Y"]["Min"] = val;
}
void Properties::Set_Mandelbrot_Window_Y_Max(double val)
{
	data["Mandelbrot"]["Window"]["Y"]["Max"] = val;
}

void Properties::Set_Mandelbrot_Point_Ci(double val)
{
	data["Mandelbrot"]["Point"]["Ci"] = val;
}
void Properties::Set_Mandelbrot_Point_Cr(double val)
{
	data["Mandelbrot"]["Point"]["Cr"] = val;
}
void Properties::Set_Mandelbrot_Point_Rr(double val)
{
	data["Mandelbrot"]["Point"]["Rr"] = val;
}