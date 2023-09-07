#pragma once

#include <thread>

#include "InputControl.h"
#include "Image.h"

class Render
{
public:
	Render();
	Render(Properties* properties, InputControl* inputControl);
	~Render();
	void Render::initialize_GL();
	GLuint Render::compileASMShader(GLenum program_type, const char* code);
	void Render::display_func();
	bool Render::initialize_buffers(int width, int height, bool display);
	void Render::render_loop();
	void Render::vao_init();
	void Render::vao_exit();
	void Render::vao_draw(GLuint mode);
	void Render::drawCube();
	unsigned char* screenBufferCPU;
	bool redraw;
	GLFWwindow* window;
	InputControl* inputControl;

	GLuint vbo[4] = { (GLuint)-1,(GLuint)-1,(GLuint)-1,(GLuint)-1 };
	GLuint vao[4] = { (GLuint)-1,(GLuint)-1,(GLuint)-1,(GLuint)-1 };
	const float vao_pos[9] =
	{
		//       x      y     z
			 0.75f, 0.75f, 0.0f,
			 0.75f,-0.75f, 0.0f,
			-0.75f,-0.75f, 0.0f,
	};
	const float vao_col[9] =
	{
		//      r   g    b
			 1.0f,0.0f,0.0f,
			 0.0f,1.0f,0.0f,
			 0.0f,0.0f,1.0f,
	};

private:
	Properties* properties;

	GLuint gl_PBO, gl_Tex, gl_Shader;
};

