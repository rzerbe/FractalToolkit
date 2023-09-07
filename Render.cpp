#include "Render.h"

Render::Render() 
{
}

Render::Render(Properties* properties, InputControl* inputControl)
{
    this->properties = properties;
    this->inputControl = inputControl;
    this->screenBufferCPU = NULL;
    this->redraw = redraw;
}

Render::~Render()
{
}

void Render::initialize_GL()
{
    inputControl->bindErrorFunction();
    if (!glfwInit())
    {
        exit(EXIT_FAILURE);
    }

    int width = properties->Get_General_Width();
    int height = properties->Get_General_Height();

    window = glfwCreateWindow(width, height, "FractalGL", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // Set GLFW callback functions
    inputControl->bind(window, properties);

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glewExperimental = GL_TRUE;
    glewInit();

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glViewport(0, 0, width, height);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    const GLubyte* renderer = glGetString(GL_RENDERER); // get renderer string
    const GLubyte* version = glGetString(GL_VERSION); // version as a string
    printf("\nRenderer: \t\t%s\n", renderer);
    printf("OpenGL Version: \t%s\n", version);
}

GLuint Render::compileASMShader(GLenum program_type, const char* code)
{
    GLuint program_id;
    glGenProgramsARB(1, &program_id);
    glBindProgramARB(program_type, program_id);
    glProgramStringARB(program_type, GL_PROGRAM_FORMAT_ASCII_ARB, (GLsizei)strlen(code), (GLubyte*)code);

    GLint error_pos;
    glGetIntegerv(GL_PROGRAM_ERROR_POSITION_ARB, &error_pos);

    if (error_pos != -1)
    {
        const GLubyte* error_string;
        error_string = glGetString(GL_PROGRAM_ERROR_STRING_ARB);
        fprintf(stderr, "Program error at position: %d\n%s\n", (int)error_pos, error_string);
        return 0;
    }

    return program_id;
}

void Render::display_func()
{
    int width = properties->Get_General_Width();
    int height = properties->Get_General_Height();

    glBindTexture(GL_TEXTURE_2D, gl_Tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, (GLvoid*)screenBufferCPU);

    glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, gl_Shader);
    glEnable(GL_FRAGMENT_PROGRAM_ARB);
    glDisable(GL_DEPTH_TEST);

    glMapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, GL_WRITE_ONLY_ARB);
    glBegin(GL_QUADS);

    // All verticies are flipped to start at upper left corner
    glTexCoord2f(0.0f, 0.0f);
    glVertex2f(0.0f, 1.0f);

    glTexCoord2f(1.0f, 0.0f);
    glVertex2f(1.0f, 1.0f);

    glTexCoord2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.0f);

    glTexCoord2f(0.0f, 1.0f);
    glVertex2f(0.0f, 0.0f);
    glEnd();

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
    glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB);
    //glDisable(GL_FRAGMENT_PROGRAM_ARB);

    glfwSwapBuffers(window);
}

bool Render::initialize_buffers(int width, int height, bool display)
{
    // Flush buffers
    if (display == true)
    {
        //delete_buffers();
    }

    // Check for minimized window
    if ((width == 0) && (height == 0))
    {
        return false;
    }

    // Allocate Buffers
    //calculation.escapeBufferCPU = new double[width * height];
    //calculation.escapeBufferSuperSampling = new double[width * height * config.SSAA * config.SSAA];
    //calculation.magnitude = new double[width * height * config.SSAA * config.SSAA];
    //screenBufferCPU = new unsigned char[width * height * 4];
    //calculation.updatePixel = new bool[width * height];
    //memset(calculation.updatePixel, true, width * height);
    std::cout << "Resolution Set to " << width << " " << height << std::endl;
    //std::cout << "Buffer Set to " << width * config.SSAA << " " << height * config.SSAA << std::endl;

    glGenTextures(1, &gl_Tex);
    glBindTexture(GL_TEXTURE_2D, gl_Tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_BGRA, GL_UNSIGNED_BYTE, (GLvoid*)screenBufferCPU);
    glEnable(GL_TEXTURE_2D);

    glGenBuffers(1, &gl_PBO);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_PBO);
    glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, width * height * sizeof(Image::pixel_t), (GLvoid*)screenBufferCPU, GL_STREAM_COPY);

    gl_Shader = compileASMShader(GL_FRAGMENT_PROGRAM_ARB, "!!ARBfp1.0\n"
        "TEX result.color, fragment.texcoord, texture[0], 2D; \n"
        "END");

    if (display == true) // Display
    {
        display_func();
        display_func();
    }
    return true;
}

void Render::render_loop()
{
    int width = properties->Get_General_Width();
    int height = properties->Get_General_Height();

    initialize_GL();
    initialize_buffers(width, height, true);

    do // Main Loop
    {
        if (redraw)
        {
            display_func();
            redraw = false;
        }
        glfwPollEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (!glfwWindowShouldClose(window));
    glfwDestroyWindow(window);
    glfwTerminate();
    //delete_buffers();
}