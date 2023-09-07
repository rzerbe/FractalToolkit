#include "Program.h"

Program::Program(Properties properties)
    : properties(properties)
{
}

Program::~Program()
{
}

void Program::display_help()
{
    printf("\nOpenGL Renderer Commands\n");
    printf("[left mouse] Change center\n");
    printf("[mouse wheel up] Zoom in\n");
    printf("[mouse wheel down] Zoom out\n");

    printf("Terminal Mode Commands\n");
    printf("\"?, help\" \t\t show menu items.\n");
    printf("\"point cr ci rr\" \t move the view to a specified point using real, imaginary, and magnification.\n");
    printf("\"save [name.bmp]\" \t save the view to a bmp file. If [name.bmp] is not provided, default from config will be used\n");
    return;
}

int Program::run()
{
    int width = properties.Get_General_Width();
    int height = properties.Get_General_Height();
    
    std::thread renderWindow;

    InputControl* inputControl = new InputControl();
    Render* render = new Render(&properties, inputControl);
    Mandelbrot mandelbrot(&properties, render);

    std::string line;
    std::string command;
    std::string field;

    renderWindow = std::thread(&Render::render_loop, render);

    display_help();
    do {
        std::getline(std::cin, line);
        command = line.substr(0, line.find(' '));
        field = line.substr(line.find(' ') + 1, line.find('\n'));

        if (command.compare("?") == 0 || command.compare("help") == 0)
        {
            display_help();
        }
        else if (command.compare("point") == 0) {
            double cr, ci, rr;

            std::stringstream stream(field);
            stream >> cr;
            stream.ignore();
            stream >> ci;
            stream.ignore();
            stream >> rr;
            stream.ignore();

            properties.Set_Mandelbrot_Point_Cr(cr);
            properties.Set_Mandelbrot_Point_Ci(ci);
            properties.Set_Mandelbrot_Point_Rr(rr);
            mandelbrot.fractal_render(1);
        }
        else if (command.compare("save") == 0) {
            std::string filename;

            std::stringstream stream(field);
            stream >> filename;
            stream.ignore();

            if (filename.empty()) {
                filename = properties.Get_Mandelbrot_Filename();
            }

            Image* image = mandelbrot.image;
            std::ofstream fout(filename, std::ios::binary);
            fout << image;

            fout.close();
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (true);

    renderWindow.join();
    exit(EXIT_SUCCESS);
   
    return 0;
}