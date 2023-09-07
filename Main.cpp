#include "Main.h"

#include "Program.h"
#include "Properties.h"

int main(int argc, char** argv)
{
    int error_code;

    try
    {
        Properties properties;
        properties.LoadProperties("properties.json");

        Program program(properties);
		error_code = program.run();
    }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;

        error_code = -1;
    }

    return error_code;
}