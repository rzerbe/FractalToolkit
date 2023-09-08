#pragma once

#include <stdio.h>
#include "Properties.h"

// Renderer
#include <GLEW/glew.h>
#include <GLFW/glfw3.h>

class InputControl {
public:
    InputControl();
    InputControl(Properties* properties);
    ~InputControl();
    void bind(GLFWwindow* applicationWindow, Properties* properties);
    void bindErrorFunction();
    void unbind();

    std::function<void(InputControl*, int, int, int)> MouseButtonCallback = [](auto self, int, int, int) { printf("Mouse Click Unbound"); };
    std::function<void(InputControl*, double, double)> ScrollCallback = [](auto self, double, double) { printf("Scroll Wheel Unbound"); };

    Properties* properties = nullptr;
private:

    static void errorCallbackStatic(int err,
        const char* msg);
};