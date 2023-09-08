#include "InputControl.h"

#define genericCallback(functionName)\
            [](GLFWwindow* window, auto... args) {\
                auto pointer = static_cast<InputControl*>(glfwGetWindowUserPointer(window));\
                if (pointer->functionName) pointer->functionName(pointer, args...);\
            }

InputControl::InputControl() {

}

InputControl::InputControl(Properties* properties)
	: properties(properties)
{

}

InputControl::~InputControl() {

}

void InputControl::bind(GLFWwindow* applicationWindow, Properties* properties) {
	/*/glfwSetKeyCallback(applicationWindow, keyCallbackStatic);
	glfwSetMouseButtonCallback(applicationWindow, mouseButtonCallbackStatic);
	glfwSetCursorPosCallback(applicationWindow, mouseCursorPosCallbackStatic);
	glfwSetWindowSizeCallback(applicationWindow, windowSizeCallbackStatic);
	glfwSetScrollCallback(applicationWindow, scrollCallbackStatic);*/

	properties->data["General"]["Width"];
	properties->data["General"]["Height"];
	this->properties = properties;

	glfwSetWindowUserPointer(applicationWindow, this);

	glfwSetMouseButtonCallback(applicationWindow, genericCallback(MouseButtonCallback));
	glfwSetScrollCallback(applicationWindow, genericCallback(ScrollCallback));
}

void InputControl::bindErrorFunction() {
	glfwSetErrorCallback(errorCallbackStatic);
}

void InputControl::errorCallbackStatic(int err, const char* msg)
{
	printf("Error: %d %s\n", err, msg);
}
