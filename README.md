# FractalToolkit
An interactive fractal visualizer. 

Supports Mandelbrot, with plans to expand to Buddabrot and other parameterized M-set equations. Displays to screen with OpenGL and takes KBM input using GLEW library. Mandelbrot render window can be moved by left click and zoomed in/out with mouse wheel. Commands to manually move the window center and save output image are available in command window.

Calculations have multithreading support can be switched between regular C, AVX, and CUDA instructions. Smooth shading, which takes the magnitude component per pixel, is implemented to provide gradual but accurate transitions between escape iteration regions. Supersampling Antialiasing multiplication factors can be set to improve render appearance near "high noise" areas. This operation requires a large amount of memory. For CPU based calculation, this is limited by your RAM. For GPU calculation (CUDA), this is limited by your VRAM.

Smooth Shading ON
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/98fc7384-6f2a-43e9-89f1-672d07026c4e)

Smooth Shading OFF
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/bbfb6a71-8a90-4c60-804a-a54fb1d6798c)

Smooth Shading ON, Supersampling Antialiasing = 1x
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/f9463b62-baf4-4693-ab2b-8c158c0df273)

Smooth Shading ON, Supersampling Antialiasing = 8x
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/5785bc63-c4c0-49f0-98dc-c4ebbacf2e49)
