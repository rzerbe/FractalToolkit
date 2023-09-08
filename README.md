# FractalToolkit
An interactive fractal visualizer. 

Supports Mandelbrot, with plans to expand to Buddabrot and other parameterized M-set equations. Displays to screen with OpenGL and takes KBM input using GLEW library. Mandelbrot render window can be moved by left click and zoomed in/out with mouse wheel. Commands to manually move the window center and save output image are available in command window.

Calculations have multithreading support can be switched between regular C, AVX, and CUDA instructions. Smooth shading, which takes the magnitude component per pixel, is implemented to provide gradual but accurate transitions between escape iteration regions. Supersampling Antialiasing multiplication factors can be set to improve render appearance near "high noise" areas. This operation requires a large amount of memory. For CPU based calculation, this is limited by your RAM. For GPU calculation (CUDA), this is limited by your VRAM.

Colors can be added in [R, G, B, A] format or removed from the palette.json file. 

Program settings are stored in properties.json.</br>
General</br>
  &emsp;Width - X resolution of window.</br>
  &emsp;Height - Y resolution of window.</br>
  &emsp;Threads - Maximum number of threads. It is recommended to set this to the number of physical cores.</br>
  &emsp;Instruction_Set - Sets the instruction set to use when calculating.</br>
    &emsp;&emsp;0 = C_64, </br>
    &emsp;&emsp;1 = AVX_64, </br>
    &emsp;&emsp;2 = AVX_32, </br>
    &emsp;&emsp;3 = CUDA</br>

Mandelbrot</br>
  &emsp;Window.X.Min/Max - Starting horizontal view port</br>
  &emsp;Window.Y.Min/Max - Starting vertical view port</br>
  &emsp;Point.Cr - Starting center real component</br>
  &emsp;Point.Ci - Starting center imaginary component</br>
  &emsp;Point.Rr - Magnification factor</br>
  &emsp;Max_Iterations - Limit to stop calculation and consider point as unescaped in escape time algorithm.</br>
  &emsp;Shading_Mode - Experimental monochromatic shader.</br>
  &emsp;Smooth_Shading - Uses magnitude component as shading factor to pick interpolated colors between two colors.</br>
  &emsp;Super_Sampling_Anti_Aliasing - Calculates a larger window space multiplied by this factor. For each pixel, a 2x factor will result in a 2x2 subpixel, 4x creates a 4x4 subpixel. The result is downsampled by the factor.</br>
  &emsp;Filename - Default name to use when using 'save' command</br>
  &emsp;Palette - The name of the palette file.</br>

Smooth Shading ON
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/98fc7384-6f2a-43e9-89f1-672d07026c4e)

Smooth Shading OFF
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/bbfb6a71-8a90-4c60-804a-a54fb1d6798c)

Smooth Shading ON, Supersampling Antialiasing = 1x
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/f9463b62-baf4-4693-ab2b-8c158c0df273)

Smooth Shading ON, Supersampling Antialiasing = 8x
![image](https://github.com/rzerbe/FractalToolkit/assets/14305489/5785bc63-c4c0-49f0-98dc-c4ebbacf2e49)
