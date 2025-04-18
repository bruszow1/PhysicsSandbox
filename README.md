The headless program can be compiled with the following:
    
    nvcc main.cu movement.cu bitmap.cpp shape.cpp utils.cu draw.cu -o assignment_adv_lib.exe -lcublas -lnppc -lnppial

The headless program creates a bitmap showing the initial state and a bitmap showing the state after a fixed number of iterations.

The GUI requires wxWidgets.
I used wx-3.2.6 in a directory of the same name.
My CMakeLists file is included, but it may not run properly on a different environment.

The GUI just refreshes the bitmap and almost all the code is drawn directly from the developer's site; I mainly included it for completion.

Right clicking advances by one frame and the enter key toggles between play and pause.