# fft fluid solver gfx demo

Click and drag to swirl fluid around.
     
I've been wanting to learn how these kinds of fluid sims work and of course the rabbit hole is deep. (computational fluid dynamics is complicated -- really?)      
I found [this paper by Jos Stam (2001)](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf) that details a "simple" and "good for beginners" FFT-based technique. 

Here I've implemented it with [MIT's FFTW FFT library](https://www.fftw.org/) and rendered the fluid field. Currently the renderer is very dumb: simply renders the x and y components of the vector field in the red and green components of the color at a given pixel. It's not a good representation of the field but is enough that you can easily see the fluid behavior.

# Building
Like any of my graphics projects, this will build on mac windows or linux via cmake. More on that in the [flgl readme](https://github.com/collebrusco/flgl).   
**However** this code depends on fftw as well. Place libfftw3.a built for your system in lib/fftw/bin as I have done with the library built for mac M1 arm64. Be sure to include the single precision FFTs.

