#ifndef FFT_FLUID_SOLVER_H
#define FFT_FLUID_SOLVER_H
#include <fftw3.h>
#include <cmath>

/*
    solver: https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf
*/

class StamFFT_FluidSolver {
private:
    const int n;
    fftwf_plan forward_u, forward_v, inv_u, inv_v;
    float* u, * v, * u0, * v0;
    float visc;
    float* buffer;

    void initFFT(int const& N, float* u_buffer, float* v_buffer);
    void alloc_buffers();
    void free_buffers();
    void destroyFFT();
    void stam_stable_solve(int const& N, 
                           float* const u, float* const v, 
                           float* const u0, float* const v0, 
                           float const& visc, float const& dt);
public:
    float const& viscosity;
    int const& dim();

    StamFFT_FluidSolver(int const& N);
    ~StamFFT_FluidSolver();

    void random_fill(float mag);

    float* buff();
    float* x_buffer();
    float* y_buffer();

    void add_force(int x, int y, int fx, int fy);
    void set_force(int x, int y, int fx, int fy);
    void get_force(int x, int y, int& fx, int& fy);

    void step(float const& dt);

    void slow_fill_pixbuff();

};

#endif
