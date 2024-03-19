#ifndef FFT_FLUID_SOLVER_H
#define FFT_FLUID_SOLVER_H
#include <fftw3.h>
#include <cmath>
#include "../lib/sw/Stopwatch.h"

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
    float force_mul;

    Stopwatch timer;
    float t_us_solver;
    float t_us_ffts;

    void initFFT(int const& N, float* u_buffer, float* v_buffer);
    void alloc_buffers();
    void free_buffers();
    void destroyFFT();
    void stam_stable_solve(int const& N, 
                           float* const u, float* const v, 
                           float* const u0, float* const v0, 
                           float const& visc, float const& dt);
public:
    float & viscosity;
    int const& dim() const;

    StamFFT_FluidSolver(int const& N);
    ~StamFFT_FluidSolver();

    void random_fill(float mag);
    void zero_field();

    float* buff() const;
    float* x_buffer() const;
    float* y_buffer() const;

    void add_force(int x, int y, int fx, int fy);
    void set_force(int x, int y, int fx, int fy);
    void get_force(int x, int y, int& fx, int& fy) const;
    void set_force_multiplier(float mul);
    float force_multiplier();

    void step(float const& dt);
    float get_prev_solver_t();
    float get_prev_fft_t();

    void slow_fill_pixbuff();

};

#endif
