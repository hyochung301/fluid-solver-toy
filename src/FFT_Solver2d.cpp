#include "FFT_Solver2d.h"

FFT_Solver2d::FFT_Solver2d(size_t n, float* buff) : N(n), buffer(buff) {}

FFT_Solver2d::~FFT_Solver2d() {}

FFTW_FFT_Solver2d::FFTW_FFT_Solver2d(size_t n, float* buff) : FFT_Solver2d(n,buff) {
    forw = fftwf_plan_dft_r2c_2d(N, N, buffer, (fftwf_complex*)buffer, FFTW_ESTIMATE);
    inv =     fftwf_plan_dft_c2r_2d(N, N, (fftwf_complex*)buffer, buffer, FFTW_ESTIMATE);
}

FFTW_FFT_Solver2d::~FFTW_FFT_Solver2d() {
    fftwf_destroy_plan(forw);
    fftwf_destroy_plan(inv);
}

void FFTW_FFT_Solver2d::cleanup() {
    fftwf_cleanup();
}

void FFTW_FFT_Solver2d::forward() {
	fftwf_execute(forw);
}

void FFTW_FFT_Solver2d::inverse() {
	fftwf_execute(inv);
}
