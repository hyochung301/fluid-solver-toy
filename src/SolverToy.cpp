#include "SolverToy.h"

SolverToy::SolverToy(size_t n) : N(n), userFFT(0) {}
SolverToy::~SolverToy() {}

void SolverToy::userLoop(float dt) {(void)dt;} // unimplemented by default

void SolverToy::use_fft(FFT_Solver2d* fft) {userFFT = fft;}
