#ifndef SOLVER_TOY
#define SOLVER_TOY

#include <stdint.h>
#include "FFT_Solver2d.h"

// Entire fluid sim & gfx demo wrapped up into an object
// usage:
// SolverToy toy(128); toy.useFFT(yourFFTsolver); toy.run();
// extend if you want to add code to each loop with userLoop

class SolverToy {
private:
	// called every frame, override if you need/want extra functionality
	virtual void userLoop(float dt);
protected:
	const size_t N;
	FFT_Solver2d* userFFT;
public:
	// solver will run n x n simulation, MUST BE POW OF TWO
	SolverToy(size_t n);
	virtual ~SolverToy();

	// specify what FFT solver the fluid sim will use
	void use_fft(FFT_Solver2d* fft);

	// enter main loop, runs until esc pressed or win closed
	void run();

};

#endif
