#ifndef FIELD_H 
#define FIELD_H 
#include <flgl/glm.h>
#include "StamFFT_FluidSolver.h"

struct Field {
	StamFFT_FluidSolver solver;

	Field(int const& N);

	glm::ivec2 mouse_pos()  const;
	glm::ivec2 mouse_delta() const;

	void step(float dt);

	float* buff() const;
};


#endif
