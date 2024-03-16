#include "Field.h"
#include "Stepper.h"
#include <flgl.h>
using namespace glm;

Field::Field(int const& N) : solver(N) {
	solver.random_fill(100.);
}

ivec2 Field::mouse_pos() const {
	vec2 mp = window.mouse.pos;
	ivec2 const& f = window.frame;
	mp.y = ((float)f.y) - mp.y;
	return ivec2 {((mp.x / (float)f.x)) * solver.dim(), ((mp.y / (float)f.y)) * solver.dim()};
}
ivec2 Field::mouse_delta() const {
	vec2 md = window.mouse.delta;
	md.y *= -1.;
	ivec2 const& f = window.frame;
	ivec2 r = ivec2 {((md.x / (float)f.x)) * solver.dim(), ((md.y / (float)f.y)) * solver.dim()};
	return r;
}

void Field::step(float dt) {
	if (window.mouse.left.down) {
		auto end   = mouse_pos();
		auto start = mouse_pos() - mouse_delta();
		auto delta = mouse_delta();

		for (auto pt : Stepper(start.x, start.y, end.x, end.y)) {
			solver.set_force(pt.x  , pt.y  , delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x+1, pt.y  , delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x  , pt.y+1, delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x-1, pt.y  , delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x  , pt.y-1, delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x+1, pt.y+1, delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x+1, pt.y-1, delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x-1, pt.y+1, delta.x * 33., delta.y * 33.);
			solver.set_force(pt.x-1, pt.y-1, delta.x * 33., delta.y * 33.);
		}
	}
	solver.step(dt);
}

float* Field::buff() const {return solver.buff();}
