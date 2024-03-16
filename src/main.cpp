#include <flgl.h>
#include <flgl/tools.h>
#include <flgl/logger.h>
#include <Stopwatch.h>
#include <random>
#include "Stepper.h"
#include "StamFFT_FluidSolver.h"
#include "Field.h"
#include "FieldRenderer.h"
LOG_MODULE(main);
using namespace glm;

#define DIM (256)

static Field field(DIM);
static FieldRenderer renderer(field);

class flog {
	std::string name;
	const unsigned int N;
	float * buff;
	unsigned int i;
public:
	flog(unsigned int n, std::string && nam) : N(n), name(nam) {
		buff = new float[N]; i = 0;
	}
	~flog() {
		delete [] buff;
	}
	void outavg() {
		float a = 0.;
		for (int j = 0; j < N; j++) {
			a += buff[j];
		}
		LOG_DBG(name + " rolling avg time: %.1fus", a/N);
	}
	void dump(float v) {
		buff[i++] = v;
		if (!(i < N)) {
			i = 0;
			outavg();
		}
	}
};

int main() {
	gl.init();
	window.create("fft", 1024, 1024);
	renderer.init();
	Stopwatch timer(SECONDS);
	renderer.buffer_texture();
	timer.start();	
	Stopwatch dtimer(SECONDS);
	dtimer.start();
	flog sv(256, "solver");	
	flog rd(256, "render");		
	field.solver.random_fill(400.);
	bool slowmo = false;
	while (!window.should_close()) {

		if (window.keyboard[GLFW_KEY_R].down)
			field.solver.random_fill(400.);
		if (window.keyboard[GLFW_KEY_S].pressed)
			slowmo = !slowmo;
		if (window.keyboard[GLFW_KEY_0].pressed)
			field.solver.zero_field();
		if (window.keyboard[GLFW_KEY_UP].pressed) {
			field.solver.viscosity *= 2.; LOG_DBG("viscosity: %e", field.solver.viscosity);
		}
		if (window.keyboard[GLFW_KEY_DOWN].pressed) {
			field.solver.viscosity /= 2.; LOG_DBG("viscosity: %e", field.solver.viscosity);
		}

		auto st = timer.read(MICROSECONDS);
		field.step(slowmo ? dtimer.stop_reset_start()/10. : dtimer.stop_reset_start());
		field.solver.slow_fill_pixbuff();
		auto en = timer.read(MICROSECONDS);
		sv.dump(en-st);

		st = timer.read(MICROSECONDS);
		renderer.buffer_texture();
		renderer.render();
		en = timer.read(MICROSECONDS);
		rd.dump(en-st);

		if (window.keyboard[GLFW_KEY_ESCAPE].down) break;
		window.update();
	}

	gl.destroy();
	return 0;
}









