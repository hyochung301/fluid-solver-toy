#include <flgl.h>
#include <flgl/tools.h>
#include <flgl/logger.h>
#include <Stopwatch.h>
#include <random>
#include "Stepper.h"
#include "StamFFT_FluidSolver.h"
#include "Field.h"
#include "FieldRenderer.h"
#include "rgFieldRenderer.h"
#include "vecFieldRenderer.h"
LOG_MODULE(main);
using namespace glm;

#define DIM (256)

static Field field(DIM);
static FieldRenderer* rendererA = 0, * rendererB = 0;

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

static void explode(int x0, int y0, int r, float force=20000.) {
	const float mag = 1000.;
	for (int y = y0-r; y < y0+r; y++) {
		int dy = y-y0;
		for (int x = x0-r; x < x0+r; x++) {
			int dx = x-x0;
			if (dx*dx+dy*dy <= r*r) {
				int dx = x-x0; int dy = y-y0;
				float fx = dx * (force/(dx*dx+dy*dy)) * field.solver.force_multiplier(); 
				float fy = dy * (force/(dx*dx+dy*dy)) * field.solver.force_multiplier();
				field.solver.add_force(x,y, fx,fy);
			}
		}
	}
}

int main() {
	gl.init();
	window.create("fluid solver toy", 1024, 1024);
	Stopwatch timer(SECONDS); timer.start();
	Stopwatch dtimer(SECONDS); dtimer.start();

	rendererA = new rgFieldRenderer(field);
	rendererB = new vecFieldRenderer(field);
	field.solver.random_fill(400.);

	rendererA->init();
	rendererB->init();

	rendererA->prepare();
	rendererB->prepare();

	flog sv(256, "solver");	
	flog rd(256, "render");		
	bool slowmo = false;
	int view = 0;
	while (!window.should_close()) {
			// controls
		if (window.keyboard[GLFW_KEY_ESCAPE].down) break;
		if (window.keyboard[GLFW_KEY_R].down)
			field.solver.random_fill(400.);
		if (window.keyboard[GLFW_KEY_S].pressed) {
			slowmo = !slowmo;
			field.solver.set_force_multiplier(1. + ((float)(slowmo*9)));
		}
		if (window.keyboard[GLFW_KEY_0].pressed)
			field.solver.zero_field();
		if (window.keyboard[GLFW_KEY_UP].pressed) {
			field.solver.viscosity *= 2.; LOG_DBG("viscosity: %e", field.solver.viscosity);
		}
		if (window.keyboard[GLFW_KEY_DOWN].pressed) {
			field.solver.viscosity /= 2.; LOG_DBG("viscosity: %e", field.solver.viscosity);
		}
		if (window.mouse.left.pressed && window.keyboard[GLFW_KEY_SPACE].down)
			explode(field.mouse_pos().x, field.mouse_pos().y, 128);
		if (window.keyboard[GLFW_KEY_V].pressed)
			view = (view+1)%3;

			// solver step
		auto st = timer.read(MICROSECONDS);
		field.step(slowmo ? dtimer.stop_reset_start()/10. : dtimer.stop_reset_start());
		field.solver.slow_fill_pixbuff();
		auto en = timer.read(MICROSECONDS);
		sv.dump(en-st);

			// render
		st = timer.read(MICROSECONDS);
		gl.clear();
		switch(view) {
		case 2:
			rendererA->prepare();
			rendererA->render();
		case 1:
			rendererB->prepare();
			rendererB->render();
			break;
		case 0:
			rendererA->prepare();
			rendererA->render();
			break;
		}
		en = timer.read(MICROSECONDS);
		rd.dump(en-st);

		window.update();
	}

	rendererA->destroy();
	delete rendererA;
	rendererB->destroy();
	delete rendererB;
	gl.destroy();
	return 0;
}









