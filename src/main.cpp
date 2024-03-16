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
LOG_MODULE(main);
using namespace glm;

#define DIM (256)

static Field field(DIM);
static rgFieldRenderer renderer(field);

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

static void explode(int x0, int y0, int r) {
	const float mag = 1000.;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mag, mag);
	for (int y = y0-r; y < y0+r; y++) {
		int dy = y-y0;
		for (int x = x0-r; x < x0+r; x++) {
			int dx = x-x0;
			if (dx*dx+dy*dy <= r*r) {
				int dx = x-x0; int dy = y-y0;
				float fx = dx * (1000./(dx*dx+dy*dy)); 
				float fy = dy * (1000./(dx*dx+dy*dy));
				LOG_DBG("%.2f,%.2f",fx+dis(gen),fy+dis(gen));
				field.solver.add_force(x,y, fx,fy);
			}
		}
	}
}

int main() {
	gl.init();
	window.create("fft", 1024, 1024);
	renderer.init();
	Stopwatch timer(SECONDS);
	renderer.prepare();
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
		if (window.mouse.left.pressed && window.keyboard[GLFW_KEY_SPACE].down)
			explode(field.mouse_pos().x, field.mouse_pos().y, 128);

		auto st = timer.read(MICROSECONDS);
		field.step(slowmo ? dtimer.stop_reset_start()/10. : dtimer.stop_reset_start());
		field.solver.slow_fill_pixbuff();
		auto en = timer.read(MICROSECONDS);
		sv.dump(en-st);

		st = timer.read(MICROSECONDS);
		renderer.prepare();
		renderer.render();
		en = timer.read(MICROSECONDS);
		rd.dump(en-st);

		if (window.keyboard[GLFW_KEY_ESCAPE].down) break;
		window.update();
	}

	gl.destroy();
	return 0;
}









