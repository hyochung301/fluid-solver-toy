#include <flgl.h>
#include <flgl/tools.h>
#include <flgl/logger.h>
#include <Stopwatch.h>
#include <random>
#include "Stepper.h"
#include "FFT_FluidSolver.h"
LOG_MODULE(main);
using namespace glm;

struct Field {
	StamFFT_FluidSolver solver;

	Field(int const& N, float const& vis) : solver(N) {
		solver.random_fill(100.);
	}
	~Field() {
	}

	ivec2 mouse_pos() {
		vec2 mp = window.mouse.pos;
		ivec2 const& f = window.frame;
		mp.y = ((float)f.y) - mp.y;
		return ivec2 {((mp.x / (float)f.x)) * solver.dim(), ((mp.y / (float)f.y)) * solver.dim()};
	}
	ivec2 mouse_delta() {
		vec2 md = window.mouse.delta;
		md.y *= -1.;
		ivec2 const& f = window.frame;
		ivec2 r = ivec2 {((md.x / (float)f.x)) * solver.dim(), ((md.y / (float)f.y)) * solver.dim()};
		return r;
	}

	void step(float dt) {
		if (window.mouse.left.down) {
			auto end   = mouse_pos();
			auto start = mouse_pos() - mouse_delta();
			auto delta = mouse_delta();

			for (auto pt : Stepper(start.x, start.y, end.x, end.y)) {
				solver.set_force(pt.x, pt.y, delta.x * 33., delta.y * 33.);
			}
		}
		solver.step(dt);
	}

	float* buff() {return solver.buff();}
};

struct FieldRenderer {
	Mesh<Vt_classic> quad;
	Shader field_shad;
	Texture field_tex;
	const int n;
	Field field;

	FieldRenderer(const int& N) : field(N, 0.0001f), n(N) {}
	void init() {
		quad = DefaultMeshes::tile<Vt_classic>();
		field_shad = Shader::from_source("passthrough_vert", "tex");
		field_tex.create();
		field_tex.bind();
		field_tex.pixelate(false);
    	field_tex.unbind();
	}
	void buffer_texture() {
		field_tex.bind();
    	field_tex.alloc(GL_RG32F, n, n, GL_RG, GL_FLOAT, field.buff());
		field_tex.unbind();
	}
	void render() {
		gl.clear();
		field_tex.bind();
		field_tex.bind_to_unit(0);
		field_shad.bind();
		field_shad.uInt("uTexslot", 0);
		gl.draw_mesh(quad);
		field_tex.unbind();
		field_shad.unbind();
	}
};

static FieldRenderer renderer(256);

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
	renderer.field.solver.random_fill(400.);
	bool slowmo = false;
	while (!window.should_close()) {

		if (window.keyboard[GLFW_KEY_R].down)
			renderer.field.solver.random_fill(400.);
		if (window.keyboard[GLFW_KEY_S].pressed)
			slowmo = !slowmo;
		if (window.keyboard[GLFW_KEY_0].pressed)
			renderer.field.solver.zero_field();
		if (window.keyboard[GLFW_KEY_UP].pressed) {
			renderer.field.solver.viscosity *= 2.; LOG_DBG("viscosity: %e", renderer.field.solver.viscosity);
		}
		if (window.keyboard[GLFW_KEY_DOWN].pressed) {
			renderer.field.solver.viscosity /= 2.; LOG_DBG("viscosity: %e", renderer.field.solver.viscosity);
		}

		auto st = timer.read(MICROSECONDS);
		renderer.field.step(slowmo ? dtimer.stop_reset_start()/10. : dtimer.stop_reset_start());
		renderer.field.solver.slow_fill_pixbuff();
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









