#include <flgl.h>
#include <fftw3.h>
#include <cmath>
#include <flgl/tools.h>
#include <flgl/logger.h>
#include <Stopwatch.h>
#include <random>
#include "Stepper.h"
LOG_MODULE(main);
using namespace glm;

/*
	solver: https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf
	in this file I try to render their implementation
*/

static fftwf_plan forward_u, forward_v, inv_u, inv_v;
static fftwf_complex *in, *out;
#define n_g (256)

static void initFFT(int N, float* u, float* v) {
    forward_u = fftwf_plan_dft_r2c_2d(n_g, n_g, u, (fftwf_complex*)u, FFTW_ESTIMATE);
    forward_v = fftwf_plan_dft_r2c_2d(n_g, n_g, v, (fftwf_complex*)v, FFTW_ESTIMATE);
    inv_u = 	fftwf_plan_dft_c2r_2d(n_g, n_g, (fftwf_complex*)u, u, FFTW_ESTIMATE);
    inv_v = 	fftwf_plan_dft_c2r_2d(n_g, n_g, (fftwf_complex*)v, v, FFTW_ESTIMATE);
}

static void destroyFFT() {
	fftwf_destroy_plan(forward_u);
    fftwf_destroy_plan(forward_v);
    fftwf_destroy_plan(inv_u);
    fftwf_destroy_plan(inv_v);
    LOG_DBG("destroyed fftw");
}

#define floor(x) ((x)>=0.0?((int)(x)):(-((int)(1-(x)))))

void log_sum(int n, float* u, float* v) {
	float su = 0, sv = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			su += (u[i+j*n]); sv += (v[i+j*n]);
		}
	}
	LOG_DBG("\tsum: %E, %E",su,sv);
}
#define ls(msg, n, u, v) LOG_DBG(msg); log_sum(n, u, v);


/*
	Inputs:
		n: 			size of n by n field
		(u, v): 	velocity of the prev time step
		(u0, v0): 	a force field defined on a grid
		visc: 		viscosity of fluid
		dt: 		time step

		(u0, v0) must be of size n*(n+2) due to fft
*/
void stable_solve(int n, float* u, float* v, float* u0, float* v0, float const& visc, float const& dt)
{
	float x, y, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    for ( i=0 ; i<n*n ; i++ ) {
        u[i] += dt*u0[i]; u0[i] = u[i];
        v[i] += dt*v0[i]; v0[i] = v[i];
    }

    for ( i=0 ; i<n ; i++ ) {
        for ( j=0 ; j<n ; j++ ) {
            x = i-dt*u0[i+n*j]*n; y = j-dt*v0[i+n*j]*n;
            i0 = floor(x); s = x-i0; i0 = (n+(i0%n))%n; i1 = (i0+1)%n;
            j0 = floor(y); t = y-j0; j0 = (n+(j0%n))%n; j1 = (j0+1)%n;
            u[i+n*j] = (1-s)*((1-t)*u0[i0+n*j0]+t*u0[i0+n*j1])+
                          s *((1-t)*u0[i1+n*j0]+t*u0[i1+n*j1]);
            v[i+n*j] = (1-s)*((1-t)*v0[i0+n*j0]+t*v0[i0+n*j1])+
                          s *((1-t)*v0[i1+n*j0]+t*v0[i1+n*j1]);
        }
    }

    for ( i=0 ; i<n ; i++ )
        for ( j=0 ; j<n ; j++ )
            { u0[i+(n+2)*j] = u[i+n*j]; v0[i+(n+2)*j] = v[i+n*j]; }

    fftwf_execute(forward_u); fftwf_execute(forward_v);

    for ( i=0 ; i<=n ; i+=2 ) {
        x = 0.5*i;
        for ( j=0 ; j<n ; j++ ) {
            y = j<=n/2 ? j : j-n;
            r = x*x+y*y;
            if ( r==0.0 ) continue;
            f = exp(-r*dt*visc);
            U[0] = u0[i  +(n+2)*j]; V[0] = v0[i  +(n+2)*j];
            U[1] = u0[i+1+(n+2)*j]; V[1] = v0[i+1+(n+2)*j];
            u0[i  +(n+2)*j] = f*( (1-x*x/r)*U[0]     -x*y/r *V[0] );
            u0[i+1+(n+2)*j] = f*( (1-x*x/r)*U[1]     -x*y/r *V[1] );
            v0[i+  (n+2)*j] = f*(   -y*x/r *U[0] + (1-y*y/r)*V[0] );
            v0[i+1+(n+2)*j] = f*(   -y*x/r *U[1] + (1-y*y/r)*V[1] );
        }
    }

    fftwf_execute(inv_u); fftwf_execute(inv_v);

    f = 1.0/(n*n);
    for ( i=0 ; i<n ; i++ )
        for ( j=0 ; j<n ; j++ )
        { 
            	u[i+n*j] = f*u0[i+(n+2)*j]; v[i+n*j] = f*v0[i+(n+2)*j]; 
    	}

    memset(u0,0,sizeof(float)*n*n);
    memset(v0,0,sizeof(float)*n*n);

}

struct Field {
	const int n;
	float* u, * v, * u0, * v0;
	const float visc;
	float* buffer;
	Field(int const& N, float const& vis) : n(N), visc(vis) {
		u = new float[n*n];
		v = new float[n*n];
		u0 = new float[n*(n+2)];
		v0 = new float[n*(n+2)];
		buffer = new float[n*n*2];
		std::random_device rd;
	    std::mt19937 gen(rd());
	    std::uniform_real_distribution<> dis(-10000000., 10000000.);
		for (int i = 0; i < n_g; i++) {
			for (int j = 0; j < n_g; j++) {
				u0[i+j*n] = dis(gen); v0[i+j*n] = dis(gen);
			}
		}
	}
	~Field() {
		delete [] u;
		delete [] v;
		delete [] u0;
		delete [] v0;
		delete [] buffer;
	}

	ivec2 mouse_pos() {
		vec2 mp = window.mouse.pos;
		ivec2 const& f = window.frame;
		mp.y = ((float)f.y) - mp.y;
		return ivec2 {((mp.x / (float)f.x)) * n, ((mp.y / (float)f.y)) * n};
	}
	ivec2 mouse_delta() {
		vec2 md = window.mouse.delta;
		md.y *= -1.;
		ivec2 const& f = window.frame;
		ivec2 r = ivec2 {((md.x / (float)f.x)) * n, ((md.y / (float)f.y)) * n};
		return r;
	}

	void step(float dt) {
		if (window.mouse.left.down) {
			auto end   = mouse_pos();
			auto start = mouse_pos() - mouse_delta();
			auto delta = mouse_delta();

			for (auto pt : Stepper(start.x, start.y, end.x, end.y)) {
				u0[pt.x+pt.y*n] = delta.x * 20.;
				v0[pt.x+pt.y*n] = delta.y * 20.;
			}
		}

		stable_solve(n, u, v, u0, v0, visc, dt);
		float buffmax = -9999999.f;
		float buffmin = 9999999.f;
		for (int i = 0; i < n; i++) {
			int x = i - (n/2);
			for (int j = 0; j < n; j++) {
				int y = j - (n/2);
				int r = 2*(i+j*n);
				int g = 2*(i+j*n) + 1;
				buffer[r] = u[i+j*n];
				buffer[g] = v[i+j*n];
				if (buffer[r]>buffmax) buffmax = buffer[r];
				if (buffer[g]>buffmax) buffmax = buffer[g];
				if (buffer[r]<buffmin) buffmin = buffer[r];
				if (buffer[g]<buffmin) buffmin = buffer[g];
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int r = 2*(i+j*n);
				int g = 2*(i+j*n) + 1;
				float imag = (1.f / ((buffmax-buffmin)==0.?1.f:(buffmax-buffmin)));
				buffer[r] -= buffmin; buffer[r] *= imag;
				buffer[g] -= buffmin; buffer[g] *= imag;
			}
		}
	}

	float* buff() {return buffer;}
};

struct FieldRenderer {
	Mesh<Vt_classic> quad;
	Shader field_shad;
	Texture field_tex;
	const int n;
	Field field;

	FieldRenderer(int const& N) : n(N), field(N, 0.0001f) {}
	void init() {
		quad = DefaultMeshes::tile<Vt_classic>();
		field_shad = Shader::from_source("passthrough_vert", "tex");
		field_tex.create();
		field_tex.bind();
		field_tex.pixelate(true);
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

static FieldRenderer renderer(n_g);

class flog {
	const unsigned int N;
	float * buff;
	unsigned int i;
public:
	flog(unsigned int n) : N(n) {
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
		LOG_DBG("solver avg: %.1fus", a/N);
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
	initFFT(n_g, renderer.field.u0, renderer.field.v0);
	gl.init();
	window.create("fft", 1024, 1024);
	renderer.init();
	Stopwatch timer(SECONDS);
	renderer.buffer_texture();
	timer.start();	
	Stopwatch dtimer(SECONDS);
	dtimer.start();
	flog db(128);		
	while (!window.should_close()) {
		// if (window.keyboard[GLFW_KEY_SPACE].pressed) {
			auto st = timer.read(MICROSECONDS);
			renderer.field.step(dtimer.stop_reset_start());
			auto en = timer.read(MICROSECONDS);
			db.dump(en-st);
			renderer.buffer_texture();
		// }
		auto mp = renderer.field.mouse_pos();
		// LOG_DBG("mouse pos: %d,%d", mp.x,mp.y);
		renderer.render();
		if (window.keyboard[GLFW_KEY_ESCAPE].down) break;
		window.update();
	}

	gl.destroy();
	destroyFFT();
	return 0;
}









