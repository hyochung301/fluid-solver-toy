#include <flgl.h>
#include <fftw3.h>
#include <cmath>
#include <flgl/tools.h>
#include <flgl/logger.h>
#include <Stopwatch.h>
#include <random>
LOG_MODULE(main);

/*
	solver: https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf
	in this file I try to render their implementation
*/

static fftwf_plan forward_u, forward_v, inv_u, inv_v;
static fftwf_complex *in, *out;
#define n_g (128)

static void initFFT(int N, float* u, float* v) {
    forward_u = fftwf_plan_dft_r2c_2d(n_g, n_g, u, (fftwf_complex*)u, FFTW_ESTIMATE);
    forward_v = fftwf_plan_dft_r2c_2d(n_g, n_g, v, (fftwf_complex*)v, FFTW_ESTIMATE);
    inv_u = 	fftwf_plan_dft_c2r_2d(n_g, n_g, (fftwf_complex*)v, v, FFTW_ESTIMATE);
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
	float x, y, x0, y0, f, r, U[2], V[2], s, t;
	int i, j, i0, j0, i1, j1;

	ls("force", n, u0, v0);

	// add force fields to velo fields
	// then make force fields equal velo field
	for (i = 0; i < n*n; i++) {
		u[i] += dt*u0[i]; u0[i] = u[i];
		v[i] += dt*v0[i]; v0[i] = v[i];
	}

	ls("after force", n, u, v);

	/*
		... "This effect is called “self-advection”: an auto-feedback mechanism
		within the fluid. To solve for this effect we use a very elegant technique 
		known as a “semi-Lagrangian” scheme in the computational fluid dynamics literature. 
		Here is how it works. For each voxel of the velocity grid we trace its midpoint backwards 
		through the velocity field over a time step dt. This point will end up somewhere else 
		in the grid. We then linearly interpolate the velocity at that point from the 
		neighboring voxels and transfer this velocity back to the departure voxel. 
		This requires two grids, one for the interpolation (u0,v0) and one to store 
		the new interpolated values (u,v). To perform the particle trace, we use a 
		simple linear approximation of the exact path. 
		The interpolation is very easy to implement because our grids are periodic. 
		Points that leave the grid simply reenter the grid from the opposite side."
	*/
	for (x=0.5/n, i=0; i < n; i++, x+=1.0/n) {
		for (y=0.5/n, j=0; j < n; j++, y+=1.0/n) {
			x0 = n*(x-dt*u0[i+n*j]) - 0.5;
			y0 = n*(y-dt*v0[i+n*j]) - 0.5;
			i0 = floor(x0); s = x0-i0; i0 = (n+(i0%n))%n; i1 = (i0+1)%n;
			j0 = floor(y0); t = y0-j0; j0 = (n+(j0%n))%n; j1 = (j0+1)%n;
			u[i+n*j] = ((1-s) * ( ((1-t) * u0[i0+n*j0]) + (t*u0[i0+n*j1]) )) + 
						  (s  * ( ((1-t) * u0[i1+n*j0]) + (t*u0[i1+n*j1]) ));
			v[i+n*j] = ((1-s) * ( ((1-t) * v0[i0+n*j0]) + (t*v0[i0+n*j1]) )) + 
						  (s  * ( ((1-t) * v0[i1+n*j0]) + (t*v0[i1+n*j1]) ));
			// if ((i == 6) || u[i+n*j] != 0.0 || v[i+n*j] != 0.0) 
			// 	LOG_DBG("%f", t);
			// u[i+n*j] = u0[i+n*j];
			// v[i+n*j] = v0[i+n*j];
		}
	}

	ls("after self advec", n, u, v);

	// copy interpolated values back to fft buffers
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			u0[i+(n+2)*j] = u[i+n*j];
			v0[i+(n+2)*j] = v[i+n*j];
		}
	}

	ls("0s pre fft", n, u0, v0);

	fftwf_execute(forward_u); fftwf_execute(forward_v); // transform velo field

	ls("0s post fft", n, u0, v0);

	// apply visc filter and make mass conserving. 
	// can make mask conserving by projecting fourier velo field
	// onto line perpendicular to wave num
	for (i = 0; i <= n; i += 2) {
		x = 0.5*i;
		for (j = 0; j < n; j++) {
			// calculate x,y,r relative to origin
			y = j<=n/2 ? j : j-n;
			r = x*x+y*y;
			if (r == 0.0) continue;
			f = exp(-r*dt*visc);	// this filter converges to true fluid eqs
			// U[0,1] is R,I of x coord, V[0,1] is R,I of y coord
			U[0] = u0[i  +(n+2)*j]; V[0] = v0[i  +(n+2)*j];
			U[1] = u0[i+1+(n+2)*j]; V[1] = v0[i+1+(n+2)*j];
			// project and apply filter
			u0[i  +(n+2)*j] = f*( (1-x*x/r)*U[0]      -x*y/r *V[0] );
			u0[i+1+(n+2)*j] = f*( (1-x*x/r)*U[1]      -x*y/r *V[1] );
			v0[i  +(n+2)*j] = f*(    -y*x/r*U[0]  + (1-y*y/r)*V[0] );
			v0[i+1+(n+2)*j] = f*(    -y*x/r*U[1]  + (1-y*y/r)*V[1] );
		}
	}

	ls("0s post filter", n, u0, v0);

	fftwf_execute(inv_u); fftwf_execute(inv_v); // inverse transform

	ls("0s post inverse", n, u0, v0);

	// normalize 
	f = 1.0/(n*n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			u[i+n*j] = f*u0[i+(n+2)*j];
			v[i+n*j] = f*v0[i+(n+2)*j];
		}
	}

	memset(u0, 0x00, n*(n+2)*sizeof(float));
	memset(v0, 0x00, n*(n+2)*sizeof(float));

	ls("returning", n, u, v);
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
	    std::uniform_real_distribution<> dis(-2., 2.);
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

	void step(float dt) {
		LOG_DBG("dt=%E",dt);
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
				// buffer[r] = x*x+y*y>50*50?1.:0.;
				// buffer[g] = x*x+y*y>50*50?1.:0.;
				if (buffer[r]>buffmax) buffmax = buffer[r];
				if (buffer[g]>buffmax) buffmax = buffer[g];
				if (buffer[r]<buffmin) buffmin = buffer[r];
				if (buffer[g]<buffmin) buffmin = buffer[g];
			}
		}
		LOG_DBG("buffmax: %E",buffmax);
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

	FieldRenderer(int const& N) : n(N), field(N, 0.001f) {}
	void init() {
		quad = DefaultMeshes::tile<Vt_classic>();
		field_shad = Shader::from_source("passthrough_vert", "tex");
		field_tex.create();
		field_tex.bind();
		field_tex.pixelate();
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

int main() {
	initFFT(n_g, renderer.field.u0, renderer.field.v0);
	gl.init();
	window.create("fft", 1024, 1024);
	renderer.init();
	Stopwatch timer(SECONDS);
	renderer.buffer_texture();
	timer.start();			
	while (!window.should_close()) {
		if (window.keyboard[GLFW_KEY_SPACE].pressed) {
			renderer.field.step(1.f);
			renderer.buffer_texture();
		}
		renderer.render();
		if (window.keyboard[GLFW_KEY_ESCAPE].down) break;
		window.update();
	}

	gl.destroy();
	destroyFFT();
	return 0;
}









