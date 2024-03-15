#include <flgl.h>
#include <fftw3.h>
#include <cmath>
#include <flgl/tools.h>

/*
	solver: https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf
	in this file I try to render their implementation
*/

static fftw_plan plan_rc, plan_cr;
static fftw_complex *in, *out;
static int n;

static void initFFT(int N) {
	n = N;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n);

    plan_rc = fftw_plan_dft_2d(n, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_cr = fftw_plan_dft_2d(n, n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

static void destroyFFT() {
	fftw_destroy_plan(plan_rc);
    fftw_destroy_plan(plan_cr);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
}

#define FFT(s,u) \
if (s==1) fftw_execute(plan_rc); \
else fftw_execute(plan_cr);

#define floor(x) ((x)>=0.0?((int)(x)):(-((int)(1-(x)))))

/*
	Inputs:
		n: 			size of n by n field
		(u, v): 	velocity of the prev time step
		(u0, v0): 	a force field defined on a grid
		visc: 		viscosity of fluid
		dt: 		time step

		(u0, v0) must be of size n*(n+2) due to fft
*/
void stable_solve(int n, float* u, float* v, float* u0, float* v0, float visc, float dt)
{
	float x, y, x0, y0, f, r, U[2], V[2], s, t;
	int i, j, i0, j0, i1, j1;

	// add force fields to velo fields
	// then make force fields equal velo field
	for (i = 0; i < n*n; i++) {
		u[i] += dt*u0[i]; u0[i] = u[i];
		v[i] += dt*v0[i]; v0[i] = v[i];
	}

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
			u[i+n*j] = (1-s) * ((1-t) * u0[i0+n*j0] + t*u0[i0+n*j1]) + 
						  s  * ((1-t) * u0[i1+n*j0] + t*u0[i1+n*j1]);
			v[i+n*j] = (1-s) * ((1-t) * v0[i0+n*j0] + t*v0[i0+n*j1]) + 
						  s  * ((1-t) * v0[i1+n*j0] + t*v0[i1+n*j1]);
		}
	}

	// copy interpolated values back to fft buffers
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			u0[i+(n+2)*j] = u[i+n*j];
			v0[i+(n+2)*j] = v[i+n*j];
		}
	}

	FFT(1,u0); FFT(1,v0); // transform velo field

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

	FFT(-1, u0); FFT(-1,v0); // inverse transform

	// normalize 
	f = 1.0/(n*n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			u[i+n*j] = f*u0[i+(n+2)*j];
			v[i+n*j] = f*v0[i+(n+2)*j];
		}
	}
}

struct FieldRenderer {
	Mesh<Vt_classic> quad;
	Shader field_shad;
	Texture field_tex;
	const int n;

	FieldRenderer(int const& N) : n(N) {}
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
    	field_tex.alloc(GL_RGB32F, n, n, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
		field_tex.unbind();
		// LOG_DBG("updating took %3f ms", stop);
	}
	void render() {
		gl.clear();
		field_tex.bind();

		field_shad.bind();
		gl.draw_mesh(quad);
		field_tex.unbind();
		field_shad.unbind();
	}
};

static FieldRenderer renderer(n);

int main() {
	initFFT(512);
	destroyFFT();
	return 0;
}









