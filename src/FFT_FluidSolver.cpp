#include "FFT_FluidSolver.h"
#include <flgl/logger.h>
#include <random>
LOG_MODULE(fft_solver)

#define floor(x) ((x)>=0.0?((int)(x)):(-((int)(1-(x)))))
	
// === private members: ===
// const int n;
// fftwf_plan forward_u, forward_v, inv_u, inv_v;
// float* u, * v, * u0, * v0;
// float visc;
// float* buffer;

StamFFT_FluidSolver::StamFFT_FluidSolver(int const& N) : n(N), viscosity(visc) {
    alloc_buffers();
    initFFT(n, u0, v0);
}
StamFFT_FluidSolver::~StamFFT_FluidSolver() {
    destroyFFT();
    free_buffers();
}

void StamFFT_FluidSolver::initFFT(int const& N, float* u_buffer, float* v_buffer) {
    forward_u = fftwf_plan_dft_r2c_2d(n, n, u, (fftwf_complex*)u, FFTW_ESTIMATE);
    forward_v = fftwf_plan_dft_r2c_2d(n, n, v, (fftwf_complex*)v, FFTW_ESTIMATE);
    inv_u =     fftwf_plan_dft_c2r_2d(n, n, (fftwf_complex*)u, u, FFTW_ESTIMATE);
    inv_v =     fftwf_plan_dft_c2r_2d(n, n, (fftwf_complex*)v, v, FFTW_ESTIMATE);
}
void StamFFT_FluidSolver::alloc_buffers() {
    u = new float[n*n];
    v = new float[n*n];
    u0 = new float[n*(n+2)];
    v0 = new float[n*(n+2)];
    buffer = new float[n*n*2];
}
void StamFFT_FluidSolver::free_buffers() {
    delete [] u;
    delete [] v;
    delete [] u0;
    delete [] v0;
    delete [] buffer;
}
void StamFFT_FluidSolver::destroyFFT() {
    fftwf_destroy_plan(forward_u);
    fftwf_destroy_plan(forward_v);
    fftwf_destroy_plan(inv_u);
    fftwf_destroy_plan(inv_v);
    fftwf_cleanup();
    LOG_DBG("destroyed fftw");
}
void StamFFT_FluidSolver::stam_stable_solve(int const& N, 
                       float* const u, float* const v, 
                       float* const u0, float* const v0, 
                       float const& visc, float const& dt)
{
    float x, y, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    for ( i=0 ; i<N*N ; i++ ) {
        u[i] += dt*u0[i]; u0[i] = u[i];
        v[i] += dt*v0[i]; v0[i] = v[i];
    }

    for ( i=0 ; i<N ; i++ ) {
        for ( j=0 ; j<N ; j++ ) {
            x = i-dt*u0[i+N*j]*N; y = j-dt*v0[i+N*j]*N;
            i0 = floor(x); s = x-i0; i0 = (N+(i0%N))%N; i1 = (i0+1)%N;
            j0 = floor(y); t = y-j0; j0 = (N+(j0%N))%N; j1 = (j0+1)%N;
            u[i+N*j] = (1-s)*((1-t)*u0[i0+N*j0]+t*u0[i0+N*j1])+
                          s *((1-t)*u0[i1+N*j0]+t*u0[i1+N*j1]);
            v[i+N*j] = (1-s)*((1-t)*v0[i0+N*j0]+t*v0[i0+N*j1])+
                          s *((1-t)*v0[i1+N*j0]+t*v0[i1+N*j1]);
        }
    }

    for ( i=0 ; i<N ; i++ )
        for ( j=0 ; j<N ; j++ )
            { u0[i+(N+2)*j] = u[i+N*j]; v0[i+(N+2)*j] = v[i+N*j]; }

    fftwf_execute(forward_u); fftwf_execute(forward_v);

    for ( i=0 ; i<=N ; i+=2 ) {
        x = 0.5*i;
        for ( j=0 ; j<N ; j++ ) {
            y = j<=N/2 ? j : j-N;
            r = x*x+y*y;
            if ( r==0.0 ) continue;
            f = exp(-r*dt*visc);
            U[0] = u0[i  +(N+2)*j]; V[0] = v0[i  +(N+2)*j];
            U[1] = u0[i+1+(N+2)*j]; V[1] = v0[i+1+(N+2)*j];
            u0[i  +(N+2)*j] = f*( (1-x*x/r)*U[0]     -x*y/r *V[0] );
            u0[i+1+(N+2)*j] = f*( (1-x*x/r)*U[1]     -x*y/r *V[1] );
            v0[i+  (N+2)*j] = f*(   -y*x/r *U[0] + (1-y*y/r)*V[0] );
            v0[i+1+(N+2)*j] = f*(   -y*x/r *U[1] + (1-y*y/r)*V[1] );
        }
    }

    fftwf_execute(inv_u); fftwf_execute(inv_v);

    f = 1.0/(N*N);
    for ( i=0 ; i<N ; i++ )
        for ( j=0 ; j<N ; j++ )
        { 
                u[i+N*j] = f*u0[i+(N+2)*j]; v[i+N*j] = f*v0[i+(N+2)*j]; 
        }

    memset(u0,0,sizeof(float)*N*N);
    memset(v0,0,sizeof(float)*N*N);

}

void StamFFT_FluidSolver::fill_random(int mag) { 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mag, mag);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u0[i+j*n] = dis(gen); v0[i+j*n] = dis(gen);
        }
    }
}


/*
	Inputs:
		n: 			size of n by n field
		(u, v): 	velocity of the prev time step
		(u0, v0): 	a force field defined on a grid
		visc: 		viscosity of fluid
		dt: 		time step

		(u0, v0) must be of size n*(n+2) due to fft
*/
