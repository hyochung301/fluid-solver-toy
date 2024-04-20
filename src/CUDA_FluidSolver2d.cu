#include "CUDA_FluidSolver2d.cuh"
#include <flgl/logger.h>
#include <random>
#include <cstring>
#include <cuda.h>
LOG_MODULE(fluid_solver)

#define cu() \
{cudaError_t err = cudaGetLastError();en++;\
while (err != cudaSuccess) {\
    printf("E %d: %s\n", en, cudaGetErrorString(err));\
    err = cudaGetLastError(); \
}}\

// === private members: ===
// const int n;
// fftwf_plan forward_u, forward_v, inv_u, inv_v;
// float* u, * v, * u0, * v0;
// float visc;
// float* buffer;


__global__ void apply_forces_kernel(float* u, float* v, float* u0, float* v0, float dt, int N) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N*N) {
		u[i] += dt*u0[i]; 
		u0[i] = u[i];
		v[i] += dt*v0[i]; 
		v0[i] = v[i];
	}
}

__global__ void self_advection_kernel(float* u, float* v, float* u0, float* v0, float dt, int N) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int i = index / N;
	int j = index % N;

	if (index < N*N) {
		float x = i - dt * u0[i + N*j] * N;
		float y = j - dt * v0[i + N*j] * N;

		int i0 = floor(x); float s = x - i0; i0 = (N + (i0 % N)) % N; int i1 = (i0 + 1) % N;
		int j0 = floor(y); float t = y - j0; j0 = (N + (j0 % N)) % N; int j1 = (j0 + 1) % N;

		u[i + N*j] = (1 - s) * ((1 - t) * u0[i0 + N*j0] + t * u0[i0 + N*j1]) +
					  s * ((1 - t) * u0[i1 + N*j0] + t * u0[i1 + N*j1]);
		v[i + N*j] = (1 - s) * ((1 - t) * v0[i0 + N*j0] + t * v0[i0 + N*j1]) +
					  s * ((1 - t) * v0[i1 + N*j0] + t * v0[i1 + N*j1]);
	}
}

__global__ void copy_to_fourier_buffers_kernel(float* u, float* v, float* u0, float* v0, int N) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int i = index / N;
	int j = index % N;

	if (index < N*N) {
		u0[2*(i + N*j)] = u[i+N*j]; // Real part of u0
		v0[2*(i + N*j)] = v[i+N*j]; // Real part of v0
		u0[2*(i + 1 + N*j)] = 0.; // Imaginary part of u0
		v0[2*(i + 1 + N*j)] = 0.; // Imaginary part of v0
	}
}

__global__ void applyLowPassFilter(float* u0, float* v0, float dt, float visc, int N) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int i = index / N;
	int j = index % N;

	if (index < N * N) {
		float x = (i <= N/2) ? i : (float)i - (float)N;
		float y = (j <= N/2) ? j : (float)j - (float)N;
		float r = x*x + y*y;
		if (r != 0.0) {
			float *uf = &(u0[2*(i + N*j)]);
			float *vf = &(v0[2*(i + N*j)]);

			float f = exp(-r * dt * visc);

			float ur = f * ( (1 - x*x/r)*uf[0] - x*y/r * vf[0] );
			float ui = f * ( (1 - x*x/r)*uf[1] - x*y/r * vf[1] );
			float vr = f * ( -y*x/r * uf[0] + (1 - y*y/r)*vf[0] );
			float vi = f * ( -y*x/r * uf[1] + (1 - y*y/r)*vf[1] );

			uf[0] = ur;
			uf[1] = ui;
			vf[0] = vr;
			vf[1] = vi;
		}
	}
}

__global__ void normalizeAndClear(float* u, float* v, float* u0, float* v0, float f, int N) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int i = index / N;
	int j = index % N;

	if (index < N * N) {
		u[i+N*j] = f * u0[2*(i + N*j)]; // Real part of u0
		v[i+N*j] = f * v0[2*(i + N*j)]; // Real part of v0

		u0[2*(i + N*j)] = 0.; // Clear u0
		v0[2*(i + N*j)] = 0.; // Clear v0
	}
}

StamFFT_FluidSolver::StamFFT_FluidSolver(int const& N) : n(N),
                                                         viscosity(visc),
                                                         force_mul(1.),
                                                         timer(MICROSECONDS) {
    LOG_DBG("constructing fluid solver...");
	visc = 0.001;
    alloc_buffers();
}
StamFFT_FluidSolver::~StamFFT_FluidSolver() {
    free_buffers();
}

void StamFFT_FluidSolver::use_ffts(FFT_Solver2d* fu, FFT_Solver2d* fv) {
    fftu = fu; fftv = fv;
}

float* StamFFT_FluidSolver::buff() const {return buffer;}
float* StamFFT_FluidSolver::x_buffer() const {return u;}
float* StamFFT_FluidSolver::y_buffer() const {return v;}
float* StamFFT_FluidSolver::fx_buffer() const {return u0;}
float* StamFFT_FluidSolver::fy_buffer() const {return v0;}
int const& StamFFT_FluidSolver::dim() const {return n;}

void StamFFT_FluidSolver::add_force(int x, int y, int fx, int fy) {
    if (!((x+y*n) < n*n) || (x < 0) || (y < 0) || (x >= n) || (y >= n)) {return;}//LOG_ERR("add force at %d,%d = idx %d OOB", x, y, x+y*n); return;}
	u0[x+y*n] += fx * force_mul;
	v0[x+y*n] += fy * force_mul;
}

void StamFFT_FluidSolver::set_force(int x, int y, int fx, int fy) {
    if (!((x+y*n) < n*n) || (x < 0) || (y < 0) || (x >= n) || (y >= n)) {return;}//LOG_ERR("add force at %d,%d = idx %d OOB", x, y, x+y*n); return;}
	u0[x+y*n] = fx * force_mul;
	v0[x+y*n] = fy * force_mul;
}

void StamFFT_FluidSolver::get_force(int x, int y, int& fx, int& fy) const {
    if (!((x+y*n) < n*n) || (x < 0) || (y < 0) || (x >= n) || (y >= n)) {return;}//LOG_ERR("add force at %d,%d = idx %d OOB", x, y, x+y*n); return;}
	fx = u0[x+y*n] / force_mul;
	fy = v0[x+y*n] / force_mul;
}

void StamFFT_FluidSolver::set_force_multiplier(float mul) {
    force_mul = mul;
}

float StamFFT_FluidSolver::force_multiplier() {
    return force_mul;
}

void StamFFT_FluidSolver::alloc_buffers() {
    int k;
    u = new float[n*n];
    v = new float[n*n];
    for (k = 0; k < n*n; k++) {u[k]=0.;v[k]=0.;}
    u0 = new float[n*n*2];
    v0 = new float[n*n*2];
    for (k = 0; k < n*n*2; k++) {u0[k]=0.;v0[k]=0.;}
    buffer = new float[n*n*2];
    for (k = 0; k < n*n*2; k++) {buffer[k]=0.;}
}

void StamFFT_FluidSolver::free_buffers() {
    delete [] u;
    delete [] v;
    delete [] u0;
    delete [] v0;
    delete [] buffer;
}

void StamFFT_FluidSolver::random_fill(float mag) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mag, mag);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u0[2*(i+j*n)] = dis(gen); v0[2*(i+j*n)] = dis(gen);
        }
    }
}

void StamFFT_FluidSolver::zero_field() {
    memset(u,0,sizeof(float)*n*n);
    memset(v,0,sizeof(float)*n*n);
}

void StamFFT_FluidSolver::slow_fill_pixbuff() {
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

void StamFFT_FluidSolver::step(float const& dt) {
	this->stam_stable_solve(n, u, v, u0, v0, visc, dt);
}

float StamFFT_FluidSolver::get_prev_solver_t() {
    return t_us_solver;
}

float StamFFT_FluidSolver::get_prev_fft_t() {
    return t_us_ffts;
}

#define TIMER_ACCUMULATION_START(grbg)  st = timer.read();
#define TIMER_ACCUMULATION_END(grbg)    en = timer.read(); t_us_ffts += en-st;

#define __floor(x) ((x)>=0.0?((int)(x)):(-((int)(1-(x)))))

/*
	This is Stam's solver itself
	It is almost the same as in his paper
	but adapted to use fftw3 and commented

	Inputs:
		n: 			size of n by n field
		(u, v): 	velocity of the prev time step
		(u0, v0): 	a force field defined on a grid
		visc: 		viscosity of fluid
		dt: 		time step
	Outputs:
		(u, v): 	next step of fluid solver
		(u0,v0):	cleared for new forces
*/
#define BUFF_R(buf, i, j) buf[2*(i  +(j*N))]
#define BUFF_I(buf, i, j) buf[2*(i+1+(j*N))]

void StamFFT_FluidSolver::stam_stable_solve(int const& N,
                       float* const u, float* const v,
                       float* const u0, float* const v0,
                       float const& visc, float const& dt)
{
    float x, y, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;
	float *dev_u, *dev_v, *dev_u0, *dev_v0;

	// grid size
	dim3 blockSize(256);
	dim3 gridSize((N * N + blockSize.x - 1) / blockSize.x);

    // debug code
    timer.reset_start(); float st; float en; st=en=t_us_ffts=t_us_solver=0.;


	// Assuming u, v, u0, and v0 are already allocated and initialized on the host
	cudaMalloc((void**)&dev_u, N * N * sizeof(float));
	cudaMalloc((void**)&dev_v, N * N * sizeof(float));
	cudaMalloc((void**)&dev_u0, N * N * 2 * sizeof(float)); // Times 2 for complex numbers (real and imaginary parts)
	cudaMalloc((void**)&dev_v0, N * N * 2 * sizeof(float)); // Same as above

	// Copy data from host to device
	cudaMemcpy(dev_u, u, N * N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v, v, N * N * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_u0, u0, N * N * 2 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v0, v0, N * N * 2 * sizeof(float), cudaMemcpyHostToDevice);


    // apply forces
	apply_forces_kernel<<<gridSize, blockSize>>>(dev_u, dev_v, dev_u0, dev_v0, dt, N);
	cudaDeviceSynchronize()

    // self-advection: semi-Lagrangian scheme
    // here, (u0,v0) are used to interpolate
    // and (u,v) stores the interpolation & result
	self_advection_kernel<<<gridSize, blockSize>>>(dev_u, dev_v, dev_u0, dev_v0, dt, N);
	cudaDeviceSynchronize();

    // copy velos into real-part of fourier buffers (i think)
	copy_to_fourier_buffers_kernel<<<gridSize, blockSize>>>(dev_u, dev_v, dev_u0, dev_v0, N);
	cudaDeviceSynchronize();

    // transform to fourier domain
    TIMER_ACCUMULATION_START(); // debug code
    fftu->forward(); fftv->forward();           // SPATIAL -> FOURIER
    TIMER_ACCUMULATION_END();   // debug code

    // apply low pass filter to simulate viscosity
    // and force field to be mass converving
    // by projecting vectors onto line perpendicular to wave #
    // which is line tan to circles centered at origin
	applyLowPassFilter<<<gridSize, blockSize>>>(dev_u0, dev_v0, dt, visc, N);
	cudaDeviceSynchronize();

    // inverse ffts back to spatial domain
    TIMER_ACCUMULATION_START(); // debug code
    fftu->inverse(); fftv->inverse();           // FOURIER -> SPATIAL
    TIMER_ACCUMULATION_END();   // debug code

    // normalize (r2c then c2r tform multiplies all by n*n)
	float f = 1.0/(N*N);
	normalizeAndClear<<<gridSize, blockSize>>>(dev_u, dev_v, dev_u0, dev_v0, 1.0f / (N * N), N);
	cudaDeviceSynchronize();

    // clear force field
    memset(u0,0,sizeof(float)*N*N*2);
    memset(v0,0,sizeof(float)*N*N*2);

	cudaMemcpy(u, dev_u, N * N * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(v, dev_v, N * N * sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(dev_u);
	cudaFree(dev_v);
	cudaFree(dev_u0);
	cudaFree(dev_v0);

    t_us_solver = timer.stop();
}

