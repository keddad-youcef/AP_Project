//AVX VECTORIZATION
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particles_s {

  f32 *x, *y, *z;
  f32 *vx, *vy, *vz;
  
} particles_t;

//
static inline f32 _hadd_ps_mm256(const __m256 p) {
  __m256 t1 = _mm256_hadd_ps(p, p);
  __m256 t2 = _mm256_hadd_ps(t1, t1);
  __m128 t3 = _mm256_extractf128_ps(t2, 1);
  __m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2), t3);
  return _mm_cvtss_f32(t4);
}

//
void init(particles_t *p, u64 n)
{

  const u64 boundary = 64;  
  //
   p->x = aligned_alloc(boundary, sizeof(f32) * n);
   p->y = aligned_alloc(boundary, sizeof(f32) * n);
   p->z = aligned_alloc(boundary, sizeof(f32) * n);

   p->vx = aligned_alloc(boundary, sizeof(f32) * n);
   p->vy = aligned_alloc(boundary, sizeof(f32) * n);
   p->vz = aligned_alloc(boundary, sizeof(f32) * n);

  

  for (u64 i = 0; i < n; i++)
    {
      //
      u64 r1 = (u64)rand();
      u64 r2 = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      
      //
      p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->y[i] = (f32)rand() / (f32)RAND_MAX;
      p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;

      //
      p->vx[i] = (f32)rand() / (f32)RAND_MAX;
      p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particles_t *p, const f32 dt, u64 n)
{
  //
  __m256 softening_asm = _mm256_set1_ps(1e-20);
  __m256 dt_asm = _mm256_set1_ps(dt);
  
  //
  for (u64 i = 0; i < n; i++)
    {
      //
       __m256 fx = _mm256_setzero_ps();
       __m256 fy = _mm256_setzero_ps();
       __m256 fz = _mm256_setzero_ps();

       __m256 dxi_asm = _mm256_set1_ps(p->x[i]);
       __m256 dyi_asm = _mm256_set1_ps(p->y[i]);
       __m256 dzi_asm = _mm256_set1_ps(p->z[i]);

        //20 floating-point operations
        for (u64 j = 0; j < n; j = j+8) {

            //Newton's law
            __m256 dx = _mm256_sub_ps(_mm256_loadu_ps(p->x+j),dxi_asm);
            __m256 dy = _mm256_sub_ps(_mm256_loadu_ps(p->y+j),dyi_asm);
            __m256 dz = _mm256_sub_ps(_mm256_loadu_ps(p->z+j),dzi_asm);

            __m256 sum = _mm256_fmadd_ps(dz, dz,_mm256_fmadd_ps(dy, dy,_mm256_fmadd_ps(dx, dx, softening_asm)));

            __m256 d_2 = _mm256_rsqrt_ps (sum);
            __m256 d_3_over_2=  _mm256_mul_ps(_mm256_mul_ps(d_2,d_2),d_2);


            // Net force
            fx = _mm256_fmadd_ps(dx, d_3_over_2, fx);
            fy = _mm256_fmadd_ps(dy, d_3_over_2, fy);
            fz = _mm256_fmadd_ps(dz, d_3_over_2, fz);
	    }

        //
        p->vx[i] += dt * _hadd_ps_mm256(fx);
        p->vy[i] += dt * _hadd_ps_mm256(fy);
        p->vz[i] += dt * _hadd_ps_mm256(fz); 

    }

  //6 floating-point operations
  for (u64 i = 0; i < n; i = i+8)
    {
      __m256 pxi = _mm256_fmadd_ps(dt_asm,_mm256_loadu_ps(p->vx+i),_mm256_loadu_ps(p->x+i));
      __m256 pyi = _mm256_fmadd_ps(dt_asm,_mm256_loadu_ps(p->vy+i),_mm256_loadu_ps(p->y+i));
      __m256 pzi = _mm256_fmadd_ps(dt_asm,_mm256_loadu_ps(p->vz+i),_mm256_loadu_ps(p->z+i));

      _mm256_storeu_ps(p->x + i, pxi);
      _mm256_storeu_ps(p->y + i, pyi);
      _mm256_storeu_ps(p->z + i, pzi);
    }
}

//
int main(int argc, char **argv)
{
  //
  const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
  const u64 steps= 10;
  const f32 dt = 0.01;

  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particles_t p;



  //
  init(&p, n);

  const u64 s = sizeof(f32) * n * 6;
  
  printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  
  //
  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
  
  //
  for (u64 i = 0; i < steps; i++)
    {
      //Measure
      const f64 start = omp_get_wtime();

      move_particles(&p, dt, n);

      const f64 end = omp_get_wtime();

      //Number of interactions/iterations
      const f32 h1 = (f32)(n) * (f32)(n - 1);

      //GFLOPS
      const f32 h2 = (20.0 * h1 + 6.0 * (f32)n) * 1e-9;
      
      if (i >= warmup)
	{
	  rate += h2 / (end - start);
	  drate += (h2 * h2) / ((end - start) * (end - start));
	}
      
      //
      printf("%5llu %10.3e %10.3e %8.1f %s\n",
	     i,
	     (end - start),
	     h1 / (end - start),
	     h2 / (end - start),
	     (i < warmup) ? "*" : "");
      
      fflush(stdout);
    }
  
  //
  rate /= (f64)(steps - warmup);
  drate = sqrtf(drate / (f64)(steps - warmup) - (rate * rate));
  
  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, drate);
  printf("-----------------------------------------------------\n");
  
  //
  free(p.x);
  free(p.y);
  free(p.z);
  free(p.vx);
  free(p.vy);
  free(p.vz);

  //
  return 0;
}
