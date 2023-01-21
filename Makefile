CC=gcc
CL=clang

CFLAGS=-march=native -g3

CC_OFLAGS=-mtune=native -Ofast -fopenmp -fopt-info-all=nbody.gcc.optrpt
CL_OFLAGS=-mtune=native -Ofast -fopenmp 

all: nbodygcc nbodyclang
nbodygcc: nbody.gcc nbody_soa.gcc nbody_mem_align.gcc nbody_pow.gcc nbody_div.gcc nbody_unroll.gcc nbody_avx.gcc nbody_parallel.gcc nbody_paravx.gcc
nbodyclang: nbody.clang nbody_soa.clang nbody_mem_align.clang nbody_pow.clang nbody_div.clang nbody_unroll.clang nbody_avx.clang nbody_parallel.clang nbody_paravx.clang


#Compile

#gcc
nbody.gcc: nbody.c
	$(CC) $(CFLAGS) -O1 -fopenmp -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm

nbody_soa.gcc: nbody_soa.c
	$(CC) $(CFLAGS) $(CC_OFLAGS) $< -o $@ -lm

nbody_mem_align.gcc: nbody_mem_align.c
	$(CC) $(CFLAGS) -finline-functions -falign-functions $(CC_OFLAGS)  $< -o $@ -lm

nbody_pow.gcc: nbody_pow.c	
	$(CC) $(CFLAGS) $(CC_OFLAGS)  $< -o $@ -lm

nbody_div.gcc: nbody_div.c	
	$(CC) $(CFLAGS) $(CC_OFLAGS)  $< -o $@ -lm

nbody_unroll.gcc: nbody_unroll.c	
	$(CC) $(CFLAGS)  $(CC_OFLAGS) -msse2 -funroll-all-loops -floop-interchange $< -o $@ -lm

nbody_avx.gcc: nbody_avx.c
	$(CC) $(CFLAGS) -mavx2 -ftree-vectorize -ftree-loop-vectorize $(CC_OFLAGS)  $< -o $@ -lm

nbody_parallel.gcc: nbody_parallel.c
	$(CC) $(CFLAGS) $(CC_OFLAGS) -mavx2 -funroll-loops -ftree-vectorize -ftree-loop-vectorize $< -o $@ -lm	

nbody_paravx.gcc: nbody_paravx.c
	$(CC) $(CFLAGS) $(CC_OFLAGS) -mavx2 -funroll-loops -ftree-vectorize -ftree-loop-vectorize $< -o $@ -lm		




#clang
nbody.clang: nbody.c
	$(CL) $(CFLAGS) -O1 -fopenmp $< -o $@ -lm

nbody_soa.clang: nbody_soa.c
	$(CL) $(CFLAGS) $(CL_OFLAGS) $< -o $@ -lm

nbody_mem_align.clang: nbody_mem_align.c
	$(CL) $(CFLAGS) -finline-functions -falign-functions $(CL_OFLAGS) $< -o $@ -lm

nbody_pow.clang: nbody_pow.c	
	$(CL) $(CFLAGS) $(CL_OFLAGS) $< -o $@ -lm

nbody_div.clang: nbody_div.c	
	$(CL) $(CFLAGS) $(CL_OFLAGS) $< -o $@ -lm

nbody_unroll.clang: nbody_unroll.c	
	$(CL) $(CFLAGS)  $(CL_OFLAGS) -Wno-pass-failed -msse2 -funroll-loops $< -o $@ -lm

nbody_avx.clang: nbody_avx.c
	$(CL) $(CFLAGS) $(CL_OFLAGS) -mavx2 -ftree-vectorize $< -o $@ -lm

nbody_parallel.clang: nbody_parallel.c
	$(CL) $(CFLAGS) $(CL_OFLAGS) -mavx2 -funroll-loops -ftree-vectorize $< -o $@ -lm

nbody_paravx.clang: nbody_paravx.c		
	$(CL) $(CFLAGS) -mavx2 $(CL_OFLAGS) -mavx2 -funroll-loops -ftree-vectorize $< -o $@ -lm



#Execute

run_all: base soa mem_align pow div unroll avx parallel paravx

base:
	taskset -c 2 ./nbody.gcc
	taskset -c 2 ./nbody.clang

soa:
	taskset -c 2 ./nbody_soa.gcc
	taskset -c 2 ./nbody_soa.clang

mem_align:	
	taskset -c 2 ./nbody_mem_align.gcc
	taskset -c 2 ./nbody_mem_align.clang

pow: 
	taskset -c 4 ./nbody_pow.gcc
	taskset -c 4 ./nbody_pow.clang

div: 
	taskset -c 4 ./nbody_div.gcc
	taskset -c 4 ./nbody_div.clang

unroll: 
	taskset -c 4 ./nbody_unroll.gcc
	taskset -c 4 ./nbody_unroll.clang

avx:
	taskset -c 4 ./nbody_avx.gcc
	taskset -c 4 ./nbody_avx.clang

parallel: 
	./nbody_parallel.gcc
	./nbody_parallel.clang

paravx:
	./nbody_paravx.gcc
	./nbody_paravx.clang	




#Clean
clean:
	rm -Rf *.gcc *.clang *.optrpt



