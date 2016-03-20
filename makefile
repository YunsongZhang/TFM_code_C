cpu: CPUcode.c mypsov3.o
	gcc -g -O3 -fno-aggressive-loop-optimizations CPUcode.c mypsov3.o -lgsl -lgslcblas -lm -lrt -o cpu 
mypsov3.o:mypsov3.c
	gcc -c -O3 -fno-aggressive-loop-optimizations -g mypsov3.c -lgsl -lgslcblas -lm -o mypsov3.o
