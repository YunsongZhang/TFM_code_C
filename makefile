cpu: CPUcode.c mypsov2.o
	gcc -g -O3 -fno-aggressive-loop-optimizations CPUcode.c mypsov2.o -lgsl -lgslcblas -lm -lrt -o cpu 
mypsov2.o:mypsov2.c
	gcc -c -O3 -fno-aggressive-loop-optimizations -g mypsov2.c -lgsl -lgslcblas -lm -o mypsov2.o
