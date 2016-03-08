cpu: CPUcode.c mypsov1.o
	gcc -g -O3 -fno-aggressive-loop-optimizations CPUcode.c mypsov1.o -lgsl -lgslcblas -lm -lrt -o cpu 
mypsov1.o:mypsov1.c
	gcc -c -O3 -fno-aggressive-loop-optimizations -g mypsov1.c -lgsl -lgslcblas -lm -o mypsov1.o
