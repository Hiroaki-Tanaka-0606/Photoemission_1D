all: lapack_test

lapack_test: lapack_test.cpp
	g++ lapack_test.cpp -L/source/lapack/lapack-3.12.0 -llapack -lrefblas -ldl -lm -lgfortran -lquadmath -static -mcmodel=medium -o lapack_test.o