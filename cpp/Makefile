all: lapack_test photoemission_1d photoemission_well

lapack_test: lapack_test.cpp
	g++ lapack_test.cpp -L/source/lapack/lapack-3.12.0 -llapack -lrefblas -ldl -lm -lgfortran -lquadmath -static -mcmodel=medium -o lapack_test.o

photoemission_1d: photoemission_1d.cpp
	g++ photoemission_1d.cpp -L/source/lapack/lapack-3.12.0 -llapack -lrefblas -ldl -lm -lgfortran -lquadmath -static -mcmodel=medium -o photoemission_1d.o

photoemission_well: photoemission_well.cpp
	g++ photoemission_well.cpp -static -o photoemission_well.o