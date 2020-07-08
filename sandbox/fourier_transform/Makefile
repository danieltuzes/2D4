all:
	g++ fourier_transform.cpp -lboost_program_options -lfftw3 -static -Wall -Wno-unknown-pragmas -std=gnu++11 -pedantic -Ofast -flto -o "fourier_transform.bin" 2>&1 | tee error.txt