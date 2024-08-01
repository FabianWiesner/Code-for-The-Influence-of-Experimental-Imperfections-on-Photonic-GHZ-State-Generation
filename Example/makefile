# This is an example make file for a macOS build system. You may change the compiler from clang++ to g++ (and delete the flag -Xclang) if you are not using macOS. 
# You may change paths to boost and openmp. 
# Further, openmp is not required - however, we recommend using it. If you prefer not to use it you should delete the flags after and including -fopenmp. 

CC=clang++
CFLAGS=-I/opt/homebrew/Cellar/boost/1.83.0/include -std=c++20 -O3  -I../ -L/opt/homebrew/Cellar/boost/1.83.0/lib -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp

# In this make file we assume you'd like to compile a file named test.cpp. 
final: test.cpp 
	$(CC) $(CFLAGS) test.cpp -o test

