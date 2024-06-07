CC=clang++
CFLAGS=-I/opt/homebrew/Cellar/boost/1.83.0/include -std=c++20 -O3  -L/opt/homebrew/Cellar/boost/1.83.0/lib -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp

final: test2.cpp 
	$(CC) $(CFLAGS) test2.cpp -o final2

test: test3.cpp 
	$(CC) $(CFLAGS) test3.cpp -o test

cluster: cluster.cpp
	clang++ -I/opt/homebrew/Cellar/boost/1.83.0/include -std=c++20 -O3  -L/opt/homebrew/Cellar/boost/1.83.0/lib --target=x86_64-linux-gnu cluster.cpp -o cluster.out

data: dataImport.cpp
	clang++ -O3 -I/usr/local/mysql-connector-c++-8.3.0/include/jdbc -L/usr/local/mysql-connector-c++-8.3.0/lib64 -Wl,-rpath,/usr/local/mysql-connector-c++-8.3.0/lib64 -lssl.3 -lmysqlcppconn -I/opt/homebrew/Cellar/boost/1.83.0/include -std=c++20 -O3 -L/opt/homebrew/Cellar/boost/1.83.0/lib -L/usr/local/mysql-connector-c++-8.3.0/lib64 dataImport.cpp -o test

clean: final2
	rm final2