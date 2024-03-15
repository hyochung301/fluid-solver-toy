
main.o: src/main.cpp
	clang++ -c -o main.o src/main.cpp -Ilib/

run: main.o
	clang++ -o run main.o lib/libfftw3.a     

clean:
	rm main.o run
