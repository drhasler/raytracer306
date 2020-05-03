all: main

main: main.cpp
	g++ -std=c++17 -fopenmp -O3 -o o main.cpp
