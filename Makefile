all: main

main: main.cpp
	g++ -std=c++17 -fopenmp -O2 -o o main.cpp
