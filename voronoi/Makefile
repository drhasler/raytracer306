all: o

OPT=-Iincl -O3 -Wall -fopenmp

o: main.cpp main.h
	g++ -std=c++17 main.cpp $(OPT) -o o

report:
	inkscape --export-pdf pics/petri.pdf pics/petri.svg
	bash -c "pandoc <(echo '# Report';\
	sed 's/\.svg/.pdf/' README.md) -o report.pdf"

animate:
	convert -delay 1x8 pics/frame*.png pics/anim.gif

