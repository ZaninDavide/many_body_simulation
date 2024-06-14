Lez06.x: Lez06.c
	gcc -o Lez06.x Lez06.c -lm

build: Lez06.x

run: 
	./Lez06.x

scatter.dat:

histogram.dat:

Lez06.dat:

scatter.png: scatter.dat scatter.gnu
	gnuplot scatter.gnu

histogram.png: histogram.dat histogram.gnu
	gnuplot histogram.gnu

Lez06.png: Lez06.dat Lez06.gnu
	gnuplot Lez06.gnu

varAvgB.png: varAvgB.dat varAvgB.gnu
	gnuplot varAvgB.gnu

plot: scatter.png histogram.png Lez06.png varAvgB.png

all: build run plot