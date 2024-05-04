Lez06.x: Lez06.c
	gcc -o Lez06.x Lez06.c

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

plot: scatter.png histogram.png Lez06.png

all: build run plot