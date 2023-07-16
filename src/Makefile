LAMC.exe: main.o \
		  mathematics.o \
		  vectors.o \
		  molecule.o \
		  potential.o \
		  monte_carlo.o \
		  structural.o \
		  properties.o \
		  files_io.o \
		  simulation_setup.o
	gcc -I ./ -o ./LAMC.exe main.o \
					mathematics.o \
					vectors.o \
					molecule.o \
					potential.o \
					monte_carlo.o \
					structural.o \
					properties.o \
					simulation_setup.o \
					files_io.o -lm

files_io.o: files_io.c files_io.h mathematics.h potential.h monte_carlo.h simulation_setup.h
	gcc -I ./ -c files_io.c -lm

mathematics.o: mathematics.c mathematics.h
	gcc -I ./ -c mathematics.c -lm

molecule.o: molecule.c molecule.h constants.h mathematics.h vectors.h
	gcc -I ./ -c molecule.c -lm

monte_carlo.o: monte_carlo.c monte_carlo.h mathematics.h vectors.h molecule.h potential.h structural.h
	gcc -I ./ -c monte_carlo.c -lm

potential.o: potential.c constants.h molecule.h simulation_setup.h vectors.h
	gcc -I ./ -c potential.c -lm

properties.o: properties.c properties.h constants.h molecule.h vectors.h simulation_setup.h structural.h potential.h
	gcc -I ./ -c properties.c -lm

simulation_setup.o: simulation_setup.c simulation_setup.h
	gcc -I ./ -c simulation_setup.c -lm

structural.o: structural.c structural.h constants.h vectors.h molecule.h potential.h monte_carlo.h properties.h structural.h files_io.h
	gcc -I ./ -c structural.c -lm

vectors.o: vectors.c vectors.h mathematics.h
	gcc -I ./ -c vectors.c -lm

main.o: main.c \
		constants.h \
		vectors.o \
		molecule.o \
		potential.o \
		monte_carlo.o \
		properties.o \
		structural.o \
		files_io.o \
		simulation_setup.h
	gcc -I ./ -c main.c -lm

clean:
	rm *.o ./LAMC.exe