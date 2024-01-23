CC=g++
CFLAGS=-I/usr/include/gsl 
LDFLAGS=-lm `gsl-config --libs`

main: anim_fun.o forces.o functions.o matrix_ctors.o matrix_fun.o matrix_LT.o set_ctors.o set_fun.o set_statistic.o set_tools.o sim_manager.o
	${CC} ${CFLAGS} -c main.cpp
	${CC} ${CFLAGS} ${LDFLAGS} anim_fun.o forces.o functions.o matrix_ctors.o matrix_fun.o matrix_LT.o set_ctors.o set_fun.o set_statistic.o set_tools.o sim_manager.o main.o -o spatcorr

anim_fun.o:
	${CC} ${CFLAGS} -c anim_fun.cpp
functions.o:
	${CC} ${CFLAGS} -c functions.cpp
matrix_ctors.o:
	${CC} ${CFLAGS} -c matrix_ctors.cpp
matrix_fun.o:
	${CC} ${CFLAGS} -c matrix_fun.cpp
matrix_LT.o:
	${CC} ${CFLAGS} -c matrix_LT.cpp
set_ctors.o:
	${CC} ${CFLAGS} -c set_ctors.cpp
set_fun.o:
	${CC} ${CFLAGS} -c set_fun.cpp
set_statistic.o:
	${CC} ${CFLAGS} -c set_statistic.cpp
set_tools.o:
	${CC} ${CFLAGS} -c set_tools.cpp
forces.o:
	${CC} ${CFLAGS} -c forces.cpp
sim_manager.o:
	${CC} ${CFLAGS} -c sim_manager.cpp
clean:
	rm -f spatcorr main.o anim_fun.o forces.o functions.o matrix_ctors.o matrix_fun.o matrix_LT.o set_ctors.o set_fun.o set_statistic.o set_tools.o sim_manager.o