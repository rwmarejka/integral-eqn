#!/bin/sh

LIBS=-lm
cc -c int.c
cc -o ex0 ex0.c int.o $LIBS
cc -o ex1 ex1.c int.o $LIBS
cc -o ex2 ex2.c int.o $LIBS
cc -o ex3 ex3.c int.o $LIBS
cc -o ex4 ex4.c int.o $LIBS
cc -o ex5 ex5.c int.o $LIBS
cc -o ex6 ex6.c int.o $LIBS
cc -o ex7 ex7.c int.o $LIBS
cc -o ex8 ex8.c int.o $LIBS
cc -o ex9 ex9.c int.o $LIBS
cc -o ex10 ex10.c int.o $LIBS
cc -o ex11 ex11.c int.o $LIBS
