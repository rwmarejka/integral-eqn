#!/bin/sh

LIBS=-lm
cc -c int.c
cc -o ex00 ex00.c int.o $LIBS
cc -o ex01 ex01.c int.o $LIBS
cc -o ex02 ex02.c int.o $LIBS
cc -o ex11 ex11.c int.o $LIBS
cc -o ex12 ex12.c int.o $LIBS
cc -o ex13 ex13.c int.o $LIBS
cc -o ex21 ex21.c int.o $LIBS
cc -o ex22 ex22.c int.o $LIBS
cc -o ex23 ex23.c int.o $LIBS
cc -o ex24 ex24.c int.o $LIBS
cc -o ex25 ex25.c int.o $LIBS
cc -o ex31 ex31.c int.o $LIBS
cc -o ex41 ex41.c int.o $LIBS
cc -o ex42 ex42.c int.o $LIBS
cc -o ex43 ex43.c int.o $LIBS
cc -o ex51 ex51.c int.o $LIBS
cc -o ex52 ex52.c int.o $LIBS