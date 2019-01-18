#!/bin/bash

gcc -c ODEmodel.c ODEmodel2.c
g++ -c SolveMain.cpp MicroEnv.cpp BMSC.cpp MIC.cpp MM.cpp CTL.cpp Treg.cpp Location.cpp par.cpp Utilities.cpp PC.cpp TMM.cpp -std=c++0x 

gfortran -c ./libf77/opkda1.f -o ./libf77/opkda1.o
gfortran -c ./libf77/opkda2.f -o ./libf77/opkda2.o
gfortran -c ./libf77/opkdmain.f -o ./libf77/opkdmain.o
ar -r ./libf77/libodepack.a ./libf77/opkda1.o ./libf77/opkda2.o ./libf77/opkdmain.o

gfortran SolveMain.o MicroEnv.o BMSC.o MIC.o MM.o CTL.o Treg.o Location.o par.o Utilities.o PC.o TMM.o ODEmodel.o ODEmodel2.o -o run_main -L./libf77 -lodepack -lm -lstdc++

rm SolveMain.o MicroEnv.o BMSC.o MIC.o MM.o CTL.o Treg.o Location.o par.o Utilities.o PC.o TMM.o ODEmodel.o ODEmodel2.o

#./run_main


