#!/bin/bash

prog:
	gcc -c PSF.c -fPIC -ggdb
	gcc -shared -lgsl -lgslcblas -lm -o confocalPSF.so *.o
	
clean:
	rm confocalPSF.so *.o
