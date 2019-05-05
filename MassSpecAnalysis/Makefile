# The world's simplest Makefile
OS := $(shell uname -s)

ifeq ($OS,Linux)
	CC = cc
	CFLAGS = -pthread -lm
endif
ifeq ($OS,Darwin)
	CC = gcc
	CFLAGS = -lpthread
endif

default:	computeParallelPeptideComposition

computePeptideComposition: computePeptideComposition.c
	$(CC) -o computePeptideComposition computePeptideComposition.c $(CFLAGS)

computeParalelPeptideComposition: computeParallelPeptideComposition.c
	$(CC) -o computeParallelPeptideComposition computeParallelPeptideComposition.c $(CFLAGS)

