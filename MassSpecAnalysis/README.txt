This directory contains the Mass Spec analysis software being
developed by the SSEAPS team. At the time of this writing, it consists of:

- computePeptideComposition, a C program and Makefile that that will
  compute the potential compositions of a peptide given the measured
  mass to 8 significant digits. The header file describes it's usage.

- summarizeMassSpec, an R function that reads in a Mass Spec file (a
  csv), plots it, finds the peaks, and then invokes
  computePeptideComposition on the found peaks and saves the
  compositions as a returned structure.

Building

invoking make will build the computePeptideComposition program

