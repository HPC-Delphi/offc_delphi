/**
 * strassen_utils.h
 *
 * Utility functions declarations for Strassen's matrix multiplication.
 *
 * (c) 2025 Daniel Galdo Galián
 */

#ifndef STRASSEN_UTILS_H
#define STRASSEN_UTILS_H

void mm_seq(double *A, double *B, double *C, int m, int k, int n);
void addMatrices(double *A, double *B, double *C, int r, int c);
void subtractMatrices(double *A, double *B, double *C, int r, int c);

#endif // STRASSEN_UTILS_H
