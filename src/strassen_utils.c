/**
 * strassen_utils.c
 *
 * Utility functions for Strassen's matrix multiplication.
 *
 * (c) 2025 Daniel Galdo Galián
 */

#include <stdio.h>
#include <stdlib.h>
#include "../include/strassen_utils.h"

#define THRESHOLD 64

void mm_seq(double *A, double *B, double *C, int m, int k, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i * n + j] = 0;
            for (int z = 0; z < k; z++)
                C[i * n + j] += A[i * k + z] * B[z * n + j];
        }
    }
}

void addMatrices(double *A, double *B, double *C, int r, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            C[i * c + j] = A[i * c + j] + B[i * c + j];
}

void subtractMatrices(double *A, double *B, double *C, int r, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            C[i * c + j] = A[i * c + j] - B[i * c + j];
}