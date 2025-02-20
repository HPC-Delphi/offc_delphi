/*
 * openmp_matrix_lib.c
 *
 * Implementation of the OpenMP Matrix Library functions.
 * Provides various matrix multiplication implementations using OpenMP.
 *
 * (c) 2025 Daniel Galdo Galián
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "../include/openmp_matrix_lib.h"
#include "../include/strassen_utils.h"

//----------------------
// Public API Functions
//----------------------

/**
 * Parallel matrix multiplication using OpenMP.
 *
 * Multiplies two matrices a and b in parallel using OpenMP, storing the result in matrix c.
 *
 * @param a Pointer to the first matrix.
 * @param b Pointer to the second matrix.
 * @param c Pointer to the result matrix.
 * @param n Size of the matrices (n x n).
 * @param t Number of threads to use.
 */
OPENMP_MATRIX_LIB_API void mm_naive(double *a, double *b, double *c, int n, int t)
{
  int i, j, k;
  double aik;

  omp_set_num_threads(t);
#pragma omp parallel for private(i, j, k, aik) schedule(static)
  for (i = 0; i < n; i++)
  {
    for (k = 0; k < n; k++)
    {
      aik = a[i * n + k];
      for (j = 0; j < n; j++)
      {
        c[i * n + j] += aik * b[k * n + j];
      }
    }
  }
}

/**
 * Matrix multiplication using the Strassen algorithm with OpenMP.
 *
 * Multiplies two matrices a and b using the Strassen algorithm and OpenMP for parallelization,
 * storing the result in matrix c. It is assumed that n is a power of 2.
 *
 * @param a Pointer to the first matrix.
 * @param b Pointer to the second matrix.
 * @param c Pointer to the result matrix.
 * @param n Size of the matrices (n x n).
 * @param t Number of threads to use.
 */
OPENMP_MATRIX_LIB_API void mm_strassen(double *a, double *b, double *c, int n, int t)
{
  strassen_multiply(a, b, c, n, t);
}
