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
#include "../include/bridge_delphi.h"
#include <immintrin.h>
#include <assert.h>

//----------------------
// Public API Functions
//----------------------

/**
 * Sequential matrix multiplication using the Gustavson algorithm.
 *
 * Multiplies two matrices a and b in parallel using OpenMP, storing the result in matrix c.
 *
 * @param a Pointer to matrix A (m x k).
 * @param b Pointer to matrix B (k x n).
 * @param c Pointer to matrix C (m x n).
 * @param m Number of rows in matrix A and C.
 * @param k Number of columns in matrix A and rows in matrix B.
 * @param n Number of columns in matrix B and C.
 */
BRIDGE_DELPHI_API void mm_ijk(double *a, double *b, double *c, int m, int k, int n)
{
  int i, j, p;
  double aip;

  for (i = 0; i < m; i++)
  {
    for (p = 0; p < k; p++)
    {
      aip = a[i * k + p];
      for (j = 0; j < n; j++)
        c[i * n + j] += aip * b[p * n + j];
    }
  }
}

/**
 * Parallel matrix multiplication using the Gustavson algorithm with OpenMP.
 *
 * Multiplies two matrices a and b in parallel using OpenMP, storing the result in matrix c.
 *
 * @param a Pointer to matrix A (m x k).
 * @param b Pointer to matrix B (k x n).
 * @param c Pointer to matrix C (m x n).
 * @param m Number of rows in matrix A and C.
 * @param k Number of columns in matrix A and rows in matrix B.
 * @param n Number of columns in matrix B and C.
 * @param t Number of threads to use.
 */
BRIDGE_DELPHI_API void mm_par_gustavson(double *a, double *b, double *c, int m, int k, int n, int t)
{
  int i, j, p;
  double aip;

  omp_set_num_threads(t);
#pragma omp parallel for private(i, j, p, aip) schedule(static)
  for (i = 0; i < m; i++)
  {
    for (p = 0; p < k; p++)
    {
      aip = a[i * k + p];
      for (j = 0; j < n; j++)
        c[i * n + j] += aip * b[p * n + j];
    }
  }
}

double _mm256_reduce_add_pd(__m256d vecA)
{
  __m256d t1 = _mm256_hadd_pd(vecA, vecA);

  __m128d hi = _mm256_extractf128_pd(t1, 1);
  __m128d lo = _mm256_castpd256_pd128(t1);

  __m128d sum128 = _mm_add_pd(lo, hi);

  return _mm_cvtsd_f64(sum128);
}

void transpose_matrix(double *b, double *b_t, int k, int n)
{
  for (int i = 0; i < k; i++)
    for (int j = 0; j < n; j++)
      b_t[j * k + i] = b[i * n + j];
}

BRIDGE_DELPHI_API void mm_vec(double *a, double *b, double *c, int m, int k, int n)
{
  double *b_t = (double *)malloc(k * n * sizeof(double));
  transpose_matrix(b, b_t, k, n);

  // Opcional: se puede incluir un assert para verificar la alineacion si se desea
  // assert(((uintptr_t)a % 32 == 0) && ((uintptr_t)b_t % 32 == 0) && ((uintptr_t)c % 32 == 0));

  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      // Initialize a vector register to accumulate partial sums.
      __m256d vecSum = _mm256_setzero_pd();

      int p = 0;

      // Inner loop over K, processing 4 elements at a time using AVX instructions.
      for (; p <= k - 4; p += 4)
      {
        __m256d vecA = _mm256_loadu_pd(&a[i * k + p]);   // Load 8 elements from row i of A (A[i, k:k+7])
        __m256d vecB = _mm256_loadu_pd(&b_t[j * k + p]); // Load 8 elements from column j of B (B[j, k:k+7])

        vecSum = _mm256_add_pd(vecSum, _mm256_mul_pd(vecA, vecB)); // Element-wise multiply vecA and vecB, accumulate into vecSum
      }

      // Reduce the accumulated vector sum into a single scalar value.
      double sum = _mm256_reduce_add_pd(vecSum);

      // Handle the remaining elements that do not fit into the AVX vectorized loop.
      for (; p < k; p++)
        sum += a[i * k + p] * b_t[j * k + p];

      // Store the final result in matrix C.
      c[i * n + j] = sum;
    }
  }
}

/*
BRIDGE_DELPHI_API void mm_vec(double *a, double *b, double *c, int m, int k, int n)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      // Initialize accumulator register to zero
      __m256d vecSum = _mm256_setzero_pd();

      int p = 0;
      // Process 4 elements at a time with AVX2 gather
      for (; p <= k - 4; p += 4)
      {
        // Load 4 contiguous elements from A[i, p..p+3]
        __m256d vecA = _mm256_loadu_pd(&a[i * k + p]);

        // Build index vector for B at positions:
        // B[p, j], B[p+1, j], B[p+2, j], B[p+3, j]
        // in row-major layout: b[row * n + col]
        __m256i idx = _mm256_set_epi64x(
            (long long)((p + 3) * n + j),
            (long long)((p + 2) * n + j),
            (long long)((p + 1) * n + j),
            (long long)((p + 0) * n + j));

        // Gather 4 elements from b at strides of n
        __m256d vecB = _mm256_i64gather_pd(b, idx, sizeof(double));

        // Fused multiply-add: vecSum += vecA * vecB
        vecSum = _mm256_fmadd_pd(vecA, vecB, vecSum);
      }

      // Horizontal add of vecSum into a scalar
      double sum = _mm256_reduce_add_pd(vecSum);

      // Handle any remaining elements
      for (; p < k; p++)
      {
        sum += a[i * k + p] * b[p * n + j];
      }

      // Store result into C[i, j]
      c[i * n + j] = sum;
    }
  }
}
*/