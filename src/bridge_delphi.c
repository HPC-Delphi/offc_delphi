#include <stdio.h>
#include <stdlib.h>
#include "../include/bridge_delphi.h"
// #include <immintrin.h>
// #include <assert.h>
// #include <omp.h>

/*
#define THRESHOLD 64
// ----------------------
// Utility Functions
// ----------------------
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
*/

//----------------------
// Public API Functions
//----------------------
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

/*
BRIDGE_DELPHI_API void mm_strassen(double *a, double *b, double *c, int m, int k, int n)
{
  if (m <= THRESHOLD || k <= THRESHOLD || n <= THRESHOLD) // Base case
  {

    mm_ijk(a, b, c, m, k, n);
    return;
  }

  // Dividir las matrices A y B en submatrices
  int m2 = m / 2; // Mitad de las filas de A
  int k2 = k / 2; // Mitad de las columnas de A y filas de B
  int n2 = n / 2; // Mitad de las columnas de B

  double *A11 = (double *)malloc(m2 * k2 * sizeof(double));
  double *A12 = (double *)malloc(m2 * (k - k2) * sizeof(double));
  double *A21 = (double *)malloc((m - m2) * k2 * sizeof(double));
  double *A22 = (double *)malloc((m - m2) * (k - k2) * sizeof(double));

  double *B11 = (double *)malloc(k2 * n2 * sizeof(double));
  double *B12 = (double *)malloc(k2 * (n - n2) * sizeof(double));
  double *B21 = (double *)malloc((k - k2) * n2 * sizeof(double));
  double *B22 = (double *)malloc((k - k2) * (n - n2) * sizeof(double));

  // Submatrices of A
  for (int i = 0; i < m2; i++)
    for (int j = 0; j < k2; j++)
      A11[i * k2 + j] = a[i * k + j];

  for (int i = 0; i < m2; i++)
    for (int j = 0; j < k - k2; j++)
      A12[i * (k - k2) + j] = a[i * k + k2 + j];

  for (int i = 0; i < m - m2; i++)
    for (int j = 0; j < k2; j++)
      A21[i * k2 + j] = a[(m2 + i) * k + j];

  for (int i = 0; i < m - m2; i++)
    for (int j = 0; j < k - k2; j++)
      A22[i * (k - k2) + j] = a[(m2 + i) * k + k2 + j];

  // Submatrices of B
  for (int i = 0; i < k2; i++)
    for (int j = 0; j < n2; j++)
      B11[i * n2 + j] = b[i * n + j];

  for (int i = 0; i < k2; i++)
    for (int j = 0; j < n - n2; j++)
      B12[i * (n - n2) + j] = b[i * n + n2 + j];

  for (int i = 0; i < k - k2; i++)
    for (int j = 0; j < n2; j++)
      B21[i * n2 + j] = b[(k2 + i) * n + j];

  for (int i = 0; i < k - k2; i++)
    for (int j = 0; j < n - n2; j++)
      B22[i * (n - n2) + j] = b[(k2 + i) * n + n2 + j];

  // Calcular los productos M1 a M7 usando Strassen
  double *M1 = (double *)calloc(m2 * n2 * sizeof(double), sizeof(double));
  double *M2 = (double *)calloc(m2 * (n - n2) * sizeof(double), sizeof(double));
  double *M3 = (double *)calloc((m - m2) * n2 * sizeof(double), sizeof(double));
  double *M4 = (double *)calloc((m - m2) * (n - n2) * sizeof(double), sizeof(double));
  double *M5 = (double *)calloc(m2 * (n - n2) * sizeof(double), sizeof(double));
  double *M6 = (double *)calloc((m - m2) * n2 * sizeof(double), sizeof(double));
  double *M7 = (double *)calloc((m - m2) * (n - n2) * sizeof(double), sizeof(double));

  // M1 = (A11 + A22) * (B11 + B22)
  double *A11_plus_A22 = (double *)malloc(m2 * k2 * sizeof(double));
  double *B11_plus_B22 = (double *)malloc(k2 * n2 * sizeof(double));
  addMatrices(A11, A22, A11_plus_A22, m2, k2);
  addMatrices(B11, B22, B11_plus_B22, k2, n2);
  mm_strassen(A11_plus_A22, B11_plus_B22, M1, m2, k2, n2);

  // M2 = (A21 + A22) * B11
  double *A21_plus_A22 = (double *)malloc((m - m2) * k2 * sizeof(double));
  addMatrices(A21, A22, A21_plus_A22, m - m2, k2);
  mm_strassen(A21_plus_A22, B11, M2, m - m2, k2, n2);

  // M3 = A11 * (B12 - B22)
  double *B12_minus_B22 = (double *)malloc(k2 * (n - n2) * sizeof(double));
  subtractMatrices(B12, B22, B12_minus_B22, k2, n - n2);
  mm_strassen(A11, B12_minus_B22, M3, m2, k2, n - n2);

  // M4 = A22 * (B21 - B11)
  double *B21_minus_B11 = (double *)malloc((k - k2) * n2 * sizeof(double));
  subtractMatrices(B21, B11, B21_minus_B11, k - k2, n2);
  mm_strassen(A22, B21_minus_B11, M4, m - m2, k - k2, n2);

  // M5 = (A11 + A12) * B22
  double *A11_plus_A12 = (double *)malloc(m2 * (k - k2) * sizeof(double));
  addMatrices(A11, A12, A11_plus_A12, m2, k - k2);
  mm_strassen(A11_plus_A12, B22, M5, m2, k - k2, n - n2);

  // M6 = (A21 - A11) * (B11 + B12)
  double *A21_minus_A11 = (double *)malloc((m - m2) * k2 * sizeof(double));
  subtractMatrices(A21, A11, A21_minus_A11, m - m2, k2);
  double *B11_plus_B12 = (double *)malloc(k2 * (n - n2) * sizeof(double));
  addMatrices(B11, B12, B11_plus_B12, k2, n - n2);
  mm_strassen(A21_minus_A11, B11_plus_B12, M6, m - m2, k2, n - n2);

  // M7 = (A12 - A22) * (B21 + B22)
  double *A12_minus_A22 = (double *)malloc(m2 * (k - k2) * sizeof(double));
  subtractMatrices(A12, A22, A12_minus_A22, m2, k - k2);
  double *B21_plus_B22 = (double *)malloc((k - k2) * n2 * sizeof(double));
  addMatrices(B21, B22, B21_plus_B22, k - k2, n2);
  mm_strassen(A12_minus_A22, B21_plus_B22, M7, m2, k - k2, n2);

  // Combinar los Cados de M1 a M7 para formar la matriz final C
  // C11 = M1 + M4 - M5 + M7
  double *C11 = (double *)malloc(m2 * n2 * sizeof(double));
  addMatrices(M1, M4, C11, m2, n2);
  subtractMatrices(C11, M5, C11, m2, n2);
  addMatrices(C11, M7, C11, m2, n2);

  // C12 = M3 + M5
  double *C12 = (double *)malloc(m2 * (n - n2) * sizeof(double));
  addMatrices(M3, M5, C12, m2, n - n2);

  // C21 = M2 + M4
  double *C21 = (double *)malloc((m - m2) * n2 * sizeof(double));
  addMatrices(M2, M4, C21, m - m2, n2);

  // C22 = M1 - M2 + M3 + M6
  double *C22 = (double *)malloc((m - m2) * (n - n2) * sizeof(double));
  addMatrices(M1, M6, C22, m - m2, n - n2);
  subtractMatrices(C22, M2, C22, m - m2, n - n2);
  addMatrices(C22, M3, C22, m - m2, n - n2);

  // Combinar C11, C12, C21, C22 en la matriz C final
  for (int i = 0; i < m2; i++)
    for (int j = 0; j < n2; j++)
      c[i * n + j] = C11[i * n2 + j];

  for (int i = 0; i < m2; i++)
    for (int j = 0; j < n - n2; j++)
      c[i * n + n2 + j] = C12[i * (n - n2) + j];

  for (int i = 0; i < m - m2; i++)
    for (int j = 0; j < n2; j++)
      c[(m2 + i) * n + j] = C21[i * n2 + j];

  for (int i = 0; i < m - m2; i++)
    for (int j = 0; j < n - n2; j++)
      c[(m2 + i) * n + n2 + j] = C22[i * (n - n2) + j];

  // Free used memory
  free(A11);
  free(A12);
  free(A21);
  free(A22);
  free(B11);
  free(B12);
  free(B21);
  free(B22);
  free(M1);
  free(M2);
  free(M3);
  free(M4);
  free(M5);
  free(M6);
  free(M7);
  free(C11);
  free(C12);
  free(C21);
  free(C22);
}
*/

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
/*
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
}*/

/*
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
*/
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