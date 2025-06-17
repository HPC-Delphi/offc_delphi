#pragma once

#ifndef BRIDGE_DELPHI_EXPORTS
#define BRIDGE_DELPHI_API __declspec(dllexport)
#else
#define BRIDGE_DELPHI_API __declspec(dllimport)
#endif

/**
 * Sequential matrix multiplication using ijk matrix multiplication algorithm.
 *
 * Multiplies two matrices a and b, storing the result in matrix c.
 *
 * @param a Pointer to matrix A (m x k).
 * @param b Pointer to matrix B (k x n).
 * @param c Pointer to matrix C (m x n).
 * @param m Number of rows in matrix A and C.
 * @param k Number of columns in matrix A and rows in matrix B.
 * @param n Number of columns in matrix B and C.
 */
BRIDGE_DELPHI_API void mm_ijk(double *a, double *b, double *c, int m, int k, int n);

/**
 * Sequential matrix multiplication using the Strassen algorithm.
 *
 * Multiplies two matrices a and b, storing the result in matrix c.
 * It is assumed that m, k, and n are powers of 2.
 *
 * @param a Pointer to matrix A (m x k).
 * @param b Pointer to matrix B (k x n).
 * @param c Pointer to matrix C (m x n).
 * @param m Number of rows in matrix A and C.
 * @param k Number of columns in matrix A and rows in matrix B.
 * @param n Number of columns in matrix B and C.
 */
// BRIDGE_DELPHI_API void mm_strassen(double *a, double *b, double *c, int m, int k, int n);

// BRIDGE_DELPHI_API void mm_par_gustavson(double *a, double *b, double *c, int m, int n, int k, int t);
// BRIDGE_DELPHI_API void mm_vec(double *a, double *b, double *c, int m, int n, int k);
