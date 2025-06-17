/*
 * bridge_delphi.h
 *
 * Header file for the OpenMP Matrix Delphi library.
 * Provides function declarations for matrix operations.
 *
 * (c) 2025 Daniel Galdo Galián
 */

#pragma once

#ifndef BRIDGE_DELPHI_EXPORTS
#define BRIDGE_DELPHI_API __declspec(dllexport)
#else
#define BRIDGE_DELPHI_API __declspec(dllimport)
#endif

/**
 * Sequential matrix multiplication using the Gustavson algorithm.
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
// BRIDGE_DELPHI_API void mm_par_gustavson(double *a, double *b, double *c, int m, int n, int k, int t);

// BRIDGE_DELPHI_API void mm_vec(double *a, double *b, double *c, int m, int n, int k);
