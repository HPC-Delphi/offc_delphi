/*
 * openmp_matrix_lib.h
 *
 * Header file for the OpenMP Matrix Library.
 * Provides function declarations for matrix operations.
 *
 * (c) 2025 Daniel Galdo Galián
 */

#pragma once

#ifndef OPENMP_MATRIX_LIB_EXPORTS
#define OPENMP_MATRIX_LIB_API __declspec(dllexport)
#else
#define OPENMP_MATRIX_LIB_API __declspec(dllimport)
#endif

/**
 * Matrix multiplication using the Gustavson algorithm with OpenMP.
 *
 * Multiplies two matrices a and b in parallel using OpenMP, storing the result in matrix c.
 *
 * @param a Pointer to the first matrix.
 * @param b Pointer to the second matrix.
 * @param c Pointer to the result matrix.
 * @param n Size of the matrices (n x n).
 * @param t Number of threads to use.
 */
OPENMP_MATRIX_LIB_API void mm_gustavson(double *a, double *b, double *c, int n, int t);

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
OPENMP_MATRIX_LIB_API void mm_strassen(double *a, double *b, double *c, int n, int t);
