/**
 * strassen_utils.c
 *
 * Utility functions for Strassen's matrix multiplication.
 * Implements the recursive Strassen algorithm and helper functions.
 *
 * (c) 2025 Daniel Galdo Galián
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "../include/strassen_utils.h"

#define STRASSEN_THRESHOLD 64 // Threshold for using naive multiplication

/**
 * @brief Adds two matrices A and B of size n x n.
 *
 * A and B can have different strides. The result is stored in matrix C.
 *
 * @param A Pointer to matrix A.
 * @param strideA Row stride of matrix A.
 * @param B Pointer to matrix B.
 * @param strideB Row stride of matrix B.
 * @param C Pointer to the result matrix.
 * @param n Number of rows and columns.
 * @param strideC Row stride of matrix C.
 */
static void add_matrix(const double *A, int strideA,
                       const double *B, int strideB,
                       double *C, int n, int strideC)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            C[i * strideC + j] = A[i * strideA + j] + B[i * strideB + j];
}

/**
 * @brief Subtracts matrix B from matrix A (C = A - B) for matrices of size n x n.
 *
 * A and B can have different strides.
 *
 * @param A Pointer to matrix A.
 * @param strideA Row stride of matrix A.
 * @param B Pointer to matrix B.
 * @param strideB Row stride of matrix B.
 * @param C Pointer to the result matrix.
 * @param n Number of rows and columns.
 * @param strideC Row stride of matrix C.
 */
static void sub_matrix(const double *A, int strideA,
                       const double *B, int strideB,
                       double *C, int n, int strideC)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            C[i * strideC + j] = A[i * strideA + j] - B[i * strideB + j];
}

/**
 * @brief Performs naive matrix multiplication.
 *
 * Multiplies two matrices A and B (size n x n) and stores the result in C.
 * Matrices can have different strides.
 *
 * @param A Pointer to matrix A.
 * @param B Pointer to matrix B.
 * @param C Pointer to result matrix C.
 * @param n Number of rows and columns.
 * @param strideA Row stride of matrix A.
 * @param strideB Row stride of matrix B.
 * @param strideC Row stride of matrix C.
 */
static void naive_multiply(const double *A, const double *B, double *C,
                           int n, int strideA, int strideB, int strideC)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            double sum = 0.0;
            for (k = 0; k < n; k++)
            {
                sum += A[i * strideA + k] * B[k * strideB + j];
            }
            C[i * strideC + j] = sum;
        }
    }
}

/**
 * @brief Recursive implementation of Strassen's algorithm.
 *
 * Multiplies two submatrices A and B (of size n x n) recursively,
 * storing the result in matrix C. The parameter 'stride' indicates the
 * row width of the original matrix.
 *
 * @param A Pointer to the top-left element of submatrix A.
 * @param B Pointer to the top-left element of submatrix B.
 * @param C Pointer to the top-left element of submatrix C.
 * @param n The size of the submatrix (n x n).
 * @param stride The row stride of the original matrix.
 */
static void strassen_recursive(const double *A, const double *B, double *C,
                               int n, int stride)
{
    if (n <= STRASSEN_THRESHOLD)
    {
        // Base case: use naive multiplication.
        naive_multiply(A, B, C, n, stride, stride, stride);
        return;
    }

    int half = n / 2;
    int size = half * half;

    // Allocate memory for the seven intermediate products.
    double *M1 = (double *)malloc(size * sizeof(double));
    double *M2 = (double *)malloc(size * sizeof(double));
    double *M3 = (double *)malloc(size * sizeof(double));
    double *M4 = (double *)malloc(size * sizeof(double));
    double *M5 = (double *)malloc(size * sizeof(double));
    double *M6 = (double *)malloc(size * sizeof(double));
    double *M7 = (double *)malloc(size * sizeof(double));

    // Define submatrices (quadrants) for A.
    const double *A11 = A;
    const double *A12 = A + half;
    const double *A21 = A + half * stride;
    const double *A22 = A + half * stride + half;

    // Define submatrices for B.
    const double *B11 = B;
    const double *B12 = B + half;
    const double *B21 = B + half * stride;
    const double *B22 = B + half * stride + half;

    // Define submatrices for C.
    double *C11 = C;
    double *C12 = C + half;
    double *C21 = C + half * stride;
    double *C22 = C + half * stride + half;

    // M1 = (A11 + A22) * (B11 + B22)
    {
        double *S1 = (double *)malloc(size * sizeof(double));
        double *S2 = (double *)malloc(size * sizeof(double));
        add_matrix(A11, stride, A22, stride, S1, half, half);
        add_matrix(B11, stride, B22, stride, S2, half, half);
#pragma omp task shared(M1)
        {
            strassen_recursive(S1, S2, M1, half, half);
            free(S1);
            free(S2);
        }
    }

    // M2 = (A21 + A22) * B11
    {
        double *S = (double *)malloc(size * sizeof(double));
        add_matrix(A21, stride, A22, stride, S, half, half);
#pragma omp task shared(M2)
        {
            strassen_recursive(S, B11, M2, half, half);
            free(S);
        }
    }

    // M3 = A11 * (B12 - B22)
    {
        double *S = (double *)malloc(size * sizeof(double));
        sub_matrix(B12, stride, B22, stride, S, half, half);
#pragma omp task shared(M3)
        {
            strassen_recursive(A11, S, M3, half, half);
            free(S);
        }
    }

    // M4 = A22 * (B21 - B11)
    {
        double *S = (double *)malloc(size * sizeof(double));
        sub_matrix(B21, stride, B11, stride, S, half, half);
#pragma omp task shared(M4)
        {
            strassen_recursive(A22, S, M4, half, half);
            free(S);
        }
    }

    // M5 = (A11 + A12) * B22
    {
        double *S = (double *)malloc(size * sizeof(double));
        add_matrix(A11, stride, A12, stride, S, half, half);
#pragma omp task shared(M5)
        {
            strassen_recursive(S, B22, M5, half, half);
            free(S);
        }
    }

    // M6 = (A21 - A11) * (B11 + B12)
    {
        double *S1 = (double *)malloc(size * sizeof(double));
        double *S2 = (double *)malloc(size * sizeof(double));
        sub_matrix(A21, stride, A11, stride, S1, half, half);
        add_matrix(B11, stride, B12, stride, S2, half, half);
#pragma omp task shared(M6)
        {
            strassen_recursive(S1, S2, M6, half, half);
            free(S1);
            free(S2);
        }
    }

    // M7 = (A12 - A22) * (B21 + B22)
    {
        double *S1 = (double *)malloc(size * sizeof(double));
        double *S2 = (double *)malloc(size * sizeof(double));
        sub_matrix(A12, stride, A22, stride, S1, half, half);
        add_matrix(B21, stride, B22, stride, S2, half, half);
#pragma omp task shared(M7)
        {
            strassen_recursive(S1, S2, M7, half, half);
            free(S1);
            free(S2);
        }
    }

#pragma omp taskwait

    // Combine intermediate matrices into the result quadrants:
    // C11 = M1 + M4 - M5 + M7
    {
        double *T1 = (double *)malloc(size * sizeof(double));
        double *T2 = (double *)malloc(size * sizeof(double));
        add_matrix(M1, half, M4, half, T1, half, half);  // T1 = M1 + M4
        sub_matrix(T1, half, M5, half, T2, half, half);  // T2 = T1 - M5
        add_matrix(T2, half, M7, half, C11, half, half); // C11 = T2 + M7
        free(T1);
        free(T2);
    }

    // C12 = M3 + M5
    add_matrix(M3, half, M5, half, C12, half, half);

    // C21 = M2 + M4
    add_matrix(M2, half, M4, half, C21, half, half);

    // C22 = M1 - M2 + M3 + M6
    {
        double *T1 = (double *)malloc(size * sizeof(double));
        double *T2 = (double *)malloc(size * sizeof(double));
        sub_matrix(M1, half, M2, half, T1, half, half);  // T1 = M1 - M2
        add_matrix(T1, half, M3, half, T2, half, half);  // T2 = T1 + M3
        add_matrix(T2, half, M6, half, C22, half, half); // C22 = T2 + M6
        free(T1);
        free(T2);
    }

    // Free memory allocated for intermediate products.
    free(M1);
    free(M2);
    free(M3);
    free(M4);
    free(M5);
    free(M6);
    free(M7);
}

/**
 * @brief Public function to perform matrix multiplication using Strassen's algorithm.
 *
 * Multiplies two matrices A and B using Strassen's algorithm and stores the result in C.
 * It is assumed that n is a power of 2.
 *
 * @param A Pointer to matrix A.
 * @param B Pointer to matrix B.
 * @param C Pointer to the result matrix.
 * @param n Size of the matrices (n x n).
 * @param t Number of threads to use.
 */
void strassen_multiply(const double *A, const double *B, double *C, int n, int t)
{
    omp_set_num_threads(t);

#pragma omp parallel
    {
#pragma omp single nowait
        {
            strassen_recursive(A, B, C, n, n);
        }
    }
}
