#ifndef STRASSEN_UTILS_H
#define STRASSEN_UTILS_H

/**
 * @brief Performs matrix multiplication using Strassen's algorithm.
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
void strassen_multiply(const double *A, const double *B, double *C, int n, int t);

#endif // STRASSEN_UTILS_H
