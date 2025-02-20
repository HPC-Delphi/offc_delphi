unit OpenMPMatrixLib;

interface

uses
  System.SysUtils;

/// <summary>
/// Performs sequential matrix multiplication.
/// </summary>
/// <remarks>
/// This function multiplies two square matrices of size N x N sequentially,
/// storing the result in the matrix C. No parallelism is used.
/// </remarks>
/// <param name="A">
/// Pointer to the first matrix (input matrix).
/// </param>
/// <param name="B">
/// Pointer to the second matrix (input matrix).
/// </param>
/// <param name="C">
/// Pointer to the result matrix (output matrix).
/// </param>
/// <param name="N">
/// The dimension of the matrices (N x N).
/// </param>
/// <example>
/// <code>
/// var
///   A, B, C: array of Double;
/// begin
///   mm_seq(@A[0], @B[0], @C[0], N);
/// end;
/// </code>
/// </example>
procedure mm_seq(A, B, C: PDouble; N: Integer); cdecl; external 'openmp_matrix_lib.dll';

/// <summary>
/// Performs parallel matrix multiplication using OpenMP.
/// </summary>
/// <remarks>
/// This function multiplies two square matrices of size N x N in parallel using OpenMP,
/// storing the result in the matrix C. The number of threads used for the computation
/// is determined by the parameter T.
/// </remarks>
/// <param name="A">
/// Pointer to the first matrix (input matrix).
/// </param>
/// <param name="B">
/// Pointer to the second matrix (input matrix).
/// </param>
/// <param name="C">
/// Pointer to the result matrix (output matrix).
/// </param>
/// <param name="N">
/// The dimension of the matrices (N x N).
/// </param>
/// <param name="T">
/// Number of threads to use in the parallel computation.
/// </param>
/// <example>
/// <code>
/// var
///   A, B, C: array of Double;
///   T: Integer;
/// begin
///   T := 4; // Use 4 threads
///   mm_par(@A[0], @B[0], @C[0], N, T);
/// end;
/// </code>
/// </example>
procedure mm_par(A, B, C: PDouble; N, T: Integer); cdecl; external 'openmp_matrix_lib.dll';

/// <summary>
/// Performs matrix multiplication using the Strassen algorithm with OpenMP.
/// </summary>
/// <remarks>
/// This function multiplies two square matrices of size N x N using the Strassen algorithm,
/// with parallelization provided by OpenMP. It is assumed that N is a power of 2.
/// </remarks>
/// <param name="A">
/// Pointer to the first matrix (input matrix).
/// </param>
/// <param name="B">
/// Pointer to the second matrix (input matrix).
/// </param>
/// <param name="C">
/// Pointer to the result matrix (output matrix).
/// </param>
/// <param name="N">
/// The dimension of the matrices (N x N). Must be a power of 2.
/// </param>
/// <param name="T">
/// Number of threads to use in the parallel computation.
/// </param>
/// <example>
/// <code>
/// var
///   A, B, C: array of Double;
///   T: Integer;
/// begin
///   T := 8; // Use 8 threads
///   mm_strassen(@A[0], @B[0], @C[0], N, T);
/// end;
/// </code>
/// </example>
procedure mm_strassen(A, B, C: PDouble; N, T: Integer); cdecl; external 'openmp_matrix_lib.dll';

implementation

end.
