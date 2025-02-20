# openmp_matrix_lib - C Library for Matrix Operations Using OpenMP

`openmp_matrix_lib` is a dynamic library written in C that provides parallelized matrix operations, leveraging OpenMP to perform calculations efficiently on multi-core systems. The library is designed primarily for use in Delphi projects, enabling the integration of matrix multiplication and other mathematical functions with high performance.

## Features

- **Parallelized Matrix Multiplication**: The library includes optimized functions for matrix multiplication using OpenMP, allowing you to take advantage of modern multi-core processing capabilities.
- **Easy Integration in Delphi**: A **wrapper** is provided for seamless integration into Delphi applications, using the functions exposed by the DLL.
- **Optimized Performance**: Thanks to OpenMP, matrix operations are accelerated, enabling efficient handling of large datasets.

## Project Structure

The project is organized as follows:

```
openmp_matrix_lib/
│
├── src/                     # C source code for the library
│   ├── openmp_matrix_lib.c  # Implementation of the library functions
│   └── utils.c              # Auxiliary functions, if needed
│
├── include/                 # Header files for the C interface
│   └── openmp_matrix_lib.h  # Declares the exported functions
│
├── wrappers/                # Wrappers for different languages
│   └── delphi/              # Delphi wrapper
│
├── LICENSE                  # The MIT License for this project
├── Makefile                 # Build script for the DLL
└── README.md                # Project documentation
```

### Folder Description

- **`src/`**: Contains the C source code of the library. The file `openmp_matrix_lib.c` implements the mathematical functions and parallel optimizations using OpenMP. If needed, you can add auxiliary functions in `utils.c`.
  
- **`include/`**: Contains the header file `openmp_matrix_lib.h`, which declares the functions exported from the DLL for use in Delphi or any other language that may need to access the library.

- **`wrappers/`**: This folder contains the wrappers for other programming languages. Currently, it includes only the Delphi wrapper, which facilitates the integration of the DLL with a Delphi project using external declarations.

- **`Makefile`**: The `Makefile` is used to compile the library as a DLL using `gcc`. Make sure your environment is configured with GCC and supports OpenMP.

## Requirements

- **GCC** with **OpenMP** support.
- A development environment compatible with **Delphi** to integrate the DLL.

## Usage Instructions

### 1. Compiling the DLL

To compile the library as a DLL, run the following command in the project's root directory:

```bash
make
```

This command will generate the DLL (`openmp_matrix_lib.dll`), which you can use in your Delphi project.

### 2. Using in Delphi

To use the library in Delphi, follow these steps:

1. Place the DLL `openmp_matrix_lib.dll` in the same directory as your Delphi project or in an accessible path.
2. Declare the DLL functions in your Delphi code using the `external` keyword. For example:

```delphi
unit MatrixLib;

interface

uses
  System.SysUtils, System.Classes;

type
  TOpenMPMatrixLib = class
  public
    class procedure MatrixMultiply(A, B, C: PDouble; N: Integer; Threads: Integer); external 'openmp_matrix_lib.dll' name 'matrix_multiply';
  end;

implementation

end.
```

3. Call the `MatrixMultiply` function from anywhere in your Delphi application.

### 3. Example of Usage in Delphi

Here is an example of how you could use the `MatrixMultiply` function in Delphi:

```delphi
procedure TForm1.Button1Click(Sender: TObject);
var
  A, B, C: array of Double;
  N, Threads: Integer;
begin
  N := 3;  // Matrix size
  Threads := 4;  // Number of OpenMP threads
  
  SetLength(A, N * N);
  SetLength(B, N * N);
  SetLength(C, N * N);
  
  // Fill matrices A and B with data
  
  // Call the matrix multiplication function
  TOpenMPMatrixLib.MatrixMultiply(@A[0], @B[0], @C[0], N, Threads);
  
  // Use the resulting matrix C in your application
end;
```

This example multiplies two matrices, `A` and `B`, and stores the result in matrix `C`. The `Threads` parameter controls the number of OpenMP threads used to parallelize the operation.

## Contributing

If you wish to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your changes (`git checkout -b feature/new-functionality`).
3. Make your changes and commit them (`git commit -am 'Add new functionality'`).
4. Push to the branch (`git push origin feature/new-functionality`).
5. Open a pull request describing the changes made.

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

---

Thank you for using **openmp_matrix_lib**! If you have any questions or suggestions, please feel free to open an issue in the repository.