# openmp_matrix_delphi - High-Performance Matrix Operations Library

`openmp_matrix_delphi` is a dynamic library written in C that provides efficient matrix operations, leveraging OpenMP for parallelization. Designed for integration with Delphi applications, this library enables high-performance matrix computations, including sequential and parallel implementations of Gustavson and Strassen algorithms.

## Features

- **Parallelized Matrix Multiplication**: Optimized implementations of matrix multiplication using OpenMP, enabling efficient utilization of multi-core processors.
- **Delphi Integration**: A Delphi wrapper is provided to seamlessly integrate the library into Delphi projects.
- **Support for Advanced Algorithms**: Includes both Gustavson and Strassen algorithms for matrix multiplication, with sequential and parallel variants.
- **Cross-Language Compatibility**: The library is implemented as a DLL, making it accessible from various programming environments.

## Development and Compilation Environment

The library has been developed and tested in the following environment:

- **Operating System**: Windows 10
- **Toolchain**: UCRT64, available in the MSYS2 collection of tools and libraries
- **Compiler**: GCC (UCRT version), a modern implementation of the C standard library for Windows, offering better compatibility and performance compared to MSVCRT
- **Delphi Wrapper Development**: RAD Studio (Delphi 12.1 Community Edition), an integrated development environment (IDE) providing tools for designing, compiling, and debugging Delphi applications

## Project Structure

The project is organized as follows:

```
openmp_matrix_delphi/
│
├── src/                     # C source code for the library
│   ├── openmp_matrix_delphi.c  # Core library implementation
│   └── strassen_utils.c     # Utility functions for Strassen's algorithm
│
├── include/                 # Header files for the library
│   ├── openmp_matrix_delphi.h  # Public API declarations
│   └── strassen_utils.h     # Utility function declarations
│
├── wrappers/                # Wrappers for other languages
│   └── OpenMPMatrix.pas     # Delphi wrapper
│
├── LICENSE                  # License information (MIT)
├── Makefile                 # Build script for the library
└── README.md                # Project documentation
```

### Folder Description

- **`src/`**: Contains the C source code for the library, including the main implementation and utility functions.
- **`include/`**: Header files defining the public API and utility functions.
- **`wrappers/`**: Contains the Delphi wrapper (`OpenMPMatrix.pas`) for integrating the library into Delphi projects.
- **`Makefile`**: A build script for compiling the library as a DLL using GCC with OpenMP support.

## Requirements

- **MSYS2**: Install the UCRT64 environment from the MSYS2 toolchain to compile the library.
- **GCC with OpenMP Support**: Ensure that the GCC compiler supports OpenMP for parallelization.
- **RAD Studio (Delphi 12.1 Community Edition)**: Required for developing and testing the Delphi wrapper.

## Compilation Instructions

To compile the library as a DLL, follow these steps:

1. Open the MSYS2 terminal and switch to the UCRT64 environment.
2. Navigate to the project root directory.
3. Run the following command:

```bash
make
```

This will generate the `openmp_matrix_delphi.dll` file, along with the corresponding import library (`openmp_matrix_delphi.a`).

## Using the Library in Delphi

### 1. Setup

- Place the compiled `openmp_matrix_delphi.dll` in the same directory as your Delphi project or in a directory included in the system's PATH.

### 2. Using the Wrapper in a Delphi Project

The Delphi wrapper (`OpenMPMatrix.pas`) is provided to simplify the integration of the library into your Delphi project. Follow these steps to use it:

1. **Add the Wrapper to Your Project**  
   Copy the `OpenMPMatrix.pas` file into your Delphi project directory and include it in your project. You can do this by adding the file to your Delphi project in the IDE.

2. **Adapt the DLL Path**  
   Open the `OpenMPMatrix.pas` file and ensure the `DLLPATH` constant points to the correct location of the `openmp_matrix_delphi.dll` file. For example:
   ```delphi
   const
     DLLPATH = 'path_to_your_dll\openmp_matrix_delphi.dll';
   ```

3. **Call the Functions**  
   Use the functions declared in the wrapper to perform matrix operations. Here's an example:

   ```delphi
   procedure TForm1.Button1Click(Sender: TObject);
   var
     A, B, C: array of Double;
     M, K, N, Threads: Integer;
   begin
     M := 4; K := 4; N := 4; // Matrix dimensions
     Threads := 4; // Number of threads for parallel execution

     SetLength(A, M * K);
     SetLength(B, K * N);
     SetLength(C, M * N);

     // Initialize matrices A and B with values
     A[0] := 1; A[1] := 2; A[2] := 3; A[3] := 4;
     B[0] := 5; B[1] := 6; B[2] := 7; B[3] := 8;

     // Call the parallel Gustavson multiplication function
     MMParGustavson(@A[0], @B[0], @C[0], M, K, N, Threads);

     // Use the result matrix C
     ShowMessage(Format('C[0]: %f', [C[0]]));
   end;
   ```

This example demonstrates how to initialize matrices, call the parallel Gustavson multiplication function, and display the result. Ensure that the `DLLPATH` in the wrapper is correctly set to avoid runtime errors.

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature/new-feature`).
3. Commit your changes (`git commit -m 'Add new feature'`).
4. Push to the branch (`git push origin feature/new-feature`).
5. Open a pull request.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.