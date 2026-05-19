# offc_delphi - High Performance Matrix Operations Library

[![Platform: Windows](https://img.shields.io/badge/Platform-Windows%2011-blue)](https://www.microsoft.com/en-us/windows/windows-11)
[![GCC](https://img.shields.io/badge/GCC-15.1.0-brightgreen)](https://winlibs.com/)

`offc_delphi` is a dynamic library written in C, providing efficient matrix multiplication operations for Delphi applications. It implements various advanced algorithms including sequential, OpenMP parallelized, SIMD vectorized (AVX2/FMA), and Strassen's matrix multiplication. This module is designed for High-Performance Computing (HPC) tasks requiring intensive algebraic computations.

## Features

- **Matrix Multiplication Algorithms**: Implements standard sequential (IKJ), OpenMP parallelized, and Strassen's recursive multiplication.
- **SIMD Vectorization**: Utilizes AVX and AVX2 intrinsics along with loop unrolling to maximize CPU throughput in vector operations.
- **Delphi Wrapper**: Includes a Delphi wrapper for seamless integration.
- **High-Performance Computing**: Enables vectorized/optimized computations for improved performance.
- **Cross-Language Compatibility**: Distributed as a DLL, accessible from various programming environments.
- **Modern Toolchain Support**: Built and tested with the latest GCC toolchain for Windows 11.

## Requirements & Development Environment

The library is developed, tested, and intended to be used in the following environment:

- **Operating System**: Windows 11 (64-bit)
- **Compiler**: GCC 15.1.0
- **Toolchain**: MinGW-w64 12.0.0 UCRT (release 1) from winlibs.com
- **Dependencies**: OpenMP (included in GCC) for multi-threading.
- **Build System**: `mingw32-make`
- **Delphi Integration**: RAD Studio (Delphi 12.1 Community Edition)

## Project Structure

```text
offc_delphi/
│
├── build/                  # Compiled binaries and intermediate files
├── include/                # Public API headers (C)
│   └── offc_delphi.h
├── src/                    # C source code
│   └── offc_delphi.c
├── wrappers/               # Delphi language wrapper
│   └── OffC.pas
├── LICENSE                 # License information
├── Makefile                # Build script for DLL
└── README.md               # Project documentation
```

## Compilation Instructions

To build the library as a DLL:

1. Ensure the GCC toolchain with OpenMP support is available in your system's PATH.
2. Open a Command Prompt or PowerShell window.
3. Navigate to the project root directory (`offc_delphi`).
4. Run:

   ```sh
   mingw32-make
   ```

This will generate `offc_delphi.dll` in the `build/` directory.

## Using the Library in Delphi

### 1. Setup
- Place `offc_delphi.dll` in the directory where the executable is created upon compiling your Delphi project, or ensure it is in a directory included in your system's PATH.

### 2. Using the Wrapper in a Delphi Project
- Add the path to `OffC.pas` in the "Search Path" under Project Options in your Delphi project.
- Use the functions declared in the wrapper to perform operations.

**Example Usage:**
```pascal
uses OffC;

// Multiply matrices A and B into C using sequential IKJ algorithm
mm_ikj_seq(@MatrixA[0], @MatrixB[0], @MatrixC[0], Rows, K, Cols);
```

## Academic Citation

If you use this software in your research, please cite it using the metadata provided in the `CITATION.cff` file located in the root of this repository.

## License

This project is licensed under the MIT License - see the LICENSE file for details.