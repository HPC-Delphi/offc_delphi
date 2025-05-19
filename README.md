# omp_delphi - High Performance Matrix Operations Library

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform: Windows](https://img.shields.io/badge/Platform-Windows%2011-blue)](https://www.microsoft.com/en-us/windows/windows-11)
[![GCC](https://img.shields.io/badge/GCC-15.1.0-brightgreen)](https://winlibs.com/)

`omp_delphi` is a high-performance dynamic library written in C, providing efficient matrix operations with parallelization via OpenMP. Designed for seamless integration with Delphi applications, it enables advanced matrix computations, including both sequential and parallel implementations of Gustavson and Strassen algorithms. The library is distributed as a DLL for easy cross-language usage.

---

## Features

- **Parallelized Matrix Multiplication**: Optimized matrix multiplication using OpenMP for multi-core performance.
- **Delphi Integration**: Includes a Delphi wrapper for straightforward use in Delphi projects.
- **Advanced Algorithms**: Implements Gustavson and Strassen algorithms (sequential and parallel variants).
- **Cross-Language Compatibility**: Distributed as a DLL, accessible from various programming environments.
- **Modern Toolchain Support**: Built and tested with the latest GCC and MinGW-w64 toolchains for Windows 11.

---

## Requirements & Development Environment

The library is developed, tested, and intended to be used in the following environment:

- **Operating System**: Windows 11 (64-bit)
- **Compiler**: GCC 15.1.0 (with POSIX threads)
- **Toolchain**: MinGW-w64 12.0.0 UCRT (release 1) from [winlibs.com](https://winlibs.com/)
- **OpenMP support** (included in the above GCC toolchain)
- **Build System**: `mingw32-make`
- **Delphi Integration**: RAD Studio (Delphi 12.1 Community Edition)

---

## Project Structure

```
omp_delphi/
│
├── build/                  # Compiled binaries and intermediate files
├── include/                # Public API and utility headers
│   ├── omp_delphi.h
│   └── strassen_utils.h
├── src/                    # C source code
│   ├── omp_delphi.c
│   └── strassen_utils.c
├── wrappers/               # Language wrappers (Delphi)
│   └── OMPMatrix.pas
├── LICENSE                 # License information (MIT)
├── Makefile                # Build script for DLL
└── README.md               # Project documentation
```

---

## Compilation Instructions

To build the library as a DLL:

1. Download and extract the [winlibs.com](https://winlibs.com/) GCC toolchain (ensure it is in your PATH).
2. Open a Command Prompt or PowerShell window.
3. Navigate to the project root directory (`omp_delphi`).
4. Run:

   ```sh
   mingw32-make
   ```

This will generate `omp_delphi.dll` in the project directory.

---

## Using the Library in Delphi

### 1. Setup
- Place `omp_delphi.dll` in your Delphi project directory or in a directory included in your system's PATH.

### 2. Using the Wrapper in a Delphi Project
- Copy `OMPMatrix.pas` from the `wrappers/` folder into your Delphi project.
- In `OMPMatrix.pas`, set the `DLLPATH` constant to the correct location of `omp_delphi.dll`:
  ```delphi
  const
    DLLPATH = 'path_to_your_dll\omp_delphi.dll';
  ```
- Use the functions declared in the wrapper to perform matrix operations. Example:

  ```delphi
  procedure TForm1.Button1Click(Sender: TObject);
  var
    A, B, C: array of Double;
    M, K, N, Threads: Integer;
  begin
    M := 4; K := 4; N := 4;
    Threads := 4;
    SetLength(A, M * K);
    SetLength(B, K * N);
    SetLength(C, M * N);
    // Initialize matrices A and B
    A[0] := 1; A[1] := 2; A[2] := 3; A[3] := 4;
    B[0] := 5; B[1] := 6; B[2] := 7; B[3] := 8;
    // Call the parallel Gustavson multiplication function
    MMParGustavson(@A[0], @B[0], @C[0], M, K, N, Threads);
    ShowMessage(Format('C[0]: %f', [C[0]]));
  end;
  ```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact

For questions, suggestions, or support, please open an issue or contact the maintainer.