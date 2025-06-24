# offc_delphi - High Performance Matrix Operations Library

[![Platform: Windows](https://img.shields.io/badge/Platform-Windows%2011-blue)](https://www.microsoft.com/en-us/windows/windows-11)
[![GCC](https://img.shields.io/badge/GCC-15.1.0-brightgreen)](https://winlibs.com/)

`offc_delphi` is a high-performance dynamic library written in C, providing efficient matrix operations. Designed for seamless integration with Delphi applications, it enables advanced matrix computations, including sequential implementations of Gustavson algorithms. The library is distributed as a DLL for easy cross-language usage.

---

## Features

- **Matrix Multiplication**: Optimized matrix multiplication for high performance.
- **Delphi Integration**: Includes a Delphi wrapper for straightforward use in Delphi projects.
- **Advanced Algorithms**: Implements Gustavson algorithm (sequential variant).
- **Cross-Language Compatibility**: Distributed as a DLL, accessible from various programming environments.
- **Modern Toolchain Support**: Built and tested with the latest GCC and MinGW-w64 toolchains for Windows 11.

---

## Requirements & Development Environment

The library is developed, tested, and intended to be used in the following environment:

- **Operating System**: Windows 11 (64-bit)
- **Compiler**: GCC 15.1.0
- **Toolchain**: MinGW-w64 12.0.0 UCRT (release 1) from [winlibs.com](https://winlibs.com/)
- **Build System**: `mingw32-make`
- **Delphi Integration**: RAD Studio (Delphi 12.1 Community Edition)

---

## Project Structure

```
offc_delphi/
│
├── build/                  # Compiled binaries and intermediate files
├── include/                # Public API and utility headers
│   ├── offc_delphi.h
├── src/                    # C source code
│   ├── offc_delphi.c
├── wrappers/               # Language wrappers (Delphi)
│   └── OffC.pas
├── LICENSE                 # License information (MIT)
├── Makefile                # Build script for DLL
└── README.md               # Project documentation
```

---

## Compilation Instructions

To build the library as a DLL:

1. Download and extract the [winlibs.com](https://winlibs.com/) GCC toolchain (ensure it is in your PATH).
2. Open a Command Prompt or PowerShell window.
3. Navigate to the project root directory (`offc_delphi`).
4. Run:

   ```sh
   mingw32-make
   ```

This will generate `offc_delphi.dll` in the project directory.

---

## Using the Library in Delphi

### 1. Setup
- Place `offc_delphi.dll` in the directory where the executable is created upon compiling your Delphi project, or ensure it is in a directory included in your system's PATH.

### 2. Using the Wrapper in a Delphi Project
- Add the path to `OffC.pas` in the "Search Path" under Project Options in your Delphi project, or copy the file into your Delphi project.
- Use the functions declared in the wrapper to perform matrix operations.