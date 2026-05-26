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

## Requirements & Development Environment

The library is developed, tested, and intended to be used in the following environment:

- **Operating System**: Windows 11 (64-bit)
- **Compiler**: GCC 15.1.0
- **Toolchain**: MinGW-w64 12.0.0 UCRT (release 1) from [winlibs.com](https://winlibs.com/)
- **Dependencies**: OpenMP for multi-threading.
- **Build System**: `mingw32-make`
- **Delphi Integration**: RAD Studio (Delphi 12.1 Community Edition)

### Toolchain Setup and Installation Steps
To compile this library, you must install the specific GCC toolchain that supports OpenMP multi-threading. Follow these steps carefully:

**1. Download the Correct Compiler Version**
* Navigate to [winlibs.com](https://winlibs.com/).

* Locate the block titled: **"GCC 15.1.0 (with POSIX threads) + MinGW-w64 12.0.0 (UCRT) - release 1"**.

* Under this block, look for the **Win64** row (specifically: *Win64 without LLVM/Clang/LLD/LLDB*).

* Click on **7-Zip archive** (recommended for smaller download size) or **Zip archive** to download the package.

   > **Note:** It is strictly required to choose the **POSIX threads** variant because GCC's OpenMP implementation (`libgomp`) relies natively on the POSIX threading standard (`pthread`) for stable multi-threaded execution.

**2. Extraction and Directory Setup**
* Extract the downloaded archive.

* Move the extracted folder (usually named `mingw64`) to a clean, global system path without spaces (e.g.: `C:\winlibs`).

**3. Add to Windows System Path**  
To call `gcc` and `mingw32-make` from any terminal, add the toolchain's `bin` directory to your system environment variables:

* Press `Win + S`, type **"env"**, and select **Edit the system environment variables**.

* In the System Properties window, click on the **Environment Variables...** button.

* Under **System variables** (bottom list), select the variable named **Path** and click **Edit...**.

* Click **New** on the right side and type the absolute path to the extracted `bin` directory:
   ```
   C:\winlibs\bin
   ```
* Click **OK** to close each of the three setup windows

**4. Verification**  
Open a new Command Prompt (`cmd`) or PowerShell window and run the following commands to verify the installation:
```sh
gcc --version
mingw32-make --version
```

If configured correctly, both commands will output their respective version information (GCC 15.1.0 and GNU Make), confirming that your development environment is ready.

## Project Structure

```text
offc_delphi/
│
├── build/                  # Compiled binaries and intermediate files
├── include/                # Public API headers (C)
│   └── offc_delphi.h
├── interface/              # Delphi wrapper
│   └── OffC.pas
├── src/                    # C source code
│   └── offc_delphi.c
├── LICENSE                 # License information
├── Makefile                # Build script for DLL
└── README.md               # Project documentation
```

## Compilation Instructions

Once your toolchain environment is properly set up, building the library is simple:
1. Open a terminal (Command Prompt or PowerShell).

2. Navigate to the root directory of this repository:
   ```sh
   cd path\to\offc_delphi
   ```

3. Run the compiler build tool:
   ```sh
   mingw32-make
   ```

This will generate `offc_delphi.dll` in the `build\` directory.

## Using the Library in Delphi

### 1. Linking the DLL (Choose One Method)
To allow your Delphi application to load and communicate with the dynamic library, you must make `offc_delphi.dll` discoverable by Windows. You can achieve this using one of the following methods:

**Method A: Add the Build Directory to the System Path (*Recommended for Development*)**  
During active development of the C library, you might need to recompile the DLL frequently. Copying the file manually after every compilation is tedious and error-prone.

1. Add the absolute path of your local `offc_delphi\build\` directory directly to your Windows system's PATH environment variable (following the same steps used for the compiler setup).

**Method B: Place the DLL in the Executable Folder (*Recommended for Distribution*)**  
If you are packaging, deploying, or sharing your compiled Delphi application:

1. Locate your Delphi project's output directory (usually named `Win64\Debug` or `Win64\Release` inside your project's folder structure).

2. Copy `offc_delphi.dll` from `offc_delphi\build\` and paste it directly into that folder, alongside your compiled `.exe` file.

### 2. Integrating the Wrapper in RAD Studio / Delphi
To use the library functions natively inside your Delphi Pascal source code, follow these steps:

1. Open your Delphi project in RAD Studio.

2. Go to the top menu and select **Project > Options...**.

3. In the left panel, navigate to **Building > Delphi Compiler**.

4. In the main area, click the **Search Path** field.

5. Click on the ellipsis `...` button and add the absolute path to the `interface\` folder inside your repository (e.g.: `C:\Projects\offc_delphi\interface`).

6. Click **OK** and then **Save**.

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