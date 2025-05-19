# Compiler and flags
CC        := gcc
CFLAGS    := -Wall -shared -O3 -lm -fopenmp

# Directories
SRC_DIR   := src
OBJ_DIR   := build
INC_DIR   := include

# Targets and artifacts
TARGET    := $(OBJ_DIR)\omp_delphi.dll

# Source and object files
STRASSEN_SRC := $(SRC_DIR)\strassen_utils.c
OMP_SRC      := $(SRC_DIR)\omp_delphi.c
SRCS         := $(STRASSEN_SRC) $(OMP_SRC)

STRASSEN_OBJ := $(OBJ_DIR)\strassen_utils.o
OMP_OBJ      := $(OBJ_DIR)\omp_delphi.o
OBJ          := $(STRASSEN_OBJ) $(OMP_OBJ)

# Header files
HEADERS      := $(INC_DIR)\strassen_utils.h $(INC_DIR)\omp_delphi.h

# Default rule
all: $(TARGET)

# Link shared library
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

# Compile each source into its object
$(OBJ): $(SRCS) $(HEADERS) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $(STRASSEN_SRC) -o $(STRASSEN_OBJ)
	$(CC) $(CFLAGS) -c $(OMP_SRC)      -o $(OMP_OBJ)

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)

# Clean up build artifacts
clean:
	rmdir /s /Q $(OBJ_DIR)
