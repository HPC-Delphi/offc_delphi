# Variables
CC = gcc
CFLAGS = -Wall -shared -O3 -lm -fopenmp
TARGET = openmp_matrix_lib.dll
IMPLIB = openmp_matrix_lib.a
SRC = src/strassen_utils.c src/openmp_matrix_lib.c
OBJ = strassen_utils.o openmp_matrix_lib.o
HEADERS = include/strassen_utils.h include/openmp_matrix_lib.h

# Rules
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) -Wl,--out-implib,$(IMPLIB)

$(OBJ): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -c $(SRC)

clean:
	del /Q $(OBJ) $(TARGET) $(IMPLIB)