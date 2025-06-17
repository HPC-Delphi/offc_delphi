CC        := gcc
CFLAGS    := -Wall -shared -O3 -lm -fopenmp -mavx -mavx2 -mfma -march=native -mtune=native

SRC_DIR   := src
OBJ_DIR   := build
INC_DIR   := include

TARGET    := $(OBJ_DIR)\bridge_delphi.dll
SRC   := $(SRC_DIR)\bridge_delphi.c
OBJ   := $(OBJ_DIR)\bridge_delphi.o
HEADERS      := $(INC_DIR)\bridge_delphi.h

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

$(OBJ): $(SRC) $(HEADERS) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)

clean:
	rmdir /s /Q $(OBJ_DIR)
