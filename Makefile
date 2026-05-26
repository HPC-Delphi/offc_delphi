CC        := gcc
CFLAGS    := -Wall -shared -O3 -lm -fopenmp -mavx -mavx2 -mfma -march=native -mtune=native

SRC_DIR   := src
OBJ_DIR   := build
INC_DIR   := include

TARGET    := $(OBJ_DIR)\offc_delphi.dll
SRC       := $(SRC_DIR)\offc_delphi.c
OBJ       := $(OBJ_DIR)\offc_delphi.o
INC       := $(INC_DIR)\offc_delphi.h

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)
	del /q $(OBJ)

$(OBJ): $(SRC) $(INC) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)

clean:
	rmdir /s /Q $(OBJ_DIR)
