#===============================================================
# Makefile for MyTunesNew
# Builds MTconvert from MTconvert.c, MTpkg.c, and PSpkg.c
#===============================================================

# --- Compiler and flags ---
CC      := clang
CFLAGS  := -O2 -std=c99 -Wall -Wextra -pedantic
LDFLAGS := -lm

# --- Source files and targets ---
SRC     := MTconvert.c MTpkg.c PSpkg.c
OBJ     := $(SRC:.c=.o)
TARGET  := MTconvert

# --- Default target ---
all: $(TARGET)

# --- Link the final executable ---
$(TARGET): $(OBJ)
	@echo "Linking $(TARGET)..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# --- Compile each .c into .o ---
%.o: %.c
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# --- Clean up build artifacts ---
clean:
	@echo "Cleaning up..."
	rm -f $(OBJ) $(TARGET)

# --- Force rebuild ---
rebuild: clean all

# --- Phony targets ---
.PHONY: all clean rebuild
