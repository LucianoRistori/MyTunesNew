#===============================================================
# Makefile for MyTunesNew
# Builds both MTconvert and MTgraph from shared sources
#===============================================================

# --- Compiler and flags ---
CC      := clang
CFLAGS  := -O2 -std=c99 -Wall -Wextra -pedantic
LDFLAGS := -lm

# --- Shared source files ---
COMMON_SRC := MTpkg.c PSpkg.c
COMMON_OBJ := $(COMMON_SRC:.c=.o)

# --- Program-specific sources ---
CONVERT_SRC := MTconvert.c
CONVERT_OBJ := $(CONVERT_SRC:.c=.o) $(COMMON_OBJ)
CONVERT_BIN := MTconvert

GRAPH_SRC   := MTgraph.c
GRAPH_OBJ   := $(GRAPH_SRC:.c=.o) $(COMMON_OBJ)
GRAPH_BIN   := MTgraph

# --- Default target ---
all: $(CONVERT_BIN) $(GRAPH_BIN)

# --- Link MTconvert ---
$(CONVERT_BIN): $(CONVERT_OBJ)
	@echo "Linking $@..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# --- Link MTgraph ---
$(GRAPH_BIN): $(GRAPH_OBJ)
	@echo "Linking $@..."
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# --- Compile each .c into .o ---
%.o: %.c
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# --- Clean up build artifacts ---
clean:
	@echo "Cleaning up..."
	rm -f *.o $(CONVERT_BIN) $(GRAPH_BIN)

# --- Force rebuild ---
rebuild: clean all

# --- Phony targets ---
.PHONY: all clean rebuild
