# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra

# Linker flags
LDFLAGS = -lm

# Target executable
TARGET = satellite

# Source files
SRCS = satellite.c

# Object files
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Link the target executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run the program and pipe output to Python script
run: $(TARGET)
	./$(TARGET) | python3 plot_orbit.py

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJS)

# Phony targets
.PHONY: all clean
