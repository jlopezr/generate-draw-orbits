# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra

# Linker flags
LDFLAGS = -lm

# Target executables
TARGETS = satellite kepler

# Source files
SAT_SRCS = satellite.c
KEP_SRCS = kepler.c

# Object files
SAT_OBJS = $(SAT_SRCS:.c=.o)
KEP_OBJS = $(KEP_SRCS:.c=.o)

# Default target
all: $(TARGETS)

# Link the target executables
satellite: $(SAT_OBJS)
	$(CC) $(CFLAGS) -o satellite $(SAT_OBJS) $(LDFLAGS)

kepler: $(KEP_OBJS)
	$(CC) $(CFLAGS) -o kepler $(KEP_OBJS) $(LDFLAGS)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run the satellite program and pipe output to Python script
run-satellite: satellite
	./satellite | python3 plot_orbit.py

# Run the kepler program and pipe output to Python script
run-kepler: kepler
	./kepler | python3 plot_orbit3d.py

# Clean up build files
clean:
	rm -f $(TARGETS) $(SAT_OBJS) $(KEP_OBJS)

# Phony targets
.PHONY: all clean run-satellite run-kepler