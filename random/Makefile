# Program name
NAME=exec

# Files
CPPFILES=main.cpp

# Tools
CC=g++
LD=g++

# Flags/Options
CFLAGS=-std=c++11 -O3 -funroll-loops -ftree-vectorize -Wno-deprecated
LDFLAGS=-fopenmp -Wl,-R -Wl,build
LIBS=

INCLUDE_DIR=
LIBRARY_DIR=

all: clean cc ld
	@echo "\nDONE!!"

cc:
	@echo "\nCOMPILING..." $(CPPFILES)
	$(CC) $(CFLAGS) $(foreach d, $(INCLUDE_DIR), -I$d) -c $(CPPFILES)

ld:
	@echo "\nLINKING..."
	$(LD) $(foreach d, $(LIBRARY_DIR), -L$d) $(LDFLAGS) -o $(NAME) *.o $(LIBS)

clean:
	@echo "\nDELETING ALL FILES !!"
	rm -f *.o $(NAME)
