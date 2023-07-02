# Define the directories
DIR1 :=.
DIR2 :=../../Source

INCLUDE_DIR := -I../../Source 
INCLUDE_DIR += -I./


# Define the source files
SRC1 := $(wildcard $(DIR1)/*.cpp)
SRC2 := $(wildcard $(DIR2)/*.cpp)

# Define the object files
OBJ1 := $(patsubst $(DIR1)/%.cpp, $(DIR1)/%.o, $(SRC1))
OBJ2 := $(patsubst $(DIR2)/%.cpp, $(DIR2)/%.o, $(SRC2))

# Define the header files
HEADER1 := $(wildcard $(DIR1)/*.H)
HEADER2 := $(wildcard $(DIR2)/*.H)

# Define the compiler and flags
CC := mpicxx
CFLAGS := -Wall -Wextra -std=c++14 $(INCLUDE_DIR) 

# Define the target executable
TARGET := riemann1d.ex

all: $(TARGET)

$(TARGET): $(OBJ1) $(OBJ2)
	$(CC) $(CFLAGS) $^ -o $@

$(DIR1)/%.o: $(DIR1)/%.cpp $(HEADER1)
	$(CC) $(CFLAGS) -c $< -o $@

$(DIR2)/%.o: $(DIR2)/%.cpp $(HEADER2)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(DIR1)/*.o
	rm -f $(DIR2)/*.o
	rm -f $(TARGET)
	
