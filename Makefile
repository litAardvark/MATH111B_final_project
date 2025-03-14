# Compiler
CXX = g++

# Compiler Flags
CXXFLAGS = -Wall -g -std=c++11

# Linker Flags
LDFLAGS = -Wl,--no-relax

# Source Files
SRCS = qmp_network.cpp
SRC2 = visualize_graph.cpp
# Object Files
OBJS = $(SRCS:.cpp=.o)
OBJ2 = $(SRC2:.cpp=.o)
# Executable Name
EXEC = qmp
EXEC2 = visgraph
# Default Target
all: $(EXEC)

# Target Viz
viz: $(EXEC2)

# Link the object files to create the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

$(EXEC2): $(OBJ2) 
	$(CXX) $(CXXFLAGS) -o $(EXEC2) $(OBJ2) $(LDFLAGS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target to remove generated files
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets (targets that don't represent files)
.PHONY: all clean
