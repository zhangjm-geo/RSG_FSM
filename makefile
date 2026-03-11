
CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra
TARGET = eikonal_fsm

SRC = Eikonal_Traveltime_RSGFSM.cpp
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



