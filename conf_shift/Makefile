sources = conf_shift.cpp count.cpp
objects = $(sources:%.cpp=%.o)
prog = conf_shift

all: $(prog)
count.o: count.h
conf_shift.o: count.h datafile.h

CXXFLAGS ?= -O2 -Wall -std=c++11
LDFLAGS ?= -static
CPPFLAGS += -I../include
LDLIBS = -lboost_program_options

$(prog): $(objects)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(objects): %.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	-rm -f $(objects) $(prog)
