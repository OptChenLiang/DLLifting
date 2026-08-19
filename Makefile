# DPLifting — DPL/DPT hybrid lifting for general knapsack set
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11 -Iinclude -fPIC -DNDEBUG
LDFLAGS = -shared
LIBS = -lm
RPATH = -Wl,-rpath,'$$ORIGIN'
PREFIX ?= $(HOME)/.local

REDUCTION ?= 1
ifeq ($(REDUCTION),1)
  CXXFLAGS += -DDPLIFTING_REDUCTION
endif

SRC = src/DPLifting.cpp
HDR = include/DPLifting.h include/dplifting_c.h
LIB = libdplifting.so
OBJ = DPLifting.o

.PHONY: all test test-all examples example install clean

all: $(LIB)

$(LIB): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

$(OBJ): $(SRC) include/DPLifting.h
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $@

test_dplifting: tests/test_dplifting.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_dplifting.cpp -L. -ldplifting $(RPATH) -o $@ $(LIBS)

example: examples/example.cpp $(LIB)
	$(CXX) $(CXXFLAGS) examples/example.cpp -L. -ldplifting $(RPATH) -o $@ $(LIBS)

examples: example

# Default gate: core (<=) + geq (>=); both failures affect exit code.
test: test_dplifting
	./test_dplifting

# Also run mixed-variable benchmark suite.
test-all: test_dplifting
	./test_dplifting --all

install: $(LIB) $(HDR)
	install -d $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 $(LIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 644 $(HDR) $(DESTDIR)$(PREFIX)/include/

clean:
	rm -f $(OBJ) $(LIB) test_dplifting example
