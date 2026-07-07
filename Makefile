# DLLifting — DL/DP hybrid lifting for knapsack cover inequalities
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11 -Iinclude -fPIC
LDFLAGS = -shared
LIBS = -lm
RPATH = -Wl,-rpath,'$$ORIGIN'

REDUCTION ?= 1
ifeq ($(REDUCTION),1)
  CXXFLAGS += -DDLLIFTING_REDUCTION
endif

SRC = src/DLLifting.cpp src/dllifting_c.cpp
HDR = include/DLLifting.h
LIB = libdllifting.so
OBJ = DLLifting.o dllifting_c.o

.PHONY: all test test-all examples example install clean

all: $(LIB)

$(LIB): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

DLLifting.o: src/DLLifting.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c src/DLLifting.cpp -o $@

dllifting_c.o: src/dllifting_c.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c src/dllifting_c.cpp -o $@

test_dllifting: tests/test_dllifting.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_dllifting.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

test_isgeq: tests/test_isgeq.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_isgeq.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

test_mixed_vars: tests/test_mixed_vars.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_mixed_vars.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

test_mixed_vars_r: tests/test_mixed_vars.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -DDLLIFTING_REDUCTION tests/test_mixed_vars.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

example: examples/example.cpp $(LIB)
	$(CXX) $(CXXFLAGS) examples/example.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

examples: example

test: test_dllifting
	./test_dllifting

test-all: test test_isgeq test_mixed_vars test_mixed_vars_r
	./test_isgeq
	./tests/run_test_mixed_vars.sh

install: $(LIB) $(HDR)
	install -d $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 $(LIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 644 $(HDR) $(DESTDIR)$(PREFIX)/include/

clean:
	rm -f $(OBJ) $(LIB) test_dllifting test_isgeq test_mixed_vars test_mixed_vars_r example
	rm -f tests/results_mixed_*.txt
