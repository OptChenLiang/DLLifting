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

.PHONY: all test examples example install clean

all: $(LIB)

$(LIB): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

DLLifting.o: src/DLLifting.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c src/DLLifting.cpp -o $@

dllifting_c.o: src/dllifting_c.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c src/dllifting_c.cpp -o $@

test_dllifting: tests/test_dllifting.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_dllifting.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

example: examples/example.cpp $(LIB)
	$(CXX) $(CXXFLAGS) examples/example.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

examples: example

test: test_dllifting
	./test_dllifting

install: $(LIB) $(HDR)
	install -d $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 $(LIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 644 $(HDR) $(DESTDIR)$(PREFIX)/include/

clean:
	rm -f $(OBJ) $(LIB) test_dllifting example
