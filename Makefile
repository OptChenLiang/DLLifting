# DLLifting — DL/DP hybrid lifting for knapsack cover inequalities
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11 -Iinclude -fPIC -DNDEBUG
LDFLAGS = -shared
LIBS = -lm
RPATH = -Wl,-rpath,'$$ORIGIN'
PREFIX ?= $(HOME)/.local

REDUCTION ?= 1
ifeq ($(REDUCTION),1)
  CXXFLAGS += -DDLLIFTING_REDUCTION
endif

SRC = src/DLLifting.cpp
HDR = include/DLLifting.h include/dllifting_c.h
LIB = libdllifting.so
OBJ = DLLifting.o

.PHONY: all test test-all examples example install clean

all: $(LIB)

$(LIB): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

$(OBJ): $(SRC) include/DLLifting.h
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $@

test_dllifting: tests/test_dllifting.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_dllifting.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

example: examples/example.cpp $(LIB)
	$(CXX) $(CXXFLAGS) examples/example.cpp -L. -ldllifting $(RPATH) -o $@ $(LIBS)

examples: example

# Default gate: core (<=) suite only — must pass for release.
test: test_dllifting
	./test_dllifting

# Extended: also run >= and mixed suites (geq may report known failures).
test-all: test_dllifting
	./test_dllifting --all

install: $(LIB) $(HDR)
	install -d $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 $(LIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 644 $(HDR) $(DESTDIR)$(PREFIX)/include/

clean:
	rm -f $(OBJ) $(LIB) test_dllifting example
