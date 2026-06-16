# Stand-alone tests for DLLifting (no external solver required)
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -Iinclude -DNDEBUG
LIBS = -lm

SRC = src/DLLifting.cpp src/dllifting_c.cpp

test_dllifting: tests/test_dllifting.cpp $(SRC) include/DLLifting.h
	$(CXX) $(CXXFLAGS) -DDLLIFTING_REDUCTION tests/test_dllifting.cpp $(SRC) -o $@ $(LIBS)

test_isgeq: tests/test_isgeq.cpp $(SRC) include/DLLifting.h
	$(CXX) $(CXXFLAGS) -DDLLIFTING_REDUCTION tests/test_isgeq.cpp $(SRC) -o $@

test_mixed_vars: tests/test_mixed_vars.cpp $(SRC) include/DLLifting.h
	$(CXX) $(CXXFLAGS) tests/test_mixed_vars.cpp $(SRC) -o $@ $(LIBS)

test_mixed_vars_r: tests/test_mixed_vars.cpp $(SRC) include/DLLifting.h
	$(CXX) $(CXXFLAGS) -DDLLIFTING_REDUCTION tests/test_mixed_vars.cpp $(SRC) -o $@ $(LIBS)

run-mixed: test_mixed_vars test_mixed_vars_r
	./tests/run_test_mixed_vars.sh

run: test_dllifting
	./test_dllifting

run-all: test_dllifting test_isgeq
	./test_dllifting
	./test_isgeq

clean:
	rm -f test_dllifting test_isgeq test_mixed_vars test_mixed_vars_r
	rm -f tests/results_mixed_*.txt

.PHONY: run run-all run-mixed clean
