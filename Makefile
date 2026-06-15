
# Stand-alone tests for DLLifting (no CPLEX required)
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -DDLLIFTING -DNDEBUG -DDLLIFTING_REDUCTION
LIBS = -lm

test_dllifting: test_dllifting.cpp DLLifting.cpp DLLifting.h
	$(CXX) $(CXXFLAGS) test_dllifting.cpp DLLifting.cpp -o $@ $(LIBS)

test_isgeq: test_isgeq.cpp DLLifting.cpp DLLifting.h
	$(CXX) $(CXXFLAGS) test_isgeq.cpp DLLifting.cpp -o $@ $(LIBS)

test_mixed_vars: test_mixed_vars.cpp DLLifting.cpp DLLifting.h
	$(CXX) -Wall -Wextra -O2 -DDLLIFTING -DNDEBUG test_mixed_vars.cpp DLLifting.cpp -o $@ $(LIBS)

test_mixed_vars_r: test_mixed_vars.cpp DLLifting.cpp DLLifting.h
	$(CXX) $(CXXFLAGS) test_mixed_vars.cpp DLLifting.cpp -o $@ $(LIBS)

run-mixed: test_mixed_vars test_mixed_vars_r
	./run_test_mixed_vars.sh

run: test_dllifting
	./test_dllifting

run-reduction: test_dllifting_reduction
	./test_dllifting_reduction

test_dllifting_reduction: test_dllifting.cpp DLLifting.cpp DLLifting.h
	$(CXX) $(CXXFLAGS) -DREDUCTION test_dllifting.cpp DLLifting.cpp -o $@ $(LIBS)

run-all: test_dllifting test_isgeq
	./test_dllifting
	./test_isgeq

clean:
	rm -f test_dllifting test_dllifting_reduction test_isgeq test_mixed_vars test_mixed_vars_r results_mixed_*.txt

.PHONY: run run-reduction run-all run-mixed clean
