#!/bin/bash

.PHONY: all %.test clean
.PRECIOUS: %.out

all: InnerListTest.test TopTreeBaseTest.test BiasedTreeTopTreeTest.test ContractionTopTreeTest.test

clean:
	rm -f *.out

%.test: %.out
	#time ./$<
	ulimit -c0; time ./$<

InnerListTest.out: ../TopTreeInternals/InnerList.hpp
TopTreeBaseTest.out: ../TopTree.hpp ../TopTreeInternals/*.hpp TopTreeIntegrity.hpp
%TopTreeTest.out: TopTreeTest.cpp ../TopTree.hpp ../%TopTree.hpp ../TopTreeInternals/*.hpp TopTreeIntegrity.hpp %TopTreeIntegrity.hpp
	g++ -std=c++17 $< -o $@ -fdiagnostics-color=always -g -DTOP_TREE_TYPE=$*TopTree

%.out: %.cpp
	g++ -std=c++17 $< -o $@ -fdiagnostics-color=always -g
