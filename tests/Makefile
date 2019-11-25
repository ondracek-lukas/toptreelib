#!/bin/bash

.PHONY: all %.test
.PRECIOUS: %.out

all: InnerListTest.test TopTreeBaseTest.test TopTreeTest.test

%.test: %.out
	ulimit -c0; time ./$<

InnerListTest.out: ../TopTreeInternals/InnerList.hpp
TopTreeBaseTest.out: ../TopTree.hpp ../TopTreeInternals/*.hpp TopTreeIntegrity.hpp
TopTreeTest.out: ../TopTree.hpp ../BiasedTreeTopTree.hpp ../TopTreeInternals/*.hpp TopTreeIntegrity.hpp BiasedTreeTopTreeIntegrity.hpp

%.out: %.cpp
	g++ -std=c++17 $< -o $@ -fdiagnostics-color=always -g
