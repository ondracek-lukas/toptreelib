#!/bin/bash

.PHONY: all TopTreeUnitTest

all: TopTreeUnitTest

TopTreeUnitTest: TopTreeUnitTest.out
	ulimit -c0; time ./TopTreeUnitTest.out

TopTreeUnitTest.out: |../TopTree.hpp ../TopTreeInternals/*.hpp TopTreeIntegrity.hpp

%.out: %.cpp
	g++ -std=c++17 $^ -o $@ -fdiagnostics-color=always -g
