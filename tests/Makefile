#!/bin/bash

.PHONY: all %.test
.PRECIOUS: %.out

all: InnerListUnitTest.test TopTreeUnitTest.test

%.test: %.out
	ulimit -c0; time ./$<

InnerListUnitTest.out: ../TopTreeInternals/InnerList.hpp
TopTreeUnitTest.out: ../TopTree.hpp ../TopTreeInternals/*.hpp TopTreeIntegrity.hpp

%.out: %.cpp
	g++ -std=c++17 $< -o $@ -fdiagnostics-color=always -g
