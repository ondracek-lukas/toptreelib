// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#define TOP_TREE_INTEGRITY
#include "../BiasedTreeTopTree.hpp"

#include <cassert>
#include <random>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>


struct PathLengthUserData {
	using EventData = TopTreeEventData<PathLengthUserData>;

	int length;

	PathLengthUserData(int length = 1) : length(length) {};
	operator int() { return length; }

	static void join(EventData eventData) {
		if (eventData.type == TopTreeClusterType::COMPRESS) {
			eventData.parent.length = eventData.children[0].length + eventData.children[1].length;
		} else {
			eventData.parent.length = eventData.children[0].length;
		}
	}

	static void split(EventData eventData) { }

};


int main() {

	BiasedTreeTopTree<PathLengthUserData> forestBTTT;
	TopTree<PathLengthUserData> &forest = forestBTTT;
	

	return 0;
}
