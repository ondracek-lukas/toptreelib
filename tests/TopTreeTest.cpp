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

std::default_random_engine rndGen(42);
std::uniform_real_distribution<float> rndDist;


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
	using Vertex = TopTreeVertex;

	BiasedTreeTopTree<PathLengthUserData> forestBTTT;
	TopTree<PathLengthUserData> &forest = forestBTTT;
	
	Vertex u = forest.newVertex();

	for (int i = 0; i < 256; i++) {
		if (rndDist(rndGen) < 0.2) {
			u = rndDist(rndGen) * (i + 1);
		}
		Vertex v = forest.newVertex();
		std::cout << std::endl << "=== " << i << ": link(" << u  << ", " << v << ") ===" << std::endl;
		forest.link(u, v, rndDist(rndGen) * 1000);
		u = v;
	}

	/*
	for (int i = 0; i < 4; i++) {
		Vertex v = forest.newVertex();
		forest.link(u, v, i+1);
		u = v;
	}
	u = 2;
	for (int i = 0; i < 3; i++) {
		Vertex v = forest.newVertex();
		forest.link(u, v, i+1);
		u = v;
	}
	std::cout << forest.getRootData().length << std::endl;

	forest.expose(4,7);
	std::cout << forest.getRootData().length << std::endl;
	*/

	return 0;
}
