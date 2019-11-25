// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#define TOP_TREE_INTEGRITY
#include "../TopTree.hpp"

#include <cassert>
#include <random>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

class APSP { // All Paths Shortest Path
	const int INF = std::numeric_limits<int>::max()/2;
	unsigned int verticesCnt;
	std::vector<int> matrix;
	int &matrixElem (int u, int v) {
		return matrix[u + v*verticesCnt];
	}

	public:

		APSP(unsigned int verticesCnt) : verticesCnt(verticesCnt) {
			matrix.resize(verticesCnt*verticesCnt, INF);
		}

		void newEdge(int u, int v, int length) {
			matrixElem(u, v) = length;
			matrixElem(v, u) = length;
		}
		int getLength(int u, int v) {
			return matrixElem(u, v);
		}
		int isConnected(int u, int v) {
			return matrixElem(u, v) != INF;
		}
		int getVerticesCnt() {
			return verticesCnt;
		}

		void compute() {
			for (int w = 0; w < verticesCnt; w++) {
				for (int u = 0; u < verticesCnt; u++) {
					for (int v = 0; v < verticesCnt; v++) {
						if (matrixElem(u, v) > matrixElem(u, w) + matrixElem(w, v)) {
							matrixElem(u, v) = matrixElem(u, w) + matrixElem(w, v);
						}
					}
				}
			}
		}

		void print() {
			for (int u = 0; u < verticesCnt; u++) {
				for (int v = 0; v < verticesCnt; v++) {
					if (matrixElem(u, v) == INF) {
						std::cout << std::setw(3) << ".";
					} else {
						std::cout << std::setw(3) << matrixElem(u, v);
					}
				}
				std::cout << std::endl;
			}
		}
};

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

struct TreeMinEdgeUserData {
	using EventData = TopTreeEventData<TreeMinEdgeUserData>;

	int minWeight;
	int delta = 0;  // to be added to all weights in descendants

	TreeMinEdgeUserData(int weight = 0) : minWeight(weight) {};
	operator int() { return minWeight; };

	void increaseAllInTree(int delta) {
		this->delta += delta;
		this->minWeight += delta;
	}

	static void join(EventData eventData) {
		eventData.parent.minWeight =
			eventData.children[ eventData.children[0].minWeight > eventData.children[1].minWeight ].minWeight;
		eventData.parent.delta = 0;
	}

	static void split(EventData eventData) {
		for (size_t i : {0, 1}) {
			eventData.children[i].delta += eventData.parent.delta;
			eventData.children[i].minWeight += eventData.parent.delta;
		}
	}

};

class TestingTopTree : public TopTree<PathLengthUserData, PathLengthUserData, TreeMinEdgeUserData> {
	using Vertex = TopTreeVertex;
	using ClusterType = TopTreeClusterType;
	using Node = typename TopTree<PathLengthUserData, PathLengthUserData, TreeMinEdgeUserData>::Node;

	public:

		APSP paths;
		APSP edges;
		struct { int weight = std::numeric_limits<int>::max(); Vertex u; Vertex v; } minima[10];

		// two random perfectly ballaced trees
		TestingTopTree(unsigned int seed, unsigned int depth1, unsigned int depth2, float rakeToAllRatio = 0.5) :
			rndGen(seed), rndDist(0,1), edges((1 << depth1) + (1 << depth2) + 2), paths((1 << depth1) + (1 << depth2) + 2)
		{
			Node *root = createSubTree(depth1, rakeToAllRatio, this->newVertex());
			this->markNodeAsRoot(root);
			assert(Integrity::treeConsistency(root));
			root = createSubTree(depth2, rakeToAllRatio, this->newVertex());
			this->markNodeAsRoot(root);
			assert(Integrity::treeConsistency(root));
			paths.compute();
			printf("Base nodes: %d\nRake nodes: %d\nCompress nodes: %d\n", baseCnt, rakeCnt, compressCnt);
		}


		virtual void link(Vertex u, Vertex v, PathLengthUserData data, PathLengthUserData data2, TreeMinEdgeUserData data3) {
			assert(0);
		}
		virtual std::tuple<PathLengthUserData, PathLengthUserData, TreeMinEdgeUserData> cut(Vertex u, Vertex v) {
			assert(0);
		}

	private:
		std::default_random_engine rndGen;
		std::uniform_real_distribution<float> rndDist;
		unsigned int baseCnt=0, rakeCnt=0, compressCnt=0;

		Node *createSubTree(unsigned int depth, float rakeToAllRatio, Vertex sharedVertex) {
			Node *node = this->newNode();
			if (depth == 0) {
				Vertex otherVertex = this->newVertex();
				node->setBoundary(sharedVertex, otherVertex);
				baseCnt++;

				auto &[pathLen, edgesCnt, treeMinEdge] = node->userData;
				pathLen = rndDist(rndGen) * 1000 + 1;
				edges.newEdge(sharedVertex, otherVertex, pathLen);
				paths.newEdge(sharedVertex, otherVertex, pathLen);

				treeMinEdge = 0;
				while (treeMinEdge == 0) {
					treeMinEdge = rndDist(rndGen) * 100000 + 1;
					for (auto m : minima) {
						if (m.weight == treeMinEdge) treeMinEdge = 0;
					}
				}

				if (minima[std::size(minima) - 1].weight > treeMinEdge) {
					int i = std::size(minima) - 1;
					for (; i > 0; i--) {
						if (minima[i-1].weight < treeMinEdge) break;
						minima[i] = minima[i-1];
					}
					minima[i] = {treeMinEdge, sharedVertex, otherVertex};
				}
			} else {
				ClusterType type;
				if (rndDist(rndGen) < rakeToAllRatio) {
					type = ClusterType::RAKE;
					rakeCnt++;
				} else {
					type = ClusterType::COMPRESS;
					compressCnt++;
				}
				Node *child1 = createSubTree(depth-1, rakeToAllRatio, sharedVertex);
				if ((type == ClusterType::COMPRESS) || (rndDist(rndGen) < 0.5)) {
					sharedVertex = child1->getOtherVertex(sharedVertex);
				}
				Node *child2 = createSubTree(depth-1, rakeToAllRatio, sharedVertex);
				node->attachChildren(type, child1, child2);
				this->vertexToNodeUpdateInner(node);
			}
			return node;
		}
};

int main() {

	int IntegrityCounter = 0; // do full test just once per some number of cycles
	TopTreeIntegrityLevel = 2;

	// 128-edge tree and 64-edge tree
	// on vertices 0-128 and 129-193
	TestingTopTree forest(42, 7, 6);

	std::cout << "Testing expose with bottom-up data aggregation (connectivity and lengths of all paths):" << std::endl;
	for (TopTreeVertex u = 0; u < forest.getVerticesCnt(); u++) {
		if (u % 10 == 0) std::cout << "  ...from vertex " << u << std::endl;
		for (TopTreeVertex v = u + 1; v < forest.getVerticesCnt(); v++) {
			TopTreeIntegrityLevel = IntegrityCounter++ % 21 ? 1 : 2;
			if (forest.expose(u, v)) {
				auto [u2,v2] = forest.getBoundary();
				auto rootData = forest.getRootData<0>();
				assert(((u == u2) && (v == v2)) || ((u == v2) && (v == u2)));
				assert(rootData.length == forest.paths.getLength(u,v));
			} else {
				assert(!forest.paths.isConnected(u, v));
			}
		}
	}

	std::cout << "Testing path search with bottom-up data aggregation (enumeration of chosen paths):" << std::endl;
	for (TopTreeVertex u = 0; u < forest.getVerticesCnt(); u+=5) {
		std::cout << "  ...from vertex " << u << std::endl;
		for (TopTreeVertex v = u + 1; v < forest.getVerticesCnt(); v+=3) {
			//std::cout << "    ...to vertex " << v << std::endl;
			TopTreeIntegrityLevel = IntegrityCounter++ % 97 ? 1 : 2;
			if (forest.expose(u, v)) {
				TopTreeVertex currVert = u;
				for (int i = 1; currVert != v; i++) {
					forest.pathSearch<1>([u, i](auto eventData) {
						bool rev = eventData.boundary[1] == u;
						return rev ^ (eventData.children[rev].length < i);
					});
					auto [boundVert1, boundVert2] = forest.getBoundary();
					assert((boundVert1 == currVert) || (boundVert2 == currVert));
					TopTreeVertex nextVert = boundVert1 == currVert ? boundVert2 : boundVert1;
					assert(forest.edges.isConnected(currVert, nextVert));
					currVert = nextVert;
					forest.expose(u, v);
				}
			}
		}
	}

	std::cout << "Testing tree search with both bottom-up data aggregation and top-down propagation\n(increasing all edges and enumaration of lowest weight edges):" << std::endl;

	TopTreeIntegrityLevel = 2;

	forest.exposeTree(0);
	forest.getRootData<2>().increaseAllInTree(10);
	forest.exposeTree(129);
	forest.getRootData<2>().increaseAllInTree(10);

	for (auto minimum : forest.minima) {
		forest.exposeTree(minimum.u);
		forest.search<2>([](auto eventData) {
			return eventData.children[1] < eventData.children[0];
		});
		auto [u, v] = forest.getBoundary();
		int weight = forest.getRootData<2>().minWeight;

		std::cout << "  " << minimum.u << " " << minimum.v << ": " << minimum.weight
			<< " -- " << u << " " << v << ": " << weight << std::endl;
		assert(((minimum.u == u) && (minimum.v == v)) || ((minimum.u == v) && (minimum.v == u)));
		assert(minimum.weight + 10 == weight);

		forest.getEdgeData<2>() = std::numeric_limits<int>::max();
	}


	return 0;
}
