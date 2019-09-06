// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* This will became a unit test in future.
 */

#define TOP_TREE_INTEGRITY_LEVEL 2
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
	int length;

	using Vertex = TopTreeVertex;
	using ClusterType = TopTreeClusterType;
	using Data = PathLengthUserData;

	PathLengthUserData(int length = 1) : length(length) {};

	void join(
			ClusterType type,
			Data &parent, Data &child1, Data &child2,
			Vertex vertex1, Vertex vertex2, Vertex innerVertex) {

		assert(child1.length > 0);
		if (type == ClusterType::COMPRESS) {
			parent.length = child1.length + child2.length;
		} else {
			parent.length = child1.length;
		}

		// deliberately destroy data in child1 to test split;
		// children modifications in join may be forbidden in future versions
		child1.length = 0;
	}

	void split(
			ClusterType type,
			Data &parent, Data &child1, Data &child2,
			Vertex vertex1, Vertex vertex2, Vertex innerVertex) {

		assert(child1.length == 0);
		// repair data destroyed by join
		if (type == ClusterType::COMPRESS) {
			child1.length = parent.length - child2.length;
		} else {
			child1.length = parent.length;
		}
	}
};

class TestingTopTree : public TopTree<PathLengthUserData, PathLengthUserData> {
	using Vertex = TopTreeVertex;
	using ClusterType = TopTreeClusterType;
	using Node = typename TopTree<PathLengthUserData, PathLengthUserData>::Node;

	public:

		APSP paths;
		APSP edges;

		// two random perfectly ballaced trees
		TestingTopTree(unsigned int seed, unsigned int depth1, unsigned int depth2, float rakeToAllRatio = 0.5) :
			rndGen(seed), rndDist(0,1), edges((1 << depth1) + (1 << depth2) + 2), paths((1 << depth1) + (1 << depth2) + 2)
		{
			Node *root = createSubTree(depth1, rakeToAllRatio, this->newVertex());
			this->markNodeAsRoot(root);
			assert(Integrity::treeConsistency(root));
			this->validateTree(root);
			root = createSubTree(depth2, rakeToAllRatio, this->newVertex());
			this->markNodeAsRoot(root);
			assert(Integrity::treeConsistency(root));
			this->validateTree(root);
			paths.compute();
			printf("Base nodes: %d\nRake nodes: %d\nCompress nodes: %d\n", baseCnt, rakeCnt, compressCnt);
		}


		virtual void link(Vertex u, Vertex v, PathLengthUserData data, PathLengthUserData data2) {
			assert(0);
		}
		virtual std::tuple<PathLengthUserData, PathLengthUserData> cut(Vertex u, Vertex v) {
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
				auto &[dataPathLen, dataEdgesCnt] = node->userData;
				dataPathLen.length = rndDist(rndGen) * 1000 + 1;
				node->setBoundary(sharedVertex, otherVertex);
				baseCnt++;
				edges.newEdge(sharedVertex, otherVertex, dataPathLen.length);
				paths.newEdge(sharedVertex, otherVertex, dataPathLen.length);
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
	TestingTopTree tree(42, 7, 6);

	std::cout << "Testing expose (connectivity and lengths of all paths):" << std::endl;
	for (TopTreeVertex u = 0; u < tree.getVerticesCnt(); u++) {
		std::cout << "  ...from vertex " << u << std::endl;
		for (TopTreeVertex v = u + 1; v < tree.getVerticesCnt(); v++) {
			if (tree.expose(u, v)) {
				auto [u2,v2] = tree.getBoundary();
				auto [rootData, rootData2] = tree.getRootData();
				assert(((u == u2) && (v == v2)) || ((u == v2) && (v == u2)));
				assert(rootData.length == tree.paths.getLength(u,v));
			} else {
				assert(!tree.paths.isConnected(u, v));
			}
		}
	}

	std::cout << "Testing path search (enumeration of chosen paths):" << std::endl;
	for (TopTreeVertex u = 0; u < tree.getVerticesCnt(); u+=21) {
		std::cout << "  ...from vertex " << u << std::endl;
		for (TopTreeVertex v = u + 1; v < tree.getVerticesCnt(); v+=5) {
			std::cout << "    ...to vertex " << v << std::endl;
			if (tree.expose(u, v)) {
				TopTreeVertex currVert = u;
				for (int i = 1; currVert != v; i++) {
					tree.pathSearch<1>([u, i](auto type, auto parent, auto child1, auto child2, auto vert1, auto vert2, auto innerVert) {
						if (vert1 == u) {
							return parent.length - child2.length < i;
						} else {
							return child2.length >= i;
						}
					});
					auto [boundVert1, boundVert2] = tree.getBoundary();
					assert((boundVert1 == currVert) || (boundVert2 == currVert));
					TopTreeVertex nextVert = boundVert1 == currVert ? boundVert2 : boundVert1;
					assert(tree.edges.isConnected(currVert, nextVert));
					currVert = nextVert;
					tree.expose(u, v);
				}
			}
		}
	}

	// TODO: test non-path search
	// TODO: implement and test deactivating data propagation


	return 0;
}
