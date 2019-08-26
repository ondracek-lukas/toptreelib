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

		if (type == TopTreeClusterType::COMPRESS) {
			parent.length = child1.length + child2.length;
		} else {
			parent.length = child1.length;
		}

	}

	void split(
			ClusterType type,
			Data &parent, Data &child1, Data &child2,
			Vertex vertex1, Vertex vertex2, Vertex innerVertex) { };
};

class TestingTopTree : public TopTree<PathLengthUserData> {
	using Vertex = TopTreeVertex;
	using ClusterType = TopTreeClusterType;
	using Node = typename TopTree<PathLengthUserData>::Node;

	/*
	public:
		// fixed small tree
		TestingTopTree() {
			Vertex u = this->newVertex();
			Vertex v = this->newVertex();
			Vertex w = this->newVertex();
			Vertex z = this->newVertex();
			Vertex y = this->newVertex();

			Node *a = this->newNode();
			this->setNodeBoundary(a, BASE, u, v);

			Node *b = this->newNode();
			this->setNodeBoundary(b, BASE, v, w);

			Node *c = this->newNode();
			this->setNodeBoundary(c, BASE, v, z);

			Node *d = this->newNode();
			this->attachSubtree(d, 0, b);
			this->attachSubtree(d, 1, c);
			this->setNodeBoundary(d, RAKE, w, v);

			Node *e = this->newNode();
			this->attachSubtree(e, 0, a);
			this->attachSubtree(e, 1, d);
			this->setNodeBoundary(e, COMPRESS, u, w);
			this->markNodeAsRoot(e);

			this->validateTree(e);
		}
		*/



	public:

		APSP apsp;

		// random perfectly ballaced tree
		TestingTopTree(unsigned int seed, unsigned int depth, float rakeToAllRatio = 0.5) :
			rndGen(seed), rndDist(0,1), apsp((1 << depth) + 1)
		{
			Node *root = createSubTree(depth, rakeToAllRatio, this->newVertex());
			this->markNodeAsRoot(root);
			assert(Integrity::treeConsistency(root));
			this->validateTree(root);
			apsp.compute();
			printf("Base nodes: %d\nRake nodes: %d\nCompress nodes: %d\n", baseCnt, rakeCnt, compressCnt);
		}


		virtual void link(Vertex u, Vertex v, PathLengthUserData data) {
			assert(0);
		}
		virtual std::tuple<PathLengthUserData> cut(Vertex u, Vertex v) {
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
				auto &[data] = node->userData;
				data.length = rndDist(rndGen) * 4 + 1;
				node->setBoundary(sharedVertex, otherVertex);
				baseCnt++;
				apsp.newEdge(sharedVertex, otherVertex, data.length);
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
	//TestingTopTree<> tree;
	//for (int i = 0; i < 5; i++)
	//	for (int j = 0; j < 5; j++)
	//		printf("%d %d %d\n", i, j, tree.sameComponent(i, j));


	/*
	{
		TestingTopTree tree(1,5);
		assert(tree.getVerticesCnt() == tree.apsp.getVerticesCnt());
		std::cout << "APSP matrix:" << std::endl;
		tree.apsp.print();
		auto [u,v] = tree.getBoundary();
		auto [rootData] = tree.getRootData();

		std::cout << "(" << u << ", " << v << ")-path's length by top tree: " << rootData.length << std::endl;
		std::cout << "(" << u << ", " << v << ")-path's length by APSP:     " << tree.apsp.getLength(u,v) << std::endl;
	}
	*/

	for (int i = 10; i < 20; i++) {
		TestingTopTree tree(i, (i < 19 ? 6 : 7));

		auto [u,v] = tree.getBoundary();
		auto [rootData] = tree.getRootData();
		std::cout
			<< "(" << u << "," << v << ") " << tree.apsp.getLength(u,v)
			<< " ~ " << rootData.length << std::endl;
		assert(rootData.length == tree.apsp.getLength(u,v));

		for (TopTreeVertex u = 0; u < tree.getVerticesCnt(); u++) {
			for (TopTreeVertex v = u + 1; v < tree.getVerticesCnt(); v++) {
				tree.expose(u, v);
				auto [u2,v2] = tree.getBoundary();
				auto [rootData] = tree.getRootData();
				std::cout
					<< "(" << u << "," << v << ") " << tree.apsp.getLength(u,v)
					<< " ~ (" << u2 << "," << v2 << ") " << rootData.length << " "
					<< std::endl;

				assert(((u == u2) && (v == v2)) || ((u == v2) && (v == u2)));
				assert(rootData.length == tree.apsp.getLength(u,v));
			}
		}

	}
	return 0;
}
