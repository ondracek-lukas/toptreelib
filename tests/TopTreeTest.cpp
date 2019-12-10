// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#define TOP_TREE_INTEGRITY
#include "../BiasedTreeTopTree.hpp"

#include <cassert>
#include <random>
#include <functional>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <limits>

std::default_random_engine rndGen(42);
std::uniform_real_distribution<float> rndDist;

class ShadowTree {
	struct Edge;
	struct Vertex {
		size_t index;
		std::list<Edge *> edges;
		Edge *rotate() {
			Edge *e = edges.front();
			if (e) {
				edges.pop_front();
				edges.push_back(e);
			}
			return e;
		}
		Edge *getNext(Edge *e) {
			if (edges.back() != e) {
				while (rotate() != e) {}
			}
			return rotate();
		}
		Vertex(size_t index) : index(index) {};
	};

	struct Edge {
		Vertex *vertices[2];
		int length;
		Edge (Vertex *u, Vertex *v, int length) : vertices{u, v}, length(length) {};
		Vertex *getOther(Vertex *v) {
			return vertices[vertices[0] == v];
		}
		bool hasVertices(Vertex *u, Vertex *v) {
			return ((vertices[0] == u) && (vertices[1] == v)) ||
			       ((vertices[1] == u) && (vertices[0] == v));
		}
	};

	std::vector<Edge *> edges;
	std::vector<Vertex *> vertices;

	Edge *findEdge(Vertex *u, Vertex *v) {
		for (auto e : u->edges) {
			if (e->hasVertices(u, v)) {
				return e;
			}
		}
		return nullptr;
	}

	bool link(Vertex *u, Vertex *v, int length) {
		if (findEdge(u, v)) return false;
		Edge *e = new Edge(u, v, length);
		u->edges.push_back(e);
		v->edges.push_back(e);
		edges.push_back(e);
		return true;
	}

	bool cut(Vertex *u, Vertex *v) {
		Edge *edge = findEdge(u, v);
		if (!edge) return false;
		edges.erase(std::find(edges.begin(), edges.end(), edge));
		for (auto vertex : {u, v}) {
			vertex->edges.erase(std::find(vertex->edges.begin(), vertex->edges.end(), edge));
		}
		delete edge;
		return true;
	}


	int pathLength(Vertex *u, Vertex *v) {
		int length = -1;
		int l = 0;
		Edge *uFirst = u->rotate();
		if (!uFirst) return -1;
		Vertex *w = u;
		Edge *e = uFirst;

		do {
			l += e->length;
			e->length *= -1;
			w = e->getOther(w);
			if (w == v) {
				length = l;
			}
			e = w->getNext(e);
		} while ((w != u) || (e != uFirst));

		assert(l == 0);
		return length;
	}

	public:
		TopTreeVertex newVertex() {
			Vertex *v = new Vertex(vertices.size());
			vertices.push_back(v);
			return v->index;
		}

		bool link(TopTreeVertex u, TopTreeVertex v, int length) {
			return link(vertices[u], vertices[v], length);
		}

		bool cut(TopTreeVertex u, TopTreeVertex v) {
			return cut(vertices[u], vertices[v]);
		}

		int pathLength(TopTreeVertex u, TopTreeVertex v) {
			return pathLength(vertices[u], vertices[v]);
		}

		size_t edgesCnt() {
			return edges.size();
		}

		size_t verticesCnt() {
			return vertices.size();
		}

		std::tuple<TopTreeVertex, TopTreeVertex, int> getEdge(size_t index) {
			return {edges[index]->vertices[0]->index, edges[index]->vertices[1]->index, edges[index]->length};
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


int main() {
	using Vertex = TopTreeVertex;
	ShadowTree shadowTree;
	BiasedTreeTopTree<PathLengthUserData> forestBTTT;
	TopTree<PathLengthUserData> &forest = forestBTTT;
	
	Vertex u = forest.newVertex();
	assert(shadowTree.newVertex() == u);
	int step = 0;

	// create one component of given number of edges
	for (int i = 0; i < 1100; i++) {
		if (rndDist(rndGen) < 0.2) {
			u = rndDist(rndGen) * (i + 1);
		}
		Vertex v = forest.newVertex();
		assert(shadowTree.newVertex() == v);
		std::cout << ++step << ": link(" << u  << ", " << v << ")" << std::endl;
		int edgeLength = rndDist(rndGen) * 1000;
		forest.link(u, v, edgeLength);
		assert(shadowTree.link(u, v, edgeLength));
		u = v;
	}

	for (int i = 0; i < 100; i++) {
		while (shadowTree.verticesCnt() - shadowTree.edgesCnt() < 31) {
			// cut random edge
			size_t edgeIndex = rndDist(rndGen) * shadowTree.edgesCnt();
			auto [u, v, length] = shadowTree.getEdge(edgeIndex);
			std::cout << ++step << ": cut(" << u  << ", " << v << ")" << std::endl;
			auto [data] = forest.cut(u, v);
			assert(data.length == length);
			assert(shadowTree.cut(u, v));
		}

		for (int j = 0; j < 10; j++) {
			// choose random pair of vertices; either link or expose them
			Vertex u, v;
			u = rndDist(rndGen) * shadowTree.verticesCnt();
			do {
				v = rndDist(rndGen) * shadowTree.verticesCnt();
			} while (u == v);
			int length = shadowTree.pathLength(u, v);
			if (length < 0) {
				std::cout << ++step << ": link(" << u  << ", " << v << ")" << std::endl;
				int edgeLength = rndDist(rndGen) * 1000;
				forest.link(u, v, edgeLength);
				assert(shadowTree.link(u, v, edgeLength));
			} else {
				std::cout << ++step << ": expose(" << u << ", " << v << ")" << std::endl;
				forest.expose(u, v);
				auto data = forest.getRootData();
				assert(data.length == length);
			}
		}
	}


	return 0;
}
