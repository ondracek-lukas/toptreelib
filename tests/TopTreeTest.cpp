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

class ShadowForest {
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
	int multiplicator = 1; // should multiply all descendants
	enum {
		HERE_BASE = 0,
		HERE,
		UP,
		DOWN
	} dataLoc = HERE_BASE;  // for data flow correctness verification

	PathLengthUserData(int length = 1) : length(length) {};
	operator int() { return length; }

	void multiplyAll(int coef) {
		length *= coef;
		multiplicator *= coef;
	}

	static void join(EventData eventData) {
		if (eventData.type == TopTreeClusterType::COMPRESS) {
			eventData.parent.length = eventData.children[0].length + eventData.children[1].length;
		} else {
			eventData.parent.length = eventData.children[0].length;
		}
		eventData.parent.multiplicator = 1;

		assert(eventData.children[0].dataLoc <= HERE);
		assert(eventData.children[1].dataLoc <= HERE);

		eventData.children[0].dataLoc = UP;
		eventData.children[1].dataLoc = UP;
		eventData.parent.dataLoc = HERE;
	}

	static void split(EventData eventData) {
		for (int i : {0, 1}) {
			eventData.children[i].multiplyAll(eventData.parent.multiplicator);
		}

		assert(eventData.children[0].dataLoc == UP);
		assert(eventData.children[1].dataLoc == UP);
		assert(eventData.parent.dataLoc == HERE);

		eventData.children[0].dataLoc = HERE;
		eventData.children[1].dataLoc = HERE;
		eventData.parent.dataLoc = DOWN;
	}
};

struct StatsUserData {
	using EventData = TopTreeEventData<StatsUserData>;

	int maxDepth = 1;
	int nodesCnt = 1;

	struct AvgMax {
		int sum = 0;
		int cnt = 0;
		int max = 0;
	};
	static std::vector<AvgMax> depthPerSize;
	static int totalLogs, totalSplits, totalJoins;

	static void join(EventData eventData) {
		eventData.parent.maxDepth = std::max(eventData.children[0].maxDepth, eventData.children[1].maxDepth) + 1;
		eventData.parent.nodesCnt = eventData.children[0].nodesCnt + eventData.children[1].nodesCnt;
		totalJoins++;
	}

	static void split(EventData eventData) {
		totalSplits++;
	}

	void log() {
		if (nodesCnt >= depthPerSize.size()) {
			depthPerSize.resize(nodesCnt + 1);
		}
		AvgMax &item = depthPerSize[nodesCnt];
		if (item.max < maxDepth) {
			item.max = maxDepth;
		}
		item.cnt++;
		item.sum += maxDepth;
		totalLogs++;
	}

	static void print() {
		std::cout << std::endl << "Stats: [tree size: avg depth, max depth, number of records]" << std::endl;
		int i = 0;
		int toI = 1;
		while (i < depthPerSize.size()) {
			int fromI = i;
			toI <<= 1;
			if (toI > depthPerSize.size()) toI = depthPerSize.size();
			int avgMin = std::numeric_limits<int>::max();
			int avgMax = 0;
			int max = 0;
			int cnt = 0;
			for (; i < toI; i++) {
				AvgMax &item = depthPerSize[i];
				cnt += item.cnt;
				int avg = item.cnt ? (item.sum / item.cnt) : 0;
				if (avgMin > avg) avgMin = avg;
				if (avgMax < avg) avgMax = avg;
				if (max < item.max) max = item.max;
			}
			std::cout << fromI << "-" << (toI-1) << ": avg " << avgMin << "-" << avgMax << ", max " << max << ", cnt " << cnt << std::endl;
		}
		std::cout << "avg splits/op: " << (totalSplits / (float) totalLogs) << std::endl;
		std::cout << "avg joins/op:  " << (totalJoins  / (float) totalLogs) << std::endl;
	}
};
std::vector<StatsUserData::AvgMax> StatsUserData::depthPerSize = std::vector<AvgMax>();
int StatsUserData::totalLogs = 0;
int StatsUserData::totalSplits = 0;
int StatsUserData::totalJoins  = 0;


int main() {
	using Vertex = TopTreeVertex;
	ShadowForest shadowForest;
	BiasedTreeTopTree<PathLengthUserData, StatsUserData> forestBTTT;
	TopTree<PathLengthUserData, StatsUserData> &forest = forestBTTT;

	Vertex u = forest.newVertex();
	assert(shadowForest.newVertex() == u);
	int step = 0;

	// create one component of given number of edges
	for (int i = 0; i < 1100; i++) {
		if (rndDist(rndGen) < 0.2) {
			u = rndDist(rndGen) * (i + 1);
		}
		Vertex v = forest.newVertex();
		assert(shadowForest.newVertex() == v);
		std::cout << ++step << ": link(" << u  << ", " << v << ")" << std::endl;
		int edgeLength = rndDist(rndGen) * 1000;
		forest.link(u, v, edgeLength, StatsUserData());
		assert(shadowForest.link(u, v, edgeLength * 2));
		forest.getRootData<1>().log();
		u = v;
	}

	forest.getRootData().multiplyAll(2);

	for (int i = 0; i < 100; i++) {
		while (shadowForest.verticesCnt() - shadowForest.edgesCnt() < 31) {
			// cut random edge
			size_t edgeIndex = rndDist(rndGen) * shadowForest.edgesCnt();
			auto [u, v, length] = shadowForest.getEdge(edgeIndex);
			std::cout << ++step << ": cut(" << u  << ", " << v << ")" << std::endl;
			auto [data, stats] = forest.cut(u, v);
			assert(data.length == length);
			assert(shadowForest.cut(u, v));
			forest.getRootData<1>().log();
		}

		for (int j = 0; j < 10; j++) {
			// choose random pair of vertices; either link or expose them
			Vertex u, v;
			u = rndDist(rndGen) * shadowForest.verticesCnt();
			do {
				v = rndDist(rndGen) * shadowForest.verticesCnt();
			} while (u == v);
			int length = shadowForest.pathLength(u, v);
			if (length < 0) {
				std::cout << ++step << ": link(" << u  << ", " << v << ")" << std::endl;
				int edgeLength = rndDist(rndGen) * 1000;
				forest.link(u, v, edgeLength, StatsUserData());
				forest.getRootData<1>().log();
				assert(shadowForest.link(u, v, edgeLength));
			} else {
				std::cout << ++step << ": expose(" << u << ", " << v << ")" << std::endl;
				forest.expose(u, v);
				auto data = forest.getRootData();
				assert(data.length == length);
				forest.getRootData<1>().log();
			}
		}
	}

	StatsUserData::print();

	return 0;
}
