// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* This will became a unit test in future.
 */

#define TOP_TREE_INTEGRITY_LEVEL 2
#include "../TopTree.hpp"

#include <cassert>

template <class... TUserData>
class TestingTopTree : public TopTree<TUserData...> {
	using Vertex = TopTreeVertex;
	using ClusterType = TopTreeClusterType;
	using Node = typename TopTree<TUserData...>::Node;

	public:
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

		virtual void link(Vertex u, Vertex v, TUserData... userData) {
			assert(0);
		}
		virtual std::tuple<TUserData...> cut(Vertex u, Vertex v) {
			assert(0);
		}
};

int main() {
	TestingTopTree<> tree;
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			printf("%d %d %d\n", i, j, tree.sameComponent(i, j));
	return 0;
}
