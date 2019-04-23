// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* Integrity tests for TopTree class.
 * Set TOP_TREE_INTEGRITY_LEVEL macro to enable:
 *   0  disabled,
 *   1  tests preserving asymptotic complexity only,
 *   2  all tests enabled.
 */

#if TOP_TREE_INTEGRITY_LEVEL > 0
#ifndef TOP_TREE_INTEGRITY
#define TOP_TREE_INTEGRITY

#include <cstdio>

template <class... TUserData>
class TopTree;

#define assert_STR(expr) #expr
#define assert(cond) { \
	if (!(cond)) { \
		printf("Assertion [" assert_STR(cond) "] failed at " __FILE__ ":%d\n\n", __LINE__); \
		fflush(stdout); \
		abort(); \
	}}

#define FAIL_HERE { printf("\n\nIntegrity test failed at " __FILE__ ":%d\n", __LINE__); return false; }
#define TEST(cond) {if (!(cond)) FAIL_HERE; }
#define TEST2(cond) { int I = 0; TEST(cond); I = 1; TEST(cond); }
#define REQUIRE_LEVEL(level) {if (TOP_TREE_INTEGRITY_LEVEL < level) return true; }

template <class... TUserData>
class TopTreeIntegrity {
	using Vertex = typename TopTree<TUserData...>::Vertex;
	using ClusterType = typename TopTree<TUserData...>::ClusterType;
	using Node = typename TopTree<TUserData...>::Node;

	public:

		// tests constraints between node's and its children's boundary
		static bool nodeChildrenBoundary(Node *node) {
			if (node->clusterType == ClusterType::BASE) return true;
			TEST2(node->children[I]);

			int boundaryInChildren[2] = {
				node->children[0]->boundary[1] == node->boundary[0],
				node->children[1]->boundary[1] == node->boundary[1]};

			TEST2(node->children[I]->boundary[boundaryInChildren[I]] == node->boundary[I]);

			if (node->clusterType == ClusterType::COMPRESS) {
				TEST(
					node->children[0]->boundary[!boundaryInChildren[0]] ==
					node->children[1]->boundary[!boundaryInChildren[1]]);
				TEST(boundaryInChildren[0] == 0); // same order as in children; is needed?
			} else {
				TEST(node->clusterType == ClusterType::RAKE);
				TEST(node->children[0]->boundary[!boundaryInChildren[0]] == node->boundary[1]);
			}

			return true;
		}

		// tests parent/children pointers consistency and nodeChildrenBoundary in each node
		static bool treeConsistency(Node *root) {
			REQUIRE_LEVEL(2);
			if (root == nullptr) return true;

			TEST(root->parent == nullptr);

			Node *node = root;
			do {
				if (node->clusterType == ClusterType::BASE) {

					TEST2(node->children[I] == nullptr);

					while ((node != root) && (node->parent->children[1] == node)) node = node->parent;
					if (node != root) node = node->parent->children[1];

				} else {

					TEST2((node->children[I]) && (node->children[I]->parent == node));
					TEST(nodeChildrenBoundary(node));

					node = node->children[0];

				}
			} while (node != root);

			return true;
		}

};

#undef TEST2
#undef TEST
#undef FAIL_HERE
#endif
#endif