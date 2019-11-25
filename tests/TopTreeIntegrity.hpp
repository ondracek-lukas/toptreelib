// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* Integrity tests for TopTree class.
 * Set TOP_TREE_INTEGRITY macro to enable tests.
 * Adjust global variable TopTreeIntegrityLevel to:
 *   0  disable most of the tests,
 *   1  enable tests preserving asymptotic complexity only,
 *   2  enable all tests (default).
 */

#ifdef TOP_TREE_INTEGRITY

#define assert_STR(expr) #expr
#define assert(cond) { \
		if (!(cond)) { \
			printf("Assertion [" assert_STR(cond) "] failed at " __FILE__ ":%d\n\n", __LINE__); \
			fflush(stdout); \
			abort(); \
		}}

#ifndef TOP_TREE_INTEGRITY_HPP
#define TOP_TREE_INTEGRITY_HPP

#include <cstdio>
#include <vector>
size_t TopTreeIntegrityLevel = 2;


namespace TopTreeInternals {
	template <class... TUserData>
	class TopTree;

#define FAIL_HERE { printf("\n\nIntegrity test failed at " __FILE__ ":%d\n", __LINE__); return false; }
#define TEST(cond) {if (!(cond)) FAIL_HERE; }
#define TEST2(cond) { int I = 0; TEST(cond); I = 1; TEST(cond); }
#define REQUIRE_LEVEL(level) {if (TopTreeIntegrityLevel < level) return true; }

	template <class... TUserData>
	class TopTreeIntegrity {
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
				} else {
					TEST(node->clusterType == ClusterType::RAKE);
					TEST(node->children[0]->boundary[!boundaryInChildren[0]] == node->boundary[1]);
				}

				return true;
			}

			// tests parent/children pointers consistency and nodeChildrenBoundary in each node;
			// also tests correctness of tree traversal
			static bool treeConsistency(Node *root) {
				REQUIRE_LEVEL(2);
				if (root == nullptr) return true;

				//std::vector<int> nodeVisit; // 0 not seen, 1 seen as child, 2 visited
				class NodeVisit : public std::vector<int> {
					public:
						bool testSet(Node *node, int oldVal, int newVal) {
							int i = node->index;
							if (this->size() <= i) this->resize(i + 1);
							if ((*this)[i] != oldVal) return false;
							(*this)[i] = newVal;
							return true;
						}
				} nodeVisit;

				TEST(root->parent == nullptr);
				TEST(nodeVisit.testSet(root, 0, 1));

				for (Node *node : root->preorder()) {
					TEST(nodeVisit.testSet(node, 1, 2));
					if (node->clusterType == ClusterType::BASE) {
						TEST2(node->children[I] == nullptr);
					} else {
						TEST2((node->children[I]) && (node->children[I]->parent == node));
						TEST(nodeChildrenBoundary(node));
						TEST2(nodeVisit.testSet(node->children[I], 0, 1));
					}
				}

				for (int visit : nodeVisit) {
					TEST(visit != 1);
				}

				return true;
			}

	};

#undef TEST2
#undef TEST
#undef FAIL_HERE
}

#endif
#endif
