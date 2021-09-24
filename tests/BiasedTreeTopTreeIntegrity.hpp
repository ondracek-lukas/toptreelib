// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* Integrity tests for BiasedTreeTopTree class,
 * it uses TopTreeIntegrity for the base class.
 *
 * Set TOP_TREE_INTEGRITY macro to enable tests.
 * Adjust global variable TopTreeIntegrityLevel to:
 *   0  disable most of the tests,
 *   1  enable tests preserving asymptotic complexity only,
 *   2  enable all tests (default).
 */

#ifdef TOP_TREE_INTEGRITY
#include "TopTreeIntegrity.hpp"

#ifndef BIASED_TREE_TOP_TREE_INTEGRITY_HPP
#define BIASED_TREE_TOP_TREE_INTEGRITY_HPP

namespace TopTreeInternals {
	template <class... TUserData>
	class BiasedTreeTopTree;

#define FAIL_HERE { printf("\n\nIntegrity test failed at " __FILE__ ":%d\n", __LINE__); return false; }
#define TEST(cond) {if (!(cond)) FAIL_HERE; }
#define TEST2(cond) { int I = 0; TEST(cond); I = 1; TEST(cond); }
#define REQUIRE_LEVEL(level) {if (TopTreeIntegrityLevel < level) return true; }

	template <class... TUserData>
	class BiasedTreeTopTreeIntegrity {
		using BaseIntegrity = typename TopTree<TUserData...>::Integrity;
		using BTTT = BiasedTreeTopTree<TUserData...>;
		using Node = typename TopTree<TUserData...>::Node;
		using ENode = typename BiasedTreeTopTree<TUserData...>::ENode;

		public:

			// tests extended node data in addition to base class consistency test
			static bool treeConsistency(ENode *root) {
				// std::cout << "-- treeConsistency --" << std::endl;
				TEST(BaseIntegrity::treeConsistency(root));
				REQUIRE_LEVEL(2);
				if (root == nullptr) return true;

				struct Accum {
					int minWeightDiff[2];
					int bottomWeight[2];
					bool reversed;

					bool btRootLeaf;
					int btBoundingLeavesRanks[2];
					int btMinReqNeighbouringLeavesRanks[2];
				};
				std::vector<Accum> accum;

				for (Node *baseClassNode : root->preorder()) {
					ENode *node        = BTTT::ext(baseClassNode);
					ENode *parent      = BTTT::ext(node->parent);
					ENode *sibling     = parent ? BTTT::ext(parent->children[ parent->children[0] == node ]) : nullptr;
					ENode *children[2] = {BTTT::ext(node->children[0]), BTTT::ext(node->children[1])};

					bool top_root               = !parent;
					bool top_leaf               = node->clusterType == BASE;
					bool rake_node              = node->rakeMaxWeight;
					bool rake_leaf              = rake_node && (top_leaf || !children[0]->rakeMaxWeight || !children[1]->rakeMaxWeight);
					bool rake_root              = rake_node && (top_root || !sibling->rakeMaxWeight);
					bool compress_root          = rake_leaf || (top_root && !rake_node);
					bool compress_leaf_root     = compress_root && top_leaf;
					bool compress_leaf_non_root = (rake_root && !top_root) || (top_leaf && !compress_root);
					bool compress_node          = compress_root || compress_leaf_root || compress_leaf_non_root || !rake_node;

					// compress_leaf ambiguity:
					//   node of two compress trees: compress_root && rake_root && !top_root
					//   trivial compress tree:      compress_root && top_leaf ~~ compress_leaf_root

					{ // resize accum if needed
						int index = node->index;
						if (!top_root && (index < parent->index)) index = parent->index;
						if (!top_leaf) {
							if (index < children[0]->index) index = children[0]->index;
							if (index < children[1]->index) index = children[1]->index;
						}
						if (accum.size() <= index) accum.resize(index + 1);
					}

					accum[node->index].btRootLeaf = !top_leaf && (rake_leaf || compress_leaf_non_root);

					// weight
					if (top_leaf) {
						TEST(node->weight == 1);
					} else {
						TEST(node->weight == children[0]->weight + children[1]->weight);
					}

					// rakeMaxWeight
					if (rake_node) {
						if (rake_leaf) {
							TEST(node->rakeMaxWeight == node->weight);
						} else {
							TEST(node->rakeMaxWeight == std::max(children[0]->rakeMaxWeight, children[1]->rakeMaxWeight));
						}
					}

					// accum -- minWeightDiff, minWeightDiffDeferred
					if (compress_node) {
						if (compress_leaf_non_root) {
							// std::cout << "(" << node->boundary[0] << " " << node->boundary[1] << ")" << std::endl;
							// for (int I : {0, 1}) std::cout << accum[node->index].minWeightDiff[I] << " == " << accum[node->index].bottomWeight[I] << " - " << node->rakeMaxWeight << std::endl;
							TEST2(accum[node->index].minWeightDiff[I] == accum[node->index].bottomWeight[I] - node->rakeMaxWeight);
						}
						if (compress_root) {
							accum[node->index].reversed = false;
							for (int i : {0, 1}) {
								accum[node->index].bottomWeight[i] = 0;
								if (rake_leaf) {
									accum[node->index].minWeightDiff[i] = node->minWeightDiffDeferred[i];
								} else {
									accum[node->index].minWeightDiff[i] = node->minWeightDiff[i];
								}
							}
							if (compress_leaf_root) {
								// std::cout << accum[node->index].minWeightDiff[0] << " == 0" << std::endl;
								// std::cout << accum[node->index].minWeightDiff[1] << " == 0" << std::endl;
								TEST2(accum[node->index].minWeightDiff[I] == 0);
							}
						}
						if (!compress_leaf_root && (!compress_leaf_non_root || compress_root)) {
							for (int i : {0, 1}) {
								bool rev = accum[node->index].reversed;
								bool childRev = rev ^ (node->boundary[i] != children[i]->boundary[i]);
								accum[children[i]->index].reversed = childRev;
								accum[children[i]->index].bottomWeight[i ^ childRev ^ rev]   = accum[node->index].bottomWeight[i];
								accum[children[i]->index].bottomWeight[!i ^ childRev ^ rev]  = accum[node->index].bottomWeight[!i] + children[!i]->weight;
								accum[children[i]->index].minWeightDiff[childRev]  = accum[node->index].minWeightDiff[rev]  + children[i]->minWeightDiff[childRev];
								accum[children[i]->index].minWeightDiff[!childRev] = accum[node->index].minWeightDiff[!rev] + children[i]->minWeightDiff[!childRev];
							}
						}
						// if (!top_leaf) {
						// 	std::cout << node->isChildReversed(0) << ", " << node->isChildReversed(1) << std::endl;
						// }
						// std::cout << node->boundary[0] << " " << node->boundary[1] << ": " << node->minWeightDiff[0] << ", " << node->minWeightDiff[1] << "; " << node->minWeightDiffDeferred[0] << ", " << node->minWeightDiffDeferred[1] << "; " << accum[node->index].minWeightDiff[0] << ", " << accum[node->index].minWeightDiff[1] << "; " << accum[node->index].reversed << std::endl;
						// std::cout << accum[node->index].bottomWeight[0] << ", " << accum[node->index].bottomWeight[1] << std::endl;
						// std::cout << std::endl;
					}

					// types of nodes and shared vertices in rake trees
					if (compress_node && (compress_root || !rake_root) && node->clusterType == RAKE) {
						TEST(!children[0]->rakeMaxWeight && children[1]->rakeMaxWeight);
						TEST((children[0]->clusterType != RAKE) || (children[0]->boundary[1] != node->boundary[1]));
					}
					if (rake_node && !rake_leaf) {
						TEST(node->clusterType == RAKE);
						if (!top_root) {
							TEST(node->boundary[1] == parent->boundary[1]);
						}
					}
				}

				// biased tree invariants,  TODO: verify
				for (Node *baseClassNode : root->postorder()) {
					ENode *node        = BTTT::ext(baseClassNode);
					Accum *nAccum      = &accum[node->index];
					ENode *parent      = !accum[node->index].btRootLeaf ? BTTT::ext(node->parent) : nullptr;
					ENode *grandparent = parent && !accum[parent->index].btRootLeaf ? BTTT::ext(parent->parent) : nullptr;
					bool  btLeaf       = (node->clusterType == BASE) || (accum[node->index].btRootLeaf);
					bool  rev          = nAccum->reversed;

					if (node->clusterType != BASE) {
						bool  childRev[2]  = {accum[node->children[0]->index].reversed, accum[node->children[1]->index].reversed};
						// 4th int. -- global bias property
						TEST2(
							accum[node->children[I^rev]->index].btBoundingLeavesRanks[!I ^ childRev[I^rev]] >=
								accum[node->children[!I^rev]->index].btMinReqNeighbouringLeavesRanks[I ^ childRev[!I^rev]]);

						if (!btLeaf) {
							for (int i : {0, 1}) {
								nAccum->btBoundingLeavesRanks[i^rev] = accum[node->children[i^rev]->index].btBoundingLeavesRanks[i ^ childRev[i^rev]];
								nAccum->btMinReqNeighbouringLeavesRanks[i^rev] = accum[node->children[i^rev]->index].btMinReqNeighbouringLeavesRanks[i ^ childRev[i^rev]];
							}
						}
					}

					if (btLeaf) {
						for (int i : {0, 1}) {
							nAccum->btBoundingLeavesRanks[i] = node->rank;
							nAccum->btMinReqNeighbouringLeavesRanks[i] = 0;
						}

						// 1st inv.
						TEST((1 << node->rank <= node->weight) && (2 << node->rank > node->weight));
					}

					// 2nd inv.
					TEST(!parent || ((accum[parent->index].btRootLeaf ? parent->rankDeferred : parent->rank) >= node->rank + btLeaf));

					// 3rd inv.
					TEST(!grandparent || ((accum[grandparent->index].btRootLeaf ? grandparent->rankDeferred : grandparent->rank) >= node->rank + 1));

					if (parent && ((accum[parent->index].btRootLeaf ? parent->rankDeferred : parent->rank) >= node->rank + 2)) { // minor node
						for (int i : {0, 1}) {
							nAccum->btMinReqNeighbouringLeavesRanks[i] = node->rank - 1;
						}
					}
				}
				return true;
			}
	};

}

#endif
#endif
