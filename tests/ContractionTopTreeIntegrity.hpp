// TopTreeLibrary  Copyright (C) 2020  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* Integrity tests for ContractionTopTree class,
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

#ifndef CONTRACTION_TOP_TREE_INTEGRITY_HPP
#define CONTRACTION_TOP_TREE_INTEGRITY_HPP

namespace TopTreeInternals {
	template <class... TUserData>
	class ContractionTopTree;

#define FAIL_HERE { printf("%sIntegrity test failed at " __FILE__ ":%d\n", IntegrityFailed ? "" : "\n\n", __LINE__); IntegrityFailed = true; return false; }
#define TEST(cond) {if (!(cond)) FAIL_HERE; }
#define TEST2(cond) { int I = 0; TEST(cond); I = 1; TEST(cond); }
#define TESTD(cond) { typename Arc::Dir D = Arc::Dir::SUCC; TEST(cond); D = Arc::Dir::PRED; TEST(cond); }
#define TEST2D(cond) { typename Arc::Dir D = Arc::Dir::SUCC; TEST2(cond); D = Arc::Dir::PRED; TEST2(cond); }
#define REQUIRE_LEVEL(level) {if (TopTreeIntegrityLevel < level) return true; }

	template <class... TUserData>
	class ContractionTopTreeIntegrity {
		using TT = TopTree<TUserData...>;
		using CTT = ContractionTopTree<TUserData...>;
		using BaseIntegrity = typename TT::Integrity;
		using Node = typename TT::Node;
		using ENode = typename CTT::ENode;
		using VirtNode = typename CTT::VirtNode;
		using Arc = typename CTT::Arc;

		public:

			// tests consistency of Euler tours in tree and maximality of moves at each level
			static bool treeConsistency(ENode *root) {
				REQUIRE_LEVEL(2);
				TEST2(root->virtNode.arcs[I].isLeaf());
				TEST(arcsConsistency(&root->virtNode));

				Arc *initArc = root->virtNode.childrenArcs[0];
				if (initArc) {
					TEST(childrenConsistency(&root->virtNode));
					TEST(childrenNeighbourhoodConsistency(&root->virtNode));
				}

				while (initArc) {
					bool baseLevel = !initArc->virtNode->childrenArcs[0];
					Arc *arc = initArc;
					do {
						TEST(arcsConsistency(arc->virtNode));
						// maximality of moves
						if (arc->virtNode->isDummy()) {
							for (Arc *arc2 : {arc, arc->twin()}) {
								TEST(  // no missing move violating maximality
									arc2->isLeaf() || !arc2->succ()->virtNode->isDummy() || (
										!arc2->canCompress() &&
										!arc2->canRake()));
							}
						}
						if (!baseLevel) {
							TEST(childrenConsistency(arc->virtNode));
							TEST(childrenNeighbourhoodConsistency(arc->virtNode));
						} else {
							TEST(!arc->virtNode->childrenArcs[0]);
							TEST(arc->virtNode->node->clusterType == BASE);
						}
						arc = arc->succ();
					} while (arc != initArc);
					
					initArc = initArc->virtNode->childrenArcs[0];
				}
				return true;
			}

			template <typename Arc::Dir D>
			static bool arcConsistency(Arc *arc) {
				REQUIRE_LEVEL(1);
				Arc *arc2 = arc->ptrs[D];
				TEST(arc2->ptrs[!D] == arc);
				return true;
			}

			static bool arcsConsistency(VirtNode *node) {
				REQUIRE_LEVEL(1);
				TEST2(node->arcs[I].virtNode == node);
				TEST2(arcConsistency<Arc::Dir::SUCC>(&node->arcs[I]));
				TEST2(arcConsistency<Arc::Dir::PRED>(&node->arcs[I]));
				return true;
			}

			static bool childrenConsistency(VirtNode *node) {
				REQUIRE_LEVEL(1);
				TEST(node->childrenArcs[0]);
				TEST(node->childrenArcs[0]->virtNode->virtParent == node);
				TEST(arcsConsistency(node->childrenArcs[0]->virtNode));
				if (node->isDummy()) {
					TEST(node != CTT::virt(node->node));
				}  else {
					TEST(node == CTT::virt(node->node));
					TEST(node->childrenArcs[1]);
					TEST(node->childrenArcs[1]->virtNode->virtParent == node);
					TEST(arcsConsistency(node->childrenArcs[1]->virtNode));
					TEST(node->isMoveValid());
				}
				return true;
			}

			static bool childrenNeighbourhoodConsistency(VirtNode *node) {
				REQUIRE_LEVEL(1);
				if (node->node->clusterType == BASE) return true;
				TEST2D(node->arcs[I].ptrs[D]->virtNode == node->lowerLevelArcs(I,D)->virtNode->virtParent);
				return true;
			}
			
			static bool nodeConsistency(VirtNode *node) {
				REQUIRE_LEVEL(1);
				TEST(arcsConsistency(node));
				TEST(childrenConsistency(node));
				return true;
			}

	};

}

#endif
#endif
