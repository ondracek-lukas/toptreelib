#ifndef BIASED_TREE_TOP_TREE
#define BIASED_TREE_TOP_TREE

/* TODO:
 *   integrity tests;
 *   splice;
 *   reroot, symmetrize;
 *   link, cut;
 *   biased-tree-balancing join/split;
 */

#include "TopTree.hpp"
#include <iostream> // XXX tmp
#include <algorithm>

#define assert_STR(expr) #expr
#define assert(cond) { \
		if (!(cond)) { \
			printf("Assertion [" assert_STR(cond) "] failed at " __FILE__ ":%d\n\n", __LINE__); \
			fflush(stdout); \
			abort(); \
		}}

namespace TopTreeInternals {

	template <class... TUserData>
	class BiasedTreeTopTree : public TopTree<TUserData...> {
		protected:
			using Integrity = typename TopTree<TUserData...>::Integrity; // TODO EIntegrity for BTTT
		private:
			struct ENode;

			template <ClusterType TTreeType = COMPRESS>
			static void debugPrintTree(ENode *node, bool reversed = false, bool isRoot = true, bool newline=true, ENode *btRoot = nullptr) {
				if (!btRoot) btRoot = node;
				if (isRoot) {
					std::cout << (TTreeType == COMPRESS ? "C[ " : "R[ ");
				}
				if (node) {
					if (node->template isLeaf<TTreeType>(btRoot)) {
						if (node->clusterType == BASE) {
							std::cout << "(" << node->boundary[reversed] << " " << node->boundary[!reversed] << ") ";
						} else {
							debugPrintTree<TTreeType == COMPRESS ? RAKE : COMPRESS>(node, false, true, false);
						}
					} else {
						for (int i : {reversed, !reversed}) {
							bool rev = false;
							if constexpr(TTreeType == COMPRESS) {
								rev = reversed ^ node->isChildReversed(i);
							}
							debugPrintTree<TTreeType>(ext(node->children[i]), rev, false, false, btRoot);
						}
					}
				} else {
					std::cout << "null ";
				}
				if (isRoot) std::cout << "] ";
				if (newline) std::cout << std::endl;
			}
			static void debugPrintRawTree(ENode *node, int indent = 0) {
				for (int i = 0; i < indent; i++) std::cout << "  ";
				char t;
				switch (node->clusterType) {
					case BASE: t = 'B'; break;
					case COMPRESS: t = 'C'; break;
					case RAKE: t = 'R'; break;
				}
				std::cout << t << "(" << node->boundary[0] << " " << node->boundary[1] << " :" << node->rakeMaxWeight << " " << node->template isLeaf<COMPRESS>() << "," << node->template isLeaf<RAKE>() << ")" << std::endl;
				if (node->clusterType != BASE) {
					debugPrintRawTree(ext(node->children[0]), indent+1);
					debugPrintRawTree(ext(node->children[1]), indent+1);
				}
			}

		public:
			BiasedTreeTopTree() {
				ENode *nodes[11];
				Vertex v = this->newVertex();
				for (auto &node : nodes) {
					node = newNode();
					Vertex u = v; v = this->newVertex();
					node->setBoundary(u, v);
				}

				ENode *tree1 = nullptr, *tree2 = nullptr;
				for (int i = 0; i < 5; i++) {
					tree1 = join<COMPRESS>(tree1, nodes[i], nullptr);       // C[ 0--5 ]
					tree2 = join<COMPRESS>(nullptr, nodes[9-i], tree2);     // C[ 5--10 ]
				}
				tree2 = join<RAKE>(nullptr, tree2, nullptr);              // R[ C[ 5--10 ] ]
				auto tree12 = join<RAKE>(tree2, tree1, nullptr);          // R[ C[ 5--10 ] C[ 0--5 ] ]
				auto tree3 = join<COMPRESS>(nullptr, nodes[10], nullptr); // C[ 10--11 ]
				auto tree123 = join<COMPRESS>(nullptr, tree12, tree3);    // C[ R[ C[ 5--10 ] C[ 0--5 ] ] 10--11 ]

				debugPrintTree(tree123);
				std::cout << std::endl;

				auto [A, B, C] = split<COMPRESS>(tree123, tree12);
				debugPrintTree<COMPRESS>(A);  // C[ null ]
				debugPrintTree<RAKE>(B);      // ~ tree12
				debugPrintTree<COMPRESS>(C);  // ~ tree3
				std::cout << std::endl;

				{
					ENode *tmpNode = newNode();
					tmpNode->weight = 1000;
					auto [D, E, F] = split<RAKE>(B, tmpNode);
					debugPrintTree<RAKE>(D);      // ~ tree12
					debugPrintTree<COMPRESS>(E);  // C[ null ]
					debugPrintTree<RAKE>(F);      // R[ null ]
					deleteNode(tmpNode);
					std::cout << std::endl;
				}

				{
					ENode *tmpNode = newNode();
					tmpNode->weight = 0;
					auto [D, E, F] = split<RAKE>(B, tmpNode);
					debugPrintTree<RAKE>(D);      // R[ null ]
					debugPrintTree<COMPRESS>(E);  // C[ null ]
					debugPrintTree<RAKE>(F);      // ~ tree12
					deleteNode(tmpNode);
					std::cout << std::endl;
				}

				auto [D, E, F] = split<RAKE>(B, tree1);
				debugPrintTree<RAKE>(D);      // ~ tree2
				debugPrintTree<COMPRESS>(E);  // ~ tree1
				debugPrintTree<RAKE>(F);      // R[ null ]
				std::cout << std::endl;

				auto [G, H, I] = split<COMPRESS>(E, nodes[2]);
				debugPrintTree<COMPRESS>(G);  // C[ 0--2 ]
				debugPrintTree<COMPRESS>(H);  // C[ 2--3 ]
				debugPrintTree<COMPRESS>(I);  // C[ 3--5 ]
				std::cout << std::endl;

				
			}

			virtual void link(Vertex u, Vertex v, TUserData... userData) {
				assert(false);
			}

			virtual std::tuple<TUserData...> cut(Vertex u, Vertex v) {
				assert(false);
			}

		private:
			using Node = typename TopTree<TUserData...>::Node;
			struct ENode : public Node { // Extended Node
				int weight, rank, rankDeferred;
				int rakeMaxWeight; // 0 iff not in rake tree
				int minWeightDiff[2]; // compress tree only
				int minWeightDiffDeferred[2];

				bool isChildReversed(int i) { // compress tree only
					return this->boundary[i] != this->children[i]->boundary[i];
				}
				bool isRakeTreeNode() {
					return rakeMaxWeight;
				}
				template <ClusterType TTreeType>
				bool isLeaf(Node *root = nullptr) {
					if constexpr(TTreeType == COMPRESS) {
						return (this->clusterType == BASE) || ((this != root) && isRakeTreeNode());
					} else { // RAKE
						return (this->clusterType != RAKE) || !ext(this->children[0])->isRakeTreeNode();
					}
				}


				void setBoundary(Vertex v1, Vertex v2); // for BASE node

				template <ClusterType TTreeType>
				void setBoundary(ClusterType type); // for COMPRESS or RAKE node

				template <ClusterType TTreeType, bool becomingLeaf>
				void setLeafData();

				template <ClusterType TTreeType>
				ENode *attachChildren(ClusterType type, ENode *child1, ENode *child2);

				template <ClusterType TTreeType>
				ENode *detach();
			};


			static ENode *ext(Node *node) { return static_cast<ENode *>(node); };

			InnerList<Node, &Node::tmpListNode> freeENodes;
			ENode *newNode() { return this->template newNodeT<ENode>(freeENodes); }
			void deleteNode(Node *node) { this->template deleteNodeT<ENode>(node, freeENodes); }


			template <ClusterType TTreeType>
			ENode *join(ENode *leftTree, ENode *middleNode, ENode *rightTree);

			template <ClusterType TTreeType>
			std::tuple<ENode *, ENode *, ENode *> split(ENode *tree, ENode *middleNode);
	};

	// --- data manipulation, split, join ---

	template <class... TUserData>
	inline void BiasedTreeTopTree<TUserData...>::ENode::
	setBoundary(Vertex v1, Vertex v2) { // for BASE node, for COMPRESS tree
		static_cast<Node *>(this)->setBoundary(v1, v2);
		// set compress tree base node data
		weight = 1;
		rakeMaxWeight = 0;
	}

	template <class... TUserData>
	template <ClusterType TTreeType>
	inline void BiasedTreeTopTree<TUserData...>::ENode::
	setBoundary(ClusterType type) { // for COMPRESS or RAKE node
		static_cast<Node *>(this)->setBoundary(type);
		// aggregate data within one rake/compress tree
		weight = ext(this->children[0])->weight + ext(this->children[1])->weight;
		if constexpr(TTreeType == RAKE) {
			rakeMaxWeight = std::max(ext(this->children[0])->rakeMaxWeight, ext(this->children[1])->rakeMaxWeight);
		} else {
			bool rev[2] = {isChildReversed(0), isChildReversed(1)};
			for (int i : {0, 1}) {
				ext(this->children[!i])->minWeightDiff[i ^ rev[!i]] += ext(this->children[i])->weight;
				minWeightDiff[i] = std::min(
					ext(this->children[0])->minWeightDiff[i ^ rev[0]],
					ext(this->children[1])->minWeightDiff[i ^ rev[1]]);
				ext(this->children[0])->minWeightDiff[i ^ rev[0]] -= minWeightDiff[i];
				ext(this->children[1])->minWeightDiff[i ^ rev[1]] -= minWeightDiff[i];
			}
			rakeMaxWeight = 0;
		}
	}

	template <class... TUserData>
	template <ClusterType TTreeType>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::ENode::
	detach() {
		// propagate data within one rake/compress tree
		if constexpr(TTreeType == COMPRESS) {
			if (this->parent && this->parent->children[0] && this->parent->children[1]) {
				bool rev[2] = {ext(this->parent)->isChildReversed(0), ext(this->parent)->isChildReversed(1)};
				for (int i : {0, 1}) {
					ext(this->parent->children[0])->minWeightDiff[i ^ rev[0]] += ext(this->parent)->minWeightDiff[i];
					ext(this->parent->children[1])->minWeightDiff[i ^ rev[1]] += ext(this->parent)->minWeightDiff[i];
					ext(this->parent->children[!i])->minWeightDiff[i ^ rev[!i]] -= ext(this->parent->children[i])->weight;
				}
			}
		}
		return ext(static_cast<Node *>(this)->detach());
	}


	template <class... TUserData>
	template <ClusterType TTreeType, bool becomingLeaf>
	inline void BiasedTreeTopTree<TUserData...>::ENode::
	setLeafData() {
		// initialize rake/compress tree leaf data or restore the original
		if constexpr(TTreeType == RAKE) {
			if constexpr(becomingLeaf) {
				minWeightDiffDeferred[0] = minWeightDiff[0];
				minWeightDiffDeferred[1] = minWeightDiff[1];
				rakeMaxWeight = weight;
			} else {
				minWeightDiff[0] = minWeightDiffDeferred[0];
				minWeightDiff[1] = minWeightDiffDeferred[1];
				rakeMaxWeight = 0;
			}
		} else {
			if constexpr(becomingLeaf) {
				minWeightDiff[0] = minWeightDiff[1] = -rakeMaxWeight;
			}
		}
	}

	template <class... TUserData>
	template <ClusterType TTreeType>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::ENode::
	attachChildren(ClusterType type, ENode *child1, ENode *child2) {
		this->attachChild(0, child1);
		this->attachChild(1, child2);
		setBoundary<TTreeType>(type);
		return this;
	}

	template <class... TUserData>
	template <ClusterType TTreeType>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::
	join(ENode *leftTree, ENode *middleNode, ENode *rightTree) {
		// TODO reimplement as biased trees (currently unbalanced)

		if (middleNode) {
			middleNode->template setLeafData<TTreeType, true>();
			ENode *node = newNode();
			ClusterType nodeType = (TTreeType == RAKE) || middleNode->isRakeTreeNode() ? RAKE : COMPRESS;
			if (leftTree) {
				node->template attachChildren<TTreeType>(nodeType, leftTree, middleNode);
				leftTree = node;
			} else if (rightTree) {
				if constexpr(TTreeType == RAKE) {
					node->template attachChildren<TTreeType>(nodeType, middleNode, rightTree); // preserve order
				} else {
					node->template attachChildren<TTreeType>(nodeType, rightTree, middleNode);
				}
				rightTree = node;
			} else {
				return middleNode;
			}
		}

		if (leftTree == nullptr) {
			assert(Integrity::treeConsistency(rightTree));
			return rightTree;
		} else if (rightTree == nullptr) {
			assert(Integrity::treeConsistency(leftTree));
			return leftTree;
		} else {
			ENode *root = newNode();
			root->template attachChildren<TTreeType>(TTreeType, leftTree, rightTree);
			assert(Integrity::treeConsistency(root));
			return root;
		}
	}

	template <class... TUserData>
	template <ClusterType TTreeType>
	std::tuple<typename BiasedTreeTopTree<TUserData...>::ENode *, typename BiasedTreeTopTree<TUserData...>::ENode *, typename BiasedTreeTopTree<TUserData...>::ENode *> BiasedTreeTopTree<TUserData...>::
	split(ENode *tree, ENode *middleNode) {
		// TODO reimplement as biased trees (currently unbalanced)

		if constexpr(TTreeType == RAKE) {
			if (!middleNode->parent) {
				// middle node is not in rake tree, find corresponding inner vertex by weights of leaves
				int weight = middleNode->weight;
				middleNode = nullptr;
				if (tree->rakeMaxWeight < weight) {
					return {tree, nullptr, nullptr};
				}
				ENode *node = tree;
				while (!node->template isLeaf<TTreeType>()) {
					if (ext(node->children[0])->rakeMaxWeight < weight) {
						middleNode = node;
						node = ext(node->children[1]);
					} else {
						node = ext(node->children[0]);
					}
				}
				if (!middleNode) {
					return {nullptr, nullptr, tree};
				}
			}
		}

		ENode *middleRet = nullptr, *left = nullptr, *right = nullptr;
		if (middleNode->template isLeaf<TTreeType>()) {
			middleRet = middleNode;
			middleNode = ext(middleRet->parent);
			TopTree<TUserData...>::releaseNode(middleNode);
			middleRet->template detach<TTreeType>();
			middleRet->template setLeafData<TTreeType, false>();
		} else {
			TopTree<TUserData...>::releaseNode(middleNode);
			left = ext(middleNode->children[0]);
			left->template detach<TTreeType>();
		}

		while (middleNode) {
			int childI = middleNode->children[0] == nullptr;
			ENode *child = ext(middleNode->children[childI]);
			child->template detach<TTreeType>();
			if (childI) {
				right  = join<TTreeType>(nullptr, child, right);
			} else {
				left   = join<TTreeType>(left,    child, nullptr);
			}
			ENode *parent = ext(middleNode->parent);
			middleNode->template detach<TTreeType>();
			deleteNode(middleNode);
			middleNode = parent;
		}

		assert(Integrity::treeConsistency(left));
		assert(Integrity::treeConsistency(right));
		
		return {left, middleRet, right};
	}

}

using TopTreeInternals::BiasedTreeTopTree;

#undef assert
#endif
