// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#ifndef BIASED_TREE_TOP_TREE_HPP
#define BIASED_TREE_TOP_TREE_HPP

/* TODO:
 *   cut;
 *   biased-tree-balancing join/split;
 */

#include "TopTree.hpp"
#include <iostream> // XXX tmp
#include <algorithm>

#ifdef TOP_TREE_INTEGRITY
#include "tests/BiasedTreeTopTreeIntegrity.hpp"
#else
namespace TopTreeInternals {
	template <class... TUserData> class BiasedTreeTopTreeIntegrity {};
}
#define assert(cond)
#endif

namespace TopTreeInternals {

	template <class... TUserData>
	class BiasedTreeTopTree : public TopTree<TUserData...> {
		protected:
			using Integrity = typename TopTree<TUserData...>::Integrity;
			using EIntegrity = BiasedTreeTopTreeIntegrity<TUserData...>; friend EIntegrity;
		private:
			struct ENode;

			template <ClusterType TTreeType = COMPRESS>
			static void debugPrintTree(ENode *node, bool reversed = false, bool isRoot = true, bool newline=true, ENode *btRoot = nullptr) {
				if (!btRoot) btRoot = node;
				if (isRoot) {
					std::cout << (TTreeType == COMPRESS ? "C" : "R") << "<";
					if (node) {
						std::cout << node->boundary[0] << " " << node->boundary[1];
					}
					std::cout << ">[ ";
				}
				if (node) {
					if (node->template isLeaf<TTreeType>(btRoot)) {
						if ((node->clusterType == BASE) && ((TTreeType == RAKE) || (!node->isRakeTreeNode()))) {
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
			BiasedTreeTopTree(bool testingTree = false) { // XXX temporary
				if (!testingTree) return;
				ENode *nodes[11];
				Vertex v = this->newVertex();
				for (auto &node : nodes) {
					node = newNode();
					Vertex u = v; v = this->newVertex();
					if (v == 11) u = 5;
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

			virtual void link(Vertex u, Vertex v, TUserData... userData);
			virtual std::tuple<TUserData...> cut(Vertex u, Vertex v);

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

				void reverse() {
					std::swap(this->boundary[0], this->boundary[1]);
					std::swap(this->children[0], this->children[1]);
					std::swap(minWeightDiffDeferred[0], minWeightDiffDeferred[1]);
					std::swap(minWeightDiff[0], minWeightDiff[1]);
				}

				bool isRakeTreeNode() {
					return rakeMaxWeight;
				}

				template <ClusterType TTreeType>
				bool isLeaf(Node *root = nullptr) { // when going top-down
					if constexpr(TTreeType == COMPRESS) {
						return (this->clusterType == BASE) || ((this != root) && isRakeTreeNode());
					} else { // RAKE
						return (this->clusterType != RAKE) || !ext(this->children[0])->isRakeTreeNode();
					}
				}

				template <ClusterType TTreeType>
				bool isRoot() { // when going bottom-up
					if constexpr(TTreeType == COMPRESS) {
						return !this->parent || (isRakeTreeNode() && isLeaf<RAKE>());
					} else {
						assert(this->parent);
						return !ext(this->parent->children[0])->isRakeTreeNode();
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

			ENode *reroot(Vertex v);
			ENode *rerootRecurse(Vertex v, InnerList<Node, &Node::tmpListNode> &splitNodes);
			ENode *symmetrize(Vertex v);
			ENode *symmetrizeRecurse(ENode *cRoot);
			std::tuple<ENode *, ENode *, ENode *, ENode *> spliceSplit(ENode *cRoot, ENode *cSplit, ENode *rSplit);
			ENode *spliceJoin(ENode *bottom, ENode *altBottom, ENode *rake, ENode *top);

			template <ClusterType TTreeType>
			ENode *findRoot(ENode *node, bool strictPred = false) {
				while (node && (strictPred || !node->template isRoot<TTreeType>())) {
					strictPred = false;
					node = ext(node->parent);
				}
				return node;
			}

			template <ClusterType TTreeType>
			ENode *findHeavyNode(ENode *root) {
				if constexpr(TTreeType == COMPRESS) {
					if (root->minWeightDiff[0] >= 0) {
						return nullptr;
					}
				}
				ENode *node = root;
				bool rev = false;
				int minWeightDiffAccum = 0;
				while (!node->template isLeaf<TTreeType>(root)) {
					int nextI = 0;
					if constexpr(TTreeType == COMPRESS) {
						minWeightDiffAccum += node->minWeightDiff[rev];
						nextI = !rev;
						bool childRev = rev ^ node->isChildReversed(nextI);
						nextI = rev ^ (minWeightDiffAccum + ext(node->children[nextI])->minWeightDiff[childRev] < 0);
						rev ^= node->isChildReversed(nextI);
					}
					node = ext(node->children[nextI]);
				}
				if constexpr(TTreeType == COMPRESS) {
					assert(node->isRakeTreeNode());
				}
				return node;
			}

			// functions for tree manipulation incl. vertexToNode update;
			// to be modified while implementing biased trees split/join
			template <ClusterType TTreeType>
			ENode *joinNodes(ClusterType clusterType, ENode *left, ENode *right) {
				ENode *node = newNode();
				node->template attachChildren<TTreeType>(clusterType, left, right);
				this->vertexToNodeUpdateInner(node);
				return node;
			}
			template <ClusterType TTreeType>
			std::tuple<ENode *, ENode *> splitNode(ENode *node) {
				this->releaseJustNode(node);
				this->markVertexAsIsolated(node->getInnerVertex());
				ENode *leftChild = ext(node->children[0]);
				ENode *rightChild = ext(node->children[1]);
				leftChild->template detach<TTreeType>();
				rightChild->template detach<TTreeType>();
				deleteNode(node);
				return {leftChild, rightChild};
			}

	};

	// --- link, cut ---

	template <class... TUserData>
	void BiasedTreeTopTree<TUserData...>::
	link(Vertex u, Vertex v, TUserData... userData) {
		this->rollback();
		ENode *uTree = ext(this->treeRoot(this->vertexToNode[u]));
		ENode *vTree = ext(this->treeRoot(this->vertexToNode[v]));
		assert((uTree != vTree) || !uTree); // interface assert, maybe use exception instead

		/*
		std::cout << std::endl;
		std::cout << "uTree: ";
		debugPrintTree(uTree);
		std::cout << "vTree: ";
		debugPrintTree(vTree);
		*/

		if (uTree) {
			uTree = reroot(u);
			this->markNodeAsNonRoot(uTree);
		}
		if (vTree) {
			vTree = reroot(v);
			this->markNodeAsNonRoot(vTree);
		}
		if (!uTree || (vTree && (uTree->weight < vTree->weight))) {
			std::tie(u, uTree, v, vTree) = std::tuple(v, vTree, u, uTree);
		}
		ENode *e = newNode();
		e->setBoundary(u, v);
		e->userData = {userData...};

		/*
		std::cout << "rerooted uTree: ";
		debugPrintTree(uTree);
		std::cout << "new edge: ";
		debugPrintTree(e);
		std::cout << "rerooted vTree: ";
		debugPrintTree(vTree);
		*/

		ENode *tree = join<COMPRESS>(uTree, e, vTree);
		this->markNodeAsRoot(tree);
		tree = symmetrize(tree->boundary[1]);
		std::cout << "Linked: ";
		debugPrintTree(tree);
		//debugPrintRawTree(tree);
	}

	template <class... TUserData>
	std::tuple<TUserData...>  BiasedTreeTopTree<TUserData...>::
	cut(Vertex u, Vertex v) {
		this->rollback();
		assert(false);
	}
	

	// --- splice, reroot, symmetrize ---

	template <class... TUserData>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::
	reroot(Vertex v) {
		ENode *node = ext(this->vertexToNode[v]);
		if (!node || (v == node->boundary[0]) || (v == node->boundary[1])) {
			return node;
		}
		// std::cout << "reroot(" << v << ")" << std::endl;

		InnerList<Node, &Node::tmpListNode> splitNodes;
		if (node->clusterType == RAKE) {
			node = ext(node->children[1]);
			while (!node->template isLeaf<RAKE>()) {
				node = ext(node->children[0]);
			}
		} else {
			ENode *origNode = node;
			for (Node *origChild : origNode->children) {
				node = ext(origChild);
				while (!node->template isLeaf<COMPRESS>()) {
					node = ext(node->children[v == node->boundary[1]]);
				}
				if (node->isRakeTreeNode()) break;
				node = origNode;
			}
			splitNodes.push(node);
			node = findRoot<COMPRESS>(origNode, false);
		}

		while (node->parent) {
			if (splitNodes.front() != node) {
				splitNodes.push(node);
			}
			node = findRoot<RAKE>(node);
			if (splitNodes.front() != node) {
				splitNodes.push(node);
			}
			node = findRoot<COMPRESS>(node, true);
		}
		this->markNodeAsNonRoot(node);
		{
			ENode *prevNode = ext(splitNodes.front());
			if (splitNodes.front() != node) {
				splitNodes.push(node);
			}
			node = prevNode;
		}
		{
			bool rev = false; // relative to node
			int subpathsWeights[2] = {0, 0};
			if (!node->template isLeaf<COMPRESS>()) {
				subpathsWeights[0] = ext(node->children[0])->weight;
				subpathsWeights[1] = ext(node->children[1])->weight;
			}
			while (node->parent) {
				int childI = node->parent->children[1] == node;
				rev ^= ext(node->parent)->isChildReversed(childI);
				node = ext(node->parent);
				subpathsWeights[rev ^ !childI] += ext(node->children[!childI])->weight;
			}
			if (subpathsWeights[rev] > subpathsWeights[!rev]) {
				node->reverse();
			}
		}

		node = rerootRecurse(v, splitNodes);
		this->markNodeAsRoot(node);
		return node;
	}

	template <class... TUserData>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::
	rerootRecurse(Vertex v, InnerList<Node, &Node::tmpListNode> &splitNodes) {
		ENode *cRoot = ext(splitNodes.pop());
		ENode *cSplit = ext(splitNodes.pop());
		if (!cSplit) {
			if (!cRoot || (cRoot->boundary[0] == v) || (cRoot->boundary[1] == v)) {
				return cRoot;
			} else {
				assert(cRoot->getInnerVertex() == v);
				cSplit = cRoot; // the split node should be used again
			}
		} else if (cSplit->isRakeTreeNode() && cSplit->template isLeaf<RAKE>() && (cSplit->parent->boundary[1] != v)) {
			splitNodes.push(cSplit); // the split node should be used again
		}
		ENode *rSplit = ext(splitNodes.front());

		bool newBottomReversed = rSplit && (cSplit->parent->boundary[1] != rSplit->boundary[1]);
		auto [oldBottom, newBottom, rake, top] = spliceSplit(cRoot, cSplit, rSplit);
		if (newBottomReversed) newBottom->reverse(); // XXX maybe maintain elsewhere
		newBottom = rerootRecurse(v, splitNodes);
		return spliceJoin(newBottom, oldBottom, rake, top);
	}

	template <class... TUserData>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::
	symmetrize(Vertex v) {
		// std::cout << "Symmetrize(" << v << ")" << std::endl;
		ENode *root = ext(this->vertexToNode[v]);
		if (!root) return nullptr;
		assert(!root->parent);
		this->markNodeAsNonRoot(root);
		if (v != root->boundary[0]) {
			root->reverse();
		}

		root = symmetrizeRecurse(root);
		this->markNodeAsRoot(root);
		return root;
	}

	template <class... TUserData>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::
	symmetrizeRecurse(ENode *cRoot) {
		ENode *cSplit = this->template findHeavyNode<COMPRESS>(cRoot);
		if (!cSplit) return cRoot;
		ENode *rSplit = this->template findHeavyNode<RAKE>(cSplit);
		auto [oldBottom, newBottom, rake, top] = spliceSplit(cRoot, cSplit, rSplit);
		oldBottom = symmetrizeRecurse(oldBottom);
		return spliceJoin(newBottom, oldBottom, rake, top);
	}


	template <class... TUserData>
	inline std::tuple<typename BiasedTreeTopTree<TUserData...>::ENode *, typename BiasedTreeTopTree<TUserData...>::ENode *,typename BiasedTreeTopTree<TUserData...>::ENode *, typename BiasedTreeTopTree<TUserData...>::ENode *> BiasedTreeTopTree<TUserData...>::
	spliceSplit(ENode *cRoot, ENode *cSplit, ENode *rSplit) {
		auto [bottom, rake, top] = split<COMPRESS>(cRoot, cSplit);
		auto [rLeft, altBottom, rRight] = split<RAKE>(rake, rSplit);
		ENode *newRake = join<RAKE>(rLeft, nullptr, rRight);
		return {bottom, altBottom, newRake, top};
	}

	template <class... TUserData>
	inline typename BiasedTreeTopTree<TUserData...>::ENode *BiasedTreeTopTree<TUserData...>::
	spliceJoin(ENode *bottom, ENode *altBottom, ENode *rake, ENode *top) {
		auto [rLeft, rMiddle, rRight] = split<RAKE>(rake, altBottom);
		ENode *newRake = join<RAKE>(rLeft, altBottom, rRight);
		return join<COMPRESS>(bottom, newRake, top);
	}


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
		/*
		std::cout << "Join: " << std::endl;
		std::cout << "  ";
		debugPrintTree<TTreeType>(leftTree);
		std::cout << "  ";
		debugPrintTree<(TTreeType == RAKE) ? COMPRESS : RAKE>(middleNode);
		std::cout << "  ";
		debugPrintTree<TTreeType>(rightTree);
		*/
		// TODO reimplement as biased trees (currently unbalanced)
		// std::cout << (TTreeType == RAKE ? "rake_join" : "compress_join") << std::endl;

		if (middleNode) {
			middleNode->template setLeafData<TTreeType, true>();
			ClusterType nodeType = (TTreeType == RAKE) || middleNode->isRakeTreeNode() ? RAKE : COMPRESS;
			if (leftTree) {
				leftTree = joinNodes<TTreeType>(nodeType, leftTree, middleNode);
			} else if (rightTree) {
				if constexpr(TTreeType == RAKE) {
					rightTree = joinNodes<TTreeType>(nodeType, middleNode, rightTree); // preserve order
				} else {
					rightTree = joinNodes<TTreeType>(nodeType, rightTree, middleNode);
				}
			} else {
				return middleNode;
			}
		}

		ENode *root;
		if (leftTree == nullptr) {
			root = rightTree;
		} else if (rightTree == nullptr) {
			root = leftTree;
		} else {
			if ((TTreeType == COMPRESS) && leftTree->isRakeTreeNode()) {
				root = joinNodes<TTreeType>(RAKE, rightTree, leftTree);
			} else if ((TTreeType == COMPRESS) && rightTree->isRakeTreeNode()) {
				root = joinNodes<TTreeType>(RAKE, leftTree, rightTree);
			} else {
				root = joinNodes<TTreeType>(TTreeType, leftTree, rightTree);
			}
		}

		// debugPrintRawTree(root);
		assert(EIntegrity::treeConsistency(root));
		return root;
	}

	template <class... TUserData>
	template <ClusterType TTreeType>
	std::tuple<typename BiasedTreeTopTree<TUserData...>::ENode *, typename BiasedTreeTopTree<TUserData...>::ENode *, typename BiasedTreeTopTree<TUserData...>::ENode *> BiasedTreeTopTree<TUserData...>::
	split(ENode *tree, ENode *middleNode) {
		// TODO reimplement as biased trees (currently unbalanced)

		/*
		std::cout << (TTreeType == COMPRESS ? "Split<COMPRESS>: " : "Split<RAKE>") << std::endl;
		std::cout << "  ";
		debugPrintTree(tree);
		std::cout << "  ";
		debugPrintTree(middleNode);
		*/

		if (!middleNode) {
			return {tree, nullptr, nullptr};
		}
		if (!tree) {
			return {nullptr, nullptr, nullptr};
		}

		if constexpr(TTreeType == RAKE) {
			if ((middleNode != tree) && !middleNode->parent) {
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

		//InnerList<Node, &Node::tmpListNode> splitNodes;
		while (middleNode->parent) {
			middleNode->tmpMark = true;
			middleNode = ext(middleNode->parent);
		}
		assert(middleNode == tree);

		bool rev = false;
		ENode *left = nullptr, *right = nullptr;
		while ((middleNode->clusterType != BASE) && (middleNode->children[0]->tmpMark || middleNode->children[1]->tmpMark)) {
			assert(!middleNode->template isLeaf<TTreeType>());
			int nextI = middleNode->children[1]->tmpMark;
			bool childRev = rev ^ middleNode->isChildReversed(nextI);
			auto [leftChild, rightChild] = splitNode<TTreeType>(middleNode);
			ENode *children[2] = {leftChild, rightChild};
			children[nextI]->tmpMark = false;

			/*
			std::cout << "Left tree: ";
			debugPrintTree<TTreeType>(left);
			std::cout << "Right tree: ";
			debugPrintTree<TTreeType>(right);
			std::cout << "Child: ";
			debugPrintTree<TTreeType>(children[!nextI]);
			*/

			if (nextI ^ rev) {
				left = join<TTreeType>(left, nullptr, children[!nextI]);
			} else {
				right = join<TTreeType>(children[!nextI], nullptr, right);
			}
			rev = childRev;
			middleNode = children[nextI];
		}
		if (!middleNode->template isLeaf<TTreeType>()) {
			auto [leftChild, rightChild] = splitNode<TTreeType>(middleNode);
			ENode *children[2] = {leftChild, rightChild};
			left  = join<TTreeType>(left, nullptr,  children[rev]);
			right = join<TTreeType>(children[!rev], nullptr, right);
			middleNode = nullptr;
		} else {
			middleNode->template setLeafData<TTreeType, false>();
		}

		assert(EIntegrity::treeConsistency(left));
		assert(EIntegrity::treeConsistency(right));

		/*
		std::cout << "  =>";
		debugPrintTree(left);
		std::cout << "  =>";
		debugPrintTree(middleNode);
		std::cout << "  =>";
		debugPrintTree(right);
		*/

		return {left, middleNode, right};
	}

}

using TopTreeInternals::BiasedTreeTopTree;

#undef assert
#endif
