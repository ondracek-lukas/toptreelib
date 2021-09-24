// TopTreeLibrary  Copyright (C) 2020  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#ifndef CONTRACTION_TOP_TREE_HPP
#define CONTRACTION_TOP_TREE_HPP

#include "TopTree.hpp"
#include <cstdint>

#ifdef TOP_TREE_INTEGRITY
#include "tests/ContractionTopTreeIntegrity.hpp"
#else
namespace TopTreeInternals {
	template <class... TUserData> class ContractionTopTreeIntegrity {};
}
#define assert(cond)
#endif

namespace TopTreeInternals {

	template <class... TUserData>
	class ContractionTopTree : public TopTree<TUserData...> {
		protected:
			using Integrity = typename TopTree<TUserData...>::Integrity;
			using EIntegrity = ContractionTopTreeIntegrity<TUserData...>; friend EIntegrity;
#ifdef TOP_TREE_INTEGRITY
		public:
			virtual void testIntegrity() { // just root tree test
				Integrity::treeConsistency(ext(this->exposedRoot));
				this->rollback();
				EIntegrity::treeConsistency(ext(this->exposedRoot));
			}
#endif

			virtual void link(Vertex u, Vertex v, TUserData... userData);
			virtual std::tuple<TUserData...> cut(Vertex u, Vertex v);

		private:


			using Node = typename TopTree<TUserData...>::Node;

			struct VirtNode;

			struct Arc { // circular dequeue, always stored in VirtNode.arcs array with its twin
				enum Dir {PRED = 0, SUCC = 1};
				private:
					friend EIntegrity;
					Arc *ptrs[2]; // pred, succ
				public:
					VirtNode *const virtNode;

					size_t arcIndex() { return virtNode->arcs != this; }
					Arc *pred() { return ptrs[PRED]; }
					Arc *succ() { return ptrs[SUCC]; }
					Arc *get(int dir) { return ptrs[dir]; }
					Arc *twin()  { return &virtNode->arcs[!arcIndex()]; }

					Arc(VirtNode *virtNode) : virtNode(virtNode) {
						ptrs[0] = ptrs[1] = twin();
					}

					// Initializes arc as leaf
					void initLeaf() {
						ptrs[SUCC] = twin();
						ptrs[SUCC]->ptrs[PRED] = this;
					}

					// Connects (possibly uninitialized) this to tour given by its new succ arc
					void connectToTour(Arc *arc) {
						ptrs[SUCC] = arc;
						twin()->ptrs[PRED] = arc->ptrs[PRED];
						arc->ptrs[PRED]->ptrs[SUCC] = twin();
						arc->ptrs[PRED] = this;
					}

					// Disconnects this from its succ leaving this->ptrs inconsistently unchanged
					void disconnectFromTour() {
						ptrs[SUCC]->ptrs[PRED] = twin()->ptrs[PRED];
						ptrs[SUCC]->ptrs[PRED]->ptrs[SUCC] = ptrs[SUCC];
					}
					
					// Sets pred/succ to given arc on lower level, marking the pointer
					template <Dir D>
					void setMarked(Arc *orig) { // call upgradeIfMarked() before using ptrs[D]
						ptrs[D] = (Arc *) ((uintptr_t) orig | 1);
					}

					// Sets marked pred/succ (to lower level) to the corresponding arc on proper level, unmarking it
					template <Dir D>
					void upgradeIfMarked() {
						if ((uintptr_t) ptrs[D] & 1) {
							ptrs[D] = (Arc *) ((uintptr_t) ptrs[D] & ~(uintptr_t)1);
							assert(ptrs[D]->virtNode->virtParent);
							assert(virtNode->node->boundary[D^arcIndex()] == ptrs[D]->virtNode->node->boundary[!D^ptrs[D]->arcIndex()]);
							assert(
									(virtNode->node->boundary[D^arcIndex()] == ptrs[D]->virtNode->virtParent->node->boundary[0]) ||
									(virtNode->node->boundary[D^arcIndex()] == ptrs[D]->virtNode->virtParent->node->boundary[1]));

							ptrs[D] = &ptrs[D]->virtNode->virtParent->arcs[
								virtNode->node->boundary[D^arcIndex()] == ptrs[D]->virtNode->virtParent->node->boundary[D]];
							ptrs[D]->ptrs[!D] = this;
							assert(EIntegrity::template arcConsistency<D>(this));
						}
					}

					// Determines whether eNode.boundary[i] is leaf where eNode.virtNode->arcs[!i] is this
					bool isLeaf() {
						return virtNode == succ()->virtNode;
					}
					bool canCompress() {
						assert(!isLeaf());
						return succ()->twin()->succ()->virtNode == virtNode;
					}
					bool canRake() {
						assert(!isLeaf());
						return succ()->isLeaf();
					}
			};

			struct ENode;

			struct VirtNode {
				Arc arcs[2] = {this, this};
					// inv: arcs[i] leads from node->boundary[i];
					// inv: node was created through virt(node->children[0])->arc[...]->succ() pointer
				ENode *node; // same as in child if dummy
				VirtNode *virtParent = nullptr; // may refer to dummy, which is skipped by node->parent
				Arc *childrenArcs[2] = {nullptr, nullptr}; // the arc through which node was created and its successor

				private:
					static inline VirtNode * const NOT_IN_LIST = (VirtNode * const) ~(uintptr_t)0;
					VirtNode *listNode = NOT_IN_LIST;
				public:
					class List : private InnerList<VirtNode, &VirtNode::listNode> {
						private:
							using base = InnerList<VirtNode, &VirtNode::listNode>;
						public:
							using base::base;
							using base::front;
							using base::next;
							using base::begin;
							using base::end;
							void append(VirtNode *node) {
								assert(node != NOT_IN_LIST);
								assert(node->listNode == NOT_IN_LIST);
								node->listNode = NULL;
								base::append(node);
							}
							bool tryAppend(VirtNode *node) {
								assert(node != NOT_IN_LIST);
								if (node->listNode != NOT_IN_LIST) {
									return false;
								}
								append(node);
								return true;
							}
							VirtNode *pop() {
								VirtNode *node = base::pop();
								assert(node != NOT_IN_LIST);
								if (node) node->listNode = NOT_IN_LIST;
								return node;
							}
					};

					bool isDummy() {
						return this != virt(node); // compares just pointers, should work even after deleting node
					}
					bool isMoveValid() {
						assert(!isDummy());
						return
							childrenArcs[0] && (childrenArcs[0]->succ() == childrenArcs[1]) &&
							(node->clusterType == COMPRESS ? childrenArcs[0]->canCompress() : childrenArcs[0]->canRake());
					}

					// Gets neighbouring arc on lower level incident to index'th arc as pred/succ
					Arc *lowerLevelArcs(int index, typename Arc::Dir direction) {
						assert(isDummy() || (node->clusterType != BASE));
						Arc *arc = childrenArcs[0];
						if (index == 0) {
							if ((direction == Arc::Dir::SUCC) && !isDummy()) {
								arc = arc->succ();
								if (node->clusterType == RAKE) arc = arc->succ();
							}
						} else {
							arc = arc->twin();
							if ((direction == Arc::Dir::PRED) && !isDummy() && (node->clusterType == COMPRESS)) {
								arc = arc->pred();
							}
						}
						return arc->get(direction);
					}
			};

			struct ENode : public Node {
				VirtNode virtNode;
				ENode() { virtNode.node = this; }
			};

			static ENode *ext(Node *node) { return static_cast<ENode *>(node); };
			static VirtNode *virt(Node *node) { return &ext(node)->virtNode; };

			InnerList<Node, &Node::tmpListNode> freeENodes;
			ENode *newNode() {
				ENode *node = this->template newNodeT<ENode>(freeENodes);
				return node;
			}
			void deleteNode(Node *node) { this->template deleteNodeT<ENode>(node, freeENodes); }

			InnerList<VirtNode, &VirtNode::virtParent> freeDummyNodes;
			VirtNode *newDummyNode() {
				if (freeDummyNodes.front()) {
					return freeDummyNodes.pop();
				}
				return new VirtNode();
			}
			void deleteDummyNode(VirtNode *node) {
				assert(!node->virtParent);
				freeDummyNodes.push(node);
			}

			std::vector<Arc *> vertexToArc; // any arc leading from the vertex, or nullptr

			VirtNode *tryJoin(Arc *arc);

			void updateClusterization(typename VirtNode::List &&I, typename VirtNode::List &&D);
	};


	template <class... TUserData>
	void ContractionTopTree<TUserData...>::
	link(Vertex u, Vertex v, TUserData... userData) {
		this->rollback();
		ENode *node = newNode();
		node->setBoundary(u, v);
		node->userData = {userData...};
		if (vertexToArc.size() <= std::max(u, v)) {
			vertexToArc.resize(std::max(u, v) + 1);
		}
		for (int i : {0, 1}) {
			if (Arc *arc = vertexToArc[node->boundary[i]]) {
				node->virtNode.arcs[!i].connectToTour(arc);
			} else {
				node->virtNode.arcs[!i].initLeaf();
				vertexToArc[node->boundary[i]] = &node->virtNode.arcs[i];
			}
		}
		typename VirtNode::List I;
		I.append(&node->virtNode);
		updateClusterization(std::move(I), {});

	}

	template <class... TUserData>
	std::tuple<TUserData...>  ContractionTopTree<TUserData...>::
	cut(Vertex u, Vertex v) {
		this->rollback();
		ENode *node = ext(TopTree<TUserData...>::baseNode(u, v));
		assert(node); // interface assert, maybe use exception instead
		typename VirtNode::List D;
		D.append(&node->virtNode);
		for (int i : {0, 1}) {
			if (node->virtNode.arcs[!i].isLeaf()) {
				this->markVertexAsIsolated(node->boundary[i]);
				assert(vertexToArc[node->boundary[i]] == &node->virtNode.arcs[i]);
				vertexToArc[node->boundary[i]] = nullptr;
			} else {
				if (vertexToArc[node->boundary[i]] == &node->virtNode.arcs[i]) {
					vertexToArc[node->boundary[i]] = node->virtNode.arcs[!i].succ();
				}
				node->virtNode.arcs[!i].disconnectFromTour();
			}
		}
		if (node->parent) {
			this->releaseNode(node->parent);
		}
		std::tuple<TUserData...> ret = std::move(node->userData);
		updateClusterization({}, std::move(D));
		return ret;
	}

	// Updates clusterization given inserted and deleted nodes on base level, see state below
	template <class... TUserData>
	void ContractionTopTree<TUserData...>::
	updateClusterization(typename VirtNode::List &&I, typename VirtNode::List &&D) {
		while (I.front() || D.front()) {
			// state:
			//   level of nodes in I,D has complete Eulerian trail (incl. I, excl. D)
			//   I-nodes have no parents
			//   D-nodes have parents and should be destroyed, trail ptrs inconsistently points to original neighbours

			typename VirtNode::List newI, newD, I2, N;
				// lists:
				//   D          nodes to be destroyed
				//   I          newly inserted nodes without parents, need creating new parents, neighbours may have invalid parents
				//   I2         nodes disconnected from their invalidated parents, need creating new parents
				//   N          neighbouring nodes with valid parents, but dummy parents may need replacement
				//   newI,newD  lists of by one higher level
				// Any node may be in at most one list at a time --
				// trying to append the node to a second list has no effect;
				// the order in which lists are created is thus very important.


			// ADD siblings of D to I2, ADD parents to newD, DETACH from parents
			for (VirtNode *node : D) {
				VirtNode *parent = node->virtParent;
				if (parent) {
					if (parent->isDummy()) {
						node->virtParent = nullptr;
						parent->childrenArcs[0] = nullptr;
					} else {
						size_t i = parent->childrenArcs[0]->virtNode == node;
						I2.tryAppend(node->virtParent->childrenArcs[i]->virtNode);
						for (Arc *&arc : parent->childrenArcs) {
							arc->virtNode->virtParent = nullptr;
							arc = nullptr;
						}
						parent->node->detachChildren();  // needed if node isDummy and node->node still exists
					}
					newD.tryAppend(parent);
				}
				if (!node->isDummy() && node->node->parent) {
					node->node->parent->detachChildren();
				}
			}

			// TEST VALIDITY of parents of arc-neighbours of I, insert invalid to newD
			// insert both children of invalid parent to I2, others to N
			for (VirtNode *node2 : I) {
				for (int i : {0, 1}) for (int d : {0, 1}) {
					VirtNode *node = node2->arcs[i].get(d)->virtNode;
					VirtNode *parent = node->virtParent;
					if (parent && parent->isDummy()) {
						N.tryAppend(node);
					} else if (!parent) {
						I2.tryAppend(node);
					} else if (!parent->isMoveValid()) {
						this->releaseNode(parent->node); // no effect if second child is already destroyed
						for (Arc *&arc : parent->childrenArcs) {
							I2.tryAppend(arc->virtNode);
							arc->virtNode->virtParent = nullptr;
							arc = nullptr;
						}
						if (node->node->parent) {
							node->node->parent->detachChildren();
						}
						newD.tryAppend(parent);
					}
				}
			}

			// ADD neighbours of D and their neighbours to N
			for (VirtNode *node : D) {
				for (int i : {0, 1}) for (int d : {0, 1}) {
					VirtNode *node2 = node->arcs[i].get(d)->virtNode;
					N.tryAppend(node2);
					for (int i : {0, 1}) for (int d : {0, 1}) {
						N.tryAppend(node2->arcs[i].get(d)->virtNode);
					}
				}
			}

			// ADD neighbours of I2 to N
			for (VirtNode *node : I2) {
				for (int i : {0, 1}) for (int d : {0, 1}) {
					N.tryAppend(node->arcs[i].get(d)->virtNode);
				}
			}


			// DELETE nodes in D
			while (VirtNode *node = D.pop()) {
				if (node->isDummy()) {
					deleteDummyNode(node);
				} else {
					deleteNode(node->node);
				}
			}
			
			// CREATE new PARENTS of I,I2,N; insert them to newI
			// temporary set their arc-ptrs to lower level nodes
			// insert unneeded dummy nodes to newD
			for (auto X : {&I, &I2, &N}) for (VirtNode *node : *X) {
				for (int i : {0, 1}) {
					VirtNode *newParent = tryJoin(&node->arcs[i]);
					if (newParent) {
						VirtNode *sibling = node->arcs[i].succ()->virtNode; // == newParent->childrenArcs[1]->virtNode
						for (VirtNode *n : {node, sibling}) {
							if (n->virtParent) {
								assert(n->virtParent->isDummy());
								n->virtParent->childrenArcs[0] = nullptr;
								newD.append(n->virtParent);
							}
							n->virtParent = newParent;
						}
						newI.append(newParent);
						assert(EIntegrity::childrenConsistency(newParent));
						break;
					}
				}
			}
			while (N.pop());
			
			// CREATE new DUMMY parents of nodes in I,I2 (unless neighbourless); insert them to newI
			// copy even their arc-ptrs (temporary)
			for (auto X : {&I, &I2}) while (VirtNode *node = (*X).pop()) {
				assert(EIntegrity::childrenNeighbourhoodConsistency(node));
				if (!node->node->parent) {
					if (node->arcs[0].isLeaf() && node->arcs[1].isLeaf()) {
						this->markNodeAsRoot(node->node);
						if (node->virtParent) {
							newD.append(node->virtParent);
							node->virtParent->childrenArcs[0] = nullptr;
							node->virtParent = nullptr;
						}
						assert(EIntegrity::treeConsistency(node->node));
					} else if (!node->virtParent) {
						VirtNode *newParent = newDummyNode();
						newParent->childrenArcs[0] = &node->arcs[0];
						node->virtParent = newParent;
						newParent->node = node->node;
						for (int i : {0, 1}) {
							newParent->arcs[i].template setMarked<Arc::Dir::PRED>(node->arcs[i].pred());
							newParent->arcs[i].template setMarked<Arc::Dir::SUCC>(node->arcs[i].succ());
						}
						newI.append(newParent);
						assert(EIntegrity::childrenConsistency(newParent));
					}
				}
			}

			// DISCONNECT newD nodes from tour
			for (VirtNode *node : newD) {
				for (int i : {0, 1}) {
					node->arcs[i].disconnectFromTour();
				}
			}
			
			// UPGRADE ARC pointers in newI to higher level if needed
			for (VirtNode *node : newI) {
				for (Arc &arc : node->arcs) {
					arc.template upgradeIfMarked<Arc::Dir::PRED>();
					arc.template upgradeIfMarked<Arc::Dir::SUCC>();
				}
				assert(EIntegrity::nodeConsistency(node));
			}

			I = std::move(newI);
			D = std::move(newD);

			assert(!newI.front());
			assert(!newD.front());
			assert(!I2.front());
			assert(!N.front());
		}

	}

	template <class... TUserData>
	typename ContractionTopTree<TUserData...>::VirtNode *ContractionTopTree<TUserData...>::
	tryJoin(Arc *arc) {
		if (
				arc->isLeaf() ||
				(arc->virtNode->virtParent && !arc->virtNode->virtParent->isDummy()) ||
				(arc->succ()->virtNode->virtParent && !arc->succ()->virtNode->virtParent->isDummy())) {
			return nullptr;
		}

		ClusterType clusterType;
		if (arc->canCompress()) {
			clusterType = COMPRESS;
		} else if (arc->canRake()) {
			clusterType = RAKE;
		} else return nullptr;

		for (Arc *arc : {arc, arc->succ()}) {
			if (arc->virtNode->node->parent) {
				assert(arc->virtNode->virtParent);
				this->releaseNode(arc->virtNode->node->parent);
				arc->virtNode->node->detach();
			}
		}

		ENode *node = newNode();
		node->attachChildren(clusterType, arc->virtNode->node, arc->succ()->virtNode->node);
		this->vertexToNodeUpdateInner(node);
		node->virtNode.childrenArcs[0] = arc;
		node->virtNode.childrenArcs[1] = arc->succ();

		VirtNode *virtNode = &node->virtNode;

		for (int i : {0, 1}) {
			virtNode->arcs[i].template setMarked<Arc::Dir::PRED>(virtNode->lowerLevelArcs(i, Arc::Dir::PRED));
			virtNode->arcs[i].template setMarked<Arc::Dir::SUCC>(virtNode->lowerLevelArcs(i, Arc::Dir::SUCC));
		}

		assert(virtNode->node == node);

		return virtNode;
	}
}

using TopTreeInternals::ContractionTopTree;

#undef assert
#endif
