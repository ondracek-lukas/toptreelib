// TopTreeLibrary  Copyright (C) 2022  Lukáš Ondráček <ondracek.lukas@gmail.com>, use under MIT license

/* The Top Tree base abstract class.
 * It should be derived by specific drivers
 * and specialized by user data via templates.
 */

#ifndef TOP_TREE_HPP
#define TOP_TREE_HPP

#if __cplusplus < 201703L
#error At least C++17 needed.
#endif

namespace TopTreeInternals {
	enum ClusterType { BASE, COMPRESS, RAKE };  // aka ::TopTreeClusterType
	typedef int Vertex; // aka ::TopTreeVertex, vertex of underlying tree
}

#ifdef TOP_TREE_INTEGRITY
#include "tests/TopTreeIntegrity.hpp"
#else
namespace TopTreeInternals {
	template <class... TUserData> class TopTreeIntegrity {};
}
#define assert(cond)
#endif

#include <tuple>
#include <vector>
#include <cstddef>

#include "TopTreeInternals/SubtreeTraversability.hpp"
#include "TopTreeInternals/InnerList.hpp"

namespace TopTreeInternals {

	template <class TUserData>
	struct EventData {  // aka ::TopTreeEventData
		private:
			struct DataRefsArray {
				private: TUserData *elems[2];
				public:
					DataRefsArray(TUserData &elem1, TUserData &elem2) : elems{&elem1, &elem2} {};
					TUserData &operator[](int i) {return *elems[i]; }
			};
		public:
			ClusterType type;
			TUserData &parent;
			DataRefsArray children;
			Vertex boundary[2];
			Vertex innerVertex;
			Vertex childrenBoundary[2][2];
			EventData(
				ClusterType type,
				TUserData &parent, TUserData &child1, TUserData &child2,
				Vertex *boundary, Vertex innerVertex, Vertex *child1Boundary, Vertex* child2Boundary) :
				type(type), parent(parent), children{child1, child2},
				boundary{boundary[0], boundary[1]}, innerVertex(innerVertex),
				childrenBoundary{{child1Boundary[0], child1Boundary[1]},{child2Boundary[0], child2Boundary[1]}} {}
	};

	template <class... TUserData>
	class TopTree {  // aka ::TopTree
		protected:
			using Integrity = TopTreeIntegrity<TUserData...>; friend Integrity;

		public:
			virtual bool expose(Vertex u, Vertex v);
			virtual bool exposeTree(Vertex v); // sets tree with v as exposed, returns false if v is isolated and no such tree exists
			virtual void link(Vertex u, Vertex v, TUserData... userData) = 0; // use move semantics for userData? XXX
			virtual std::tuple<TUserData...> cut(Vertex u, Vertex v) = 0;

			template <int I = 0, bool onExposedPath = false, class TSelect>
			void search(TSelect select); // not virtual; redeclare using std:function if needed

			template <int I = 0, class TSelect>
			void pathSearch(TSelect select) { search<I, true>(select); };


			Vertex newVertex() {
				vertexToNode.emplace_back(nullptr);
				return vertexToNode.size()-1;
			};

			int getVerticesCnt() {
				return vertexToNode.size();
			}

			// Direct connectivity queries may be removed in future versions
			// in favour of expose.
			bool sameComponent(Vertex u, Vertex v);

			// Returns reference to I-th TUserData of root node.
			// Returned root data may be accessed and modified
			// until next call to any method.
			template <int I = 0>
			auto &getExposedData();

			// Getting edge data should be used only if we know
			// that current boundary vertices are connected by an edge.
			// It is a shortcut for performing cut and link on that edge
			// and may be removed in future versions.
			template <int I = 0>
			auto &getEdgeData();

			std::pair<Vertex, Vertex> getBoundary();
			Vertex getBoundary(int index);


		protected:  // to be accessible by driver

			struct Node;

			struct NodeInPath {
				Node *node = nullptr;
				bool onPath = true;
			};

			struct Node : SubtreeTraversability<Node> {  // node in top tree
				size_t index; // < nodesAllocated, to be used in vectors with extension data; maybe not needed? XXX
				Vertex boundary[2];
				Node *children[2] = {nullptr, nullptr};
				Node *parent = nullptr;
				// constraints:
					// boundary[i] should be in children[i]->boundary;
					// children[0] is on path even in RAKE node
				ClusterType clusterType;
				std::tuple<TUserData...> userData;
				bool userDataValid[sizeof...(TUserData)] = {};
				Node *tmpListNode = nullptr;
				Node *rollbackListNode = nullptr;
				bool tmpMark = false;

				Vertex getOtherVertex(Vertex vertex) {
					assert((boundary[0] == vertex) || (boundary[1] == vertex));
					return boundary[ boundary[0] == vertex ];
				}

				Vertex getInnerVertex() { // or maybe store directly?
					assert(Integrity::nodeChildrenBoundary(this));
					return children[1]->getOtherVertex(boundary[1]);
				}

				Vertex getSharedVertex() {
					assert(Integrity::nodeChildrenBoundary(this));
					return children[0]->getOtherVertex(boundary[0]);
				}

				std::tuple<NodeInPath, NodeInPath> getChildren(Vertex firstVertex = nullptr, bool rev = false);

				template <size_t I>
				EventData<typename std::tuple_element<I, std::tuple<TUserData...>>::type> getEventData();

				template <size_t I>
				void userDataSplit();

				template <size_t I>
				void userDataJoin();

				// For base cluster use setBoundary;
				// for non-base cluster use either attachChildren, or twice attachChild followed by setBoundary.

				Node *attachChildren(ClusterType type, Node *child1, Node *child2);
				Node *attachChild(int childIndex, Node* child);
				Node *detachChildren();
				Node *detach();

				void setBoundary(ClusterType type); // for COMPRESS or RAKE node
				void setBoundary(Vertex v1, Vertex v2); // for BASE node

			};

			Node *treeRoot(Node *node) {
				if (!node) return nullptr;
				while (node->parent) node = node->parent;
				return node;
			}

			Node *baseNode(Vertex u, Vertex v) {
				Node *node = this->vertexToNode[v];
				if (!node) return nullptr;
				if ((node->boundary[0] != u) && (node->boundary[1] != u)) {
					std::swap(u, v);
					node = this->vertexToNode[v];
					if (!(node && ((node->boundary[0] == u) || (node->boundary[1] == u)))) return nullptr;; // interface assert
				}
				if (node->clusterType != BASE) {
					node = node->children[node->boundary[1] == u];
					while (node->clusterType != BASE) {
						if (node->clusterType != RAKE) return nullptr;
						node = node->children[0];
					}
				}
				assert(
						((node->boundary[0] == u) && (node->boundary[1] == v)) ||
						((node->boundary[1] == u) && (node->boundary[0] == v)));
				return node;
			}

			// Nodes of the original tree with child used in the temporary tree;
			// the child's parent pointer points into the temporary tree;
			// original roots should have tmpMark set.
			InnerList<Node, &Node::rollbackListNode> rollbackList;

			// To be called after implicit expose or search;
			// can be extended by same actions for amortization in non-worst-case version.
			// Returns whether anything has been unrolled.
			virtual bool rollback();

			Node *exposedRoot; // root of the exposed tree

			std::vector<Node*> vertexToNode; // the node having given vertex as inner vertex; root node otherwise
				// XXX or use struct Vertex instead?


			// To correctly handle user data when tree is being changed,
			// driver should always call releaseNode before changing its children.

			// Calls split on node and all ancestors with valid user data top-down
			// marking them invalid.
			void releaseNode(Node *node);

			// Calls split only on given node assuming invalid data in ancestors.
			void releaseJustNode(Node *node);

			// As releaseNode but only for user data of indices given by template arguments.
			template <size_t I, size_t... Is>
			void releaseNodeData(Node *node, std::index_sequence<I, Is...> = std::index_sequence<I, Is...>());

			// As releaseJustNode but only for user data of indices given by template arguments.
			template <std::size_t... Is>
			void releaseJustNodeData(Node *node, std::index_sequence<Is...> = std::index_sequence<Is...>());

			// Calls join on all descendants with invalid user data bottom-up
			// marking them valid, only for user data of index given by template argument.
			template <size_t I>
			void validateUserData(Node *root, InnerList<Node, &Node::rollbackListNode> &splitList);
			template <size_t I>
			void validateUserData(Node *root) { validateUserData<I>(root, rollbackList); }


			// -- tree structure manipulation --

			// Methods newNode() and deleteNode(Node) can be redefined (masked) in driver
			// to allow allocation of derived node instances;
			// the templated counterparts should be called from them.
			// All nodes created by the base class are deleted in rollback method.

			size_t nodesAllocated = 0;
			InnerList<Node, &Node::tmpListNode> freeNodes;

			template <class TNode>
			TNode *newNodeT(InnerList<Node, &Node::tmpListNode> &freeTNodes) {
				TNode *node = static_cast<TNode *>(freeTNodes.pop());
				if (!node) {
					node = new TNode;
					node->index = nodesAllocated++;
				}
				return node;
			}
			Node *newNode() { return newNodeT<Node>(freeNodes); }

			template <class TNode> // templating may be used later
			void deleteNodeT(Node *node, InnerList<Node, &Node::tmpListNode> &freeTNodes) {
				assert(node && !node->parent && !node->children[0] && !node->children[1]);
				freeTNodes.push(node);
			}
			void deleteNode(Node *node) { deleteNodeT<Node>(node, freeNodes); }

			void vertexToNodeUpdateInner(Node *node) {
				vertexToNode[node->getInnerVertex()] = node;
			}

			void markVertexAsIsolated(Vertex vertex) {
				vertexToNode[vertex] = nullptr;
			}

			void markNodeAsRoot(Node *node) {
				vertexToNode[node->boundary[0]] = node;
				vertexToNode[node->boundary[1]] = node;
				exposedRoot = node;
			}

			void markNodeAsNonRoot(Node *node) { // or mark isolated vertices instead
				vertexToNode[node->boundary[0]] = nullptr;
				vertexToNode[node->boundary[1]] = nullptr;
			}

		private:
			NodeInPath joinNodesInPath(NodeInPath np1, NodeInPath np2);

		public:
#ifdef TOP_TREE_INTEGRITY
			virtual void testIntegrity() { // XXX just root tree test
				assert(Integrity::treeConsistency(exposedRoot));
			}
#endif

			virtual ~TopTree() {
				while (Node *node = freeNodes.pop()) delete node;
			}

	};



	// === TopTree methods' implementation ===

	// --- default implementation of top tree structure maintenance  ---

	template <class... TUserData>
	inline typename TopTree<TUserData...>::NodeInPath TopTree<TUserData...>::joinNodesInPath(NodeInPath np1, NodeInPath np2) {
		if (np1.node && np2.node) {
			NodeInPath nRet = { .node = newNode(), .onPath = true };
			if (np1.onPath) {
				nRet.node->attachChildren(np2.onPath ? COMPRESS : RAKE, np1.node, np2.node);
			} else {
				nRet.node->attachChildren(RAKE, np2.node, np1.node);
				nRet.onPath = np2.onPath;
			}
			return nRet;
		} else if (np1.node) {
			return np1;
		} else {
			return np2;
		}
	}

	template <class... TUserData>
	bool TopTree<TUserData...>::expose(Vertex u, Vertex v) {
		rollback();
		if (!sameComponent(u, v) || !vertexToNode[u]) return false;

		Vertex vL, vPostpR;
		NodeInPath npL, npR, npM, npPostpM, npPostpR;

		// mark all predecessors of u and v, to be processed
		for (Vertex w : {u, v}) {
			Node *n = vertexToNode[w];
			assert(n);
			if ((n->boundary[0] == w) || (n->boundary[1] == w)) continue;
			while (!n->tmpMark) {
				n->tmpMark = true;
				if (!n->parent) {
					assert(!npM.node); // otherwise u, v are in different trees
					npM.node = n;
					break;
				}
				n = n->parent;
			}
		}

		if (!npM.node) {
			exposedRoot = vertexToNode[u];
			return true;
		}

		/* We are now maintaining a path composed of clusters (and vertices) in the order:
		 *   npL, vpL, npM = (npML, npMR), npR, npPostpM, vPostpR, npPostpR;
		 * some of np nodes may be raked to the path (if !.onPath), i.e. sharing only one vertex.
		 * We are decomposing a middle node npM into
		 *   its unmarked descendants to be added to left side npL or right side npR,
		 *   its marked descendants to become new npM, or to be postponed to npPostpM;
		 * other joining of clusters in path may occur,
		 * usually keeping the right boundary vertex of npR unchanged --
		 * having found one of u, v, it is kept there.
		 */

		vL = npM.node->boundary[0];

		while (npM.node) {
			assert(npM.node->clusterType != BASE);
			npM.node->tmpMark = false;
			rollbackList.append(npM.node);

			if ((vL == u) || (vL == v)) {
				vL = npM.node->getOtherVertex(vL);
				std::tie(npL, npR) = std::tuple(npR, npL);
			}

			auto [npML, npMR] = npM.node->getChildren(vL);
			Vertex npM_inner = npM.node->getInnerVertex();

			if ((npM_inner != u) && (npM_inner != v)) {

				// join unmarked subcluster to side, or postpone right subcluster if both marked
				if (!npML.node->tmpMark) {
					assert(npMR.node->tmpMark);
					npL = joinNodesInPath(npL, npML);
					vL = npM.node->getSharedVertex();
					npM = npMR;
				} else if (!npMR.node->tmpMark) {
					npR = joinNodesInPath(npR, npMR);
					npM = npML;
				} else {
					npPostpR = npR;  npR = {};
					vPostpR = npM.node->getOtherVertex(vL);
					npPostpM = npMR;
					npM = npML;
				}
			} else {
				if (!npML.node->tmpMark && !npMR.node->tmpMark) {

					// join npL, npML, npMR to npR keeping the inner vertex in boundary
					if (npML.onPath) {
						npL = joinNodesInPath(npL, npML);
						npML = {};
						npM = npMR;
					}
					if (npMR.onPath) {
						npR = joinNodesInPath(npR, npMR);
						npMR = {};
						npM = npML;
					}
					npL.onPath = false;
					npR = joinNodesInPath(npR, npL);
					npL = {};
					npM.onPath = true;
					npR = joinNodesInPath(npR, npM);

					// restore postponed if any
					npM = npPostpM; npPostpM = {};
					vL = vPostpR;
					npL = npPostpR;

				} else {
					assert(npML.node->tmpMark ^ npMR.node->tmpMark);
					assert(!npPostpM.node);
					if (npMR.node->tmpMark) {
						// reverse path
						std::tie(npL, npML, npMR, npR) = std::tuple(npR, npMR, npML, npL);
						vL = npM.node->getOtherVertex(vL);
					}

					// join npMR to npR keeping the inner vertex in boundary
					bool npMR_origOnPath = npMR.onPath;
					npR.onPath = npMR.onPath;
					npMR.onPath = true;
					npR = joinNodesInPath(npR, npMR);
					npR.onPath = !npMR_origOnPath;
					npM = npML;
				}
			}


			if (!npM.onPath) {
				assert(npM.node);

				// make npM onPath
				npL.onPath = false;
				npR = joinNodesInPath(npR, npL);
				npL = {};
				npM.onPath = true;
				vL = npM.node->getOtherVertex(vL);
			}
		}
		rollbackList.back()->tmpMark=true;
		exposedRoot = npR.node;

		assert(Integrity::treeConsistency(exposedRoot));
		return true;
	};

	template <class... TUserData>
	bool TopTree<TUserData...>::exposeTree(Vertex v) {
		rollback();
		if (!vertexToNode[v]) return false;
		exposedRoot = treeRoot(vertexToNode[v]);
		return true;
	};

	template <class... TUserData>
	template <int I, bool onExposedPath, class TSelect>
	inline void TopTree<TUserData...>::search(TSelect select) {
		if (!exposedRoot || (exposedRoot->clusterType == BASE)) {
			return;
		}

		if (!onExposedPath) {
			rollback();
		}

		InnerList<Node, &Node::rollbackListNode> tmpRollbackList = std::move(rollbackList);
		Node *tmpRollbackListCurNode = nullptr;
		Node *rollbackListCurNode = nullptr;


		// We maintain path: npL = (npLL, npRR), vM, npR = (npRL, npRR)
		// Inv.: npL either equals npLR, or is temporary; symmetrically for npR, npRL
		NodeInPath npL, npLL, npLR, npR, npRL, npRR;
		Vertex vM;

		npL.node = npLR.node = exposedRoot;
		vM = npLR.node->boundary[0];

		while (npLR.node->clusterType != BASE) {

			// npL, npR are not set at this point; npRL is empty; npLR should be split

			tmpRollbackList.insertAfter(tmpRollbackListCurNode, npLR.node);
			tmpRollbackListCurNode = npLR.node;

			std::tie(npLR, npRL) = npLR.node->getChildren(vM, true);
			if (npRL.onPath) {
				vM = npRL.node->getOtherVertex(vM);
			}

			npL = joinNodesInPath(npLL, npLR);
			npR = joinNodesInPath(npRL, npRR);

			bool reverse;
			if (onExposedPath && (!npLR.onPath || !npRL.onPath)) {
				reverse = !npLR.onPath;
			} else {
				NodeInPath npRoot = joinNodesInPath(npL, npR);
				assert(Integrity::treeConsistency(npRoot.node));
				validateUserData<I>(npRoot.node, tmpRollbackList);

				rollbackList.insertAfter(rollbackListCurNode, std::move(tmpRollbackList));
				rollbackListCurNode = tmpRollbackListCurNode;
				tmpRollbackListCurNode = nullptr;

				exposedRoot = npRoot.node;
				if (rollbackListCurNode) rollbackListCurNode->tmpMark = true; // clean up for the case of exception in select

				reverse = select(npRoot.node->template getEventData<I>());
				reverse ^= npRoot.node->children[0] != npL.node;

				if (rollbackListCurNode) rollbackListCurNode->tmpMark = false;

				releaseJustNode(npRoot.node);
				npRoot.node->detachChildren();
				deleteNode(npRoot.node);
			}

			if (reverse) {
				std::tie(npLL, npLR, npL, npR, npRL, npRR) = std::tuple(npRR, npRL, npR, npL, npLR, npLL);
			}

			if (npL.node != npLR.node) {
				releaseJustNode(npL.node);
				npL.node->detachChildren();
				deleteNode(npL.node);
			}
			npL = {};
			npRR = npR;
			npR = {}; npRL = {};

			if (!npLR.onPath) {
				assert(!onExposedPath);
				npLL.onPath = false;
				npRR = joinNodesInPath(npRR, npLL);
				npLL = {};
				npLR.onPath = true;
			}
		}

		if (tmpRollbackList.front()) {
			rollbackList.insertAfter(rollbackListCurNode, std::move(tmpRollbackList));
			rollbackListCurNode = tmpRollbackListCurNode;
			tmpRollbackListCurNode = nullptr;
		}

		if (rollbackListCurNode) {
			rollbackListCurNode->tmpMark = true;
		}

		npLL.onPath = false;
		npRR.onPath = false;
		npL = joinNodesInPath(npLL, npLR);
		NodeInPath npRoot = joinNodesInPath(npL, npRR);

		exposedRoot = npRoot.node;
		assert(Integrity::treeConsistency(exposedRoot));
	}

	template <class... TUserData>
	bool TopTree<TUserData...>::rollback() {
		if (Node *prevRoot = rollbackList.front()) {
			while (Node *node = rollbackList.pop()) {
				for (int i : {0, 1}) {
					assert(node->children[i]);
					if (node->children[i]->parent != node) {
						if (node->children[i]->parent) {
							// release and detach child from temporal tree
							releaseNode(node->children[i]->parent);
							node->children[i]->detach();
						}

						// attach child back to the original tree
						node->attachChild(i, node->children[i]);
					}
				}

				// destroy temporary tree if the previous one is restored
				if (node->tmpMark) {
					node->tmpMark = false;
					{
						Node *prevNode = nullptr;
						for (Node *node : exposedRoot->postorder()) {
							if (prevNode) {
								prevNode->detach();
								deleteNode(prevNode);
							}
							prevNode = node;
						}
						if (prevNode) {
							prevNode->detach();
							deleteNode(prevNode);
						}
					}
					exposedRoot = prevRoot;
					prevRoot = rollbackList.front();
					assert(Integrity::treeConsistency(exposedRoot));
				}
			}

			return true;
		}
		return false;
	}


	// --- queries ---

	template <class... TUserData>
	inline bool TopTree<TUserData...>::sameComponent(Vertex u, Vertex v) {
		if (u == v) return true;
		rollback();
		Node *uNode = vertexToNode[u];
		Node *vNode = vertexToNode[v];
		return uNode && vNode && (treeRoot(uNode) == treeRoot(vNode));
	}

	template <class... TUserData>
	inline std::pair<Vertex, Vertex> TopTree<TUserData...>::getBoundary() {
		assert(exposedRoot);
		return {exposedRoot->boundary[0], exposedRoot->boundary[1]};
	}

	template <class... TUserData>
	inline Vertex TopTree<TUserData...>::getBoundary(int index) {
		assert(exposedRoot);
		return exposedRoot->boundary[index];
	}

	template <class... TUserData>
	template <int I>
	inline auto &TopTree<TUserData...>::getExposedData() {
		assert(exposedRoot);
		validateUserData<I>(exposedRoot);
		return std::get<I>(exposedRoot->userData);
	}

	template <class... TUserData>
	template <int I>
	inline auto &TopTree<TUserData...>::getEdgeData() {
		assert(exposedRoot);
		Node *node = exposedRoot;
		while (node->clusterType != BASE) {
			assert(node->clusterType == RAKE);
			releaseJustNodeData<I>(node);
			node = node->children[0];
		}

		return std::get<I>(node->userData);
	}


	// --- internal user data manipulation ---

	template <class... TUserData>
	inline void TopTree<TUserData...>::releaseNode(Node *node) {
		assert(node->clusterType != ClusterType::BASE);
		if constexpr(sizeof...(TUserData) > 0) {
			releaseNodeData(node, std::index_sequence_for<TUserData...>());
		}
	}

	template <class... TUserData>
	inline void TopTree<TUserData...>::releaseJustNode(Node *node) {
		assert(node->clusterType != ClusterType::BASE);
		if constexpr(sizeof...(TUserData) > 0) {
			releaseJustNodeData(node, std::index_sequence_for<TUserData...>());
		}
	}

	template <class... TUserData>
	template <std::size_t... Is>
	inline void TopTree<TUserData...>::releaseJustNodeData(Node *node, std::index_sequence<Is...>) {
		(void)std::initializer_list<int>{
			(node->userDataValid[Is] ? (node->template userDataSplit<Is>(), 0) : 0)...
		};
	}

	template <class... TUserData>
	template <size_t I, size_t... Is>
	inline void TopTree<TUserData...>::releaseNodeData(Node *node, std::index_sequence<I, Is...>) {
		assert(node->clusterType != ClusterType::BASE);
		Node *origNode = node;
		InnerList<Node, &Node::tmpListNode> nodes;
		while (node && node->userDataValid[I]) {
			nodes.push(node);
			node = node->parent;
		}
		while (Node *node = nodes.pop()) {
			assert(node->userDataValid[I]);
			node->template userDataSplit<I>();
		}
		if constexpr(sizeof...(Is) > 0) {
			releaseNodeData<Is...>(origNode);
		}
	}

	template <class... TUserData>
	template <size_t I>
	inline void TopTree<TUserData...>::validateUserData(Node *root, InnerList<Node, &Node::rollbackListNode> &splitList) {
		assert(root);
		if (root->userDataValid[I]) return;
		for (Node *node = splitList.front(); node; node = splitList.next(node)) {
			if (node->userDataValid[I]) {
				node->template userDataSplit<I>();
			}
		}
		for (Node *node : root->postorder([](Node *node){return !node->userDataValid[I];})) {
			assert(node->clusterType != ClusterType::BASE);
			node->template userDataJoin<I>();
		}
	}


	// === TopTree::Node methods' implementation ===

	// --- tree structure manipulation ---

	template <class... TUserData>
	inline std::tuple<typename TopTree<TUserData...>::NodeInPath, typename TopTree<TUserData...>::NodeInPath> TopTree<TUserData...>::Node::getChildren(Vertex firstVertex, bool rev) {
		assert(clusterType != BASE);
		rev ^= (boundary[0] != firstVertex);
		NodeInPath ret[2] = {
			{ .node = children[0], .onPath = true },
			{ .node = children[1], .onPath = (clusterType == COMPRESS) }};
		return std::tuple(ret[rev], ret[!rev]);
	}


	template <class... TUserData>
	inline typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::attachChildren(ClusterType type, Node *child1, Node *child2) {
		attachChild(0, child1);
		attachChild(1, child2);
		setBoundary(type);
		return this;
	}

	template <class... TUserData>
	inline typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::attachChild(int childIndex, Node* child) {
		this->children[childIndex] = child;
		child->parent = this;
		return this;
	}

	template <class... TUserData>
	inline typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::detach() {
		Node *p = this->parent;
		if (p) {
			assert(p->children[ p->children[0] != this ] == this);
			assert(!p->userDataValid[0]); // XXX
			p->children[ p->children[0] != this ] = nullptr;
			this->parent = nullptr;
		}
		return p;
	}

	template <class... TUserData>
	inline typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::detachChildren() {
		for (size_t i : {0, 1}) {
			if (this->children[i]) {
				this->children[i]->detach();
			}
		}
		return this;
	}


	// --- user data manipulation ---

	template <class... TUserData>
	template <size_t I>
	inline EventData<typename std::tuple_element<I, std::tuple<TUserData...>>::type> TopTree<TUserData...>::Node::getEventData() {
		return EventData<typename std::tuple_element<I, std::tuple<TUserData...>>::type>(
			clusterType,
			std::get<I>(this->userData),
			std::get<I>(children[0]->userData),
			std::get<I>(children[1]->userData),
			boundary,
			this->getInnerVertex(),
			children[0]->boundary, children[1]->boundary);
	}

	template <class... TUserData>
	template <size_t I>
	inline void TopTree<TUserData...>::Node::userDataSplit() {
		assert(userDataValid[I]);
		assert(!this->parent || !this->parent->userDataValid[I]);
		std::tuple_element<I, decltype(userData)>::type::split(getEventData<I>());
		userDataValid[I] = false;
	}

	template <class... TUserData>
	template <size_t I>
	inline void TopTree<TUserData...>::Node::userDataJoin() {
		assert(!userDataValid[I]);
		assert(children[0]->userDataValid[I] && children[1]->userDataValid[I]);
		std::tuple_element<I, decltype(userData)>::type::join(getEventData<I>());
		userDataValid[I] = true;
	}


	// --- cluster type and boundary ---

	template <class... TUserData>
	inline void TopTree<TUserData...>::Node::setBoundary(ClusterType type) { // for COMPRESS or RAKE node
		this->clusterType = type;
		auto [child1, child2] = this->children;
		assert((type != BASE) && child1 && child2);

		Vertex sharedVertex = child2->boundary[0];
		if ((child1->boundary[0] != sharedVertex) && (child1->boundary[1] != sharedVertex)) {
			sharedVertex = child2->boundary[1];
			assert((child1->boundary[0] == sharedVertex) || (child1->boundary[1] == sharedVertex));
		}

		this->boundary[0] = child1->getOtherVertex(sharedVertex);
		this->boundary[1] = type == COMPRESS ?
			child2->getOtherVertex(sharedVertex) : sharedVertex;

		for (bool &valid : this->userDataValid) valid = false;

		assert(Integrity::nodeChildrenBoundary(this));
	}

	template <class... TUserData>
	inline void TopTree<TUserData...>::Node::setBoundary(Vertex v1, Vertex v2) { // for BASE node
		assert((this->children[0] == nullptr) && (this->children[1] == nullptr));
		this->clusterType = BASE;
		this->boundary[0] = v1;
		this->boundary[1] = v2;
		for (bool &valid : this->userDataValid) valid = true;
	}

}

using TopTreeClusterType = TopTreeInternals::ClusterType;
using TopTreeVertex = TopTreeInternals::Vertex;
using TopTreeInternals::TopTree;

template <class TUserData>
using TopTreeEventData = TopTreeInternals::EventData<TUserData>;

#undef assert
#endif
