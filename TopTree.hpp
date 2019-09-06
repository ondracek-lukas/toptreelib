// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

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
	enum ClusterType { BASE, COMPRESS, RAKE };
	typedef int Vertex; // vertex of underlying tree
}

#if TOP_TREE_INTEGRITY_LEVEL > 0
#include "tests/TopTreeIntegrity.hpp"
#else
template <class... TUserData> class TopTreeIntegrity {};
#define assert(cond)
#endif

#include <tuple>
#include <vector>
#include <list>

#include "TopTreeInternals/SubtreeTraversability.hpp"
#include "TopTreeInternals/InnerList.hpp"

namespace TopTreeInternals {

	template <class... TUserData>
	class TopTree {
		protected:
			using Integrity = TopTreeIntegrity<TUserData...>; friend Integrity;

		public:
			virtual bool expose(Vertex u, Vertex v);
			virtual bool exposeTree(Vertex v); // set tree with v as main
			virtual void link(Vertex u, Vertex v, TUserData... userData) = 0;
			virtual std::tuple<TUserData...> cut(Vertex u, Vertex v) = 0;


			template <int TUserDataIndex = 0, bool onExposedPath = false, class TSelect>
			void search(TSelect select); // not virtual; redeclare using std:function if needed

			template <int TUserDataIndex = 0, class TSelect>
			void pathSearch(TSelect select) { search<TUserDataIndex, true>(select); };


			Vertex newVertex() {
				vertexToNode.emplace_back(nullptr);
				return vertexToNode.size()-1;
			};

			int getVerticesCnt() {
				return vertexToNode.size();
			}

			bool sameComponent(Vertex u, Vertex v) { // XXX requires rollback
				if (u == v) return true;
				Node *uNode = vertexToNode[u];
				Node *vNode = vertexToNode[v];
				return uNode && vNode && (treeRoot(uNode) == treeRoot(vNode));
			}

			std::tuple<TUserData...> getRootData() {
				assert(exposedRoot);
				return exposedRoot->userData;
			}
			//void setRootData(const TUserData... &userData); // (or maybe set one userData at a time)

			std::pair<Vertex, Vertex> getBoundary() {
				assert(exposedRoot);
				return {exposedRoot->boundary[0], exposedRoot->boundary[1]};
			}

			// User can request which data of TUserData... should be updated and which not.
			// The TopTree changes the state as soon as it is consistent
			// and then calls split/join only on active data.
			//void setActiveData(bool.../*somehow make number of parameters the same as of TUserData*/ activeUserData);

			// ...

		protected:  // to be accessible by driver

			struct Node;

			struct NodeInPath {
				Node *node = nullptr;
				bool onPath = true;
			};

			struct Node : SubtreeTraversability<Node> {  // node in top tree
				size_t index; // < nodesAllocated, to be used in vectors with extension data
				Vertex boundary[2];
				Node *children[2] = {nullptr, nullptr};
				Node *parent = nullptr;
				// constraints:
					// boundary[i] should be in children[i]->boundary;
				ClusterType clusterType;
				std::tuple<TUserData...> userData;
				bool userDataChanged[sizeof...(TUserData)];
				bool userDataValid = false;
				Node *tmpListNode = nullptr;
				Node *rollbackListNode = nullptr;
				bool tmpMark = false;
				// Node *temporalListOfNodes; // several such lists for different purposes
				// ...

				Vertex getOtherVertex(Vertex vertex) {
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

				template <std::size_t... Is>
				void userDataSplitSeq(std::index_sequence<Is...>);

				void userDataSplit() {
					userDataSplitSeq(std::index_sequence_for<TUserData...>{});
				}

				template <std::size_t... Is>
				void userDataJoinSeq(std::index_sequence<Is...>);

				void userDataJoin() {
					userDataJoinSeq(std::index_sequence_for<TUserData...>{});
				}

				// For base cluster use setBoundary;
				// for non-base cluster use either attachChildren, or twice attachChild and setBoundary.

				Node *attachChildren(ClusterType type, Node *child1, Node *child2);
				Node *attachChild(int childIndex, Node* child);
				Node *detachChildren();
				Node *detach();

				void setBoundary(ClusterType type); // for COMPRESS or RAKE node
				void setBoundary(Vertex v1, Vertex v2); // for BASE node

			};

			Node *treeRoot(Node *node) {
				assert(node);
				while (node->parent) node = node->parent;
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

			std::vector<Node*> vertexToNode; // topmost nodes having given vertex as inner vertex, or root

			bool activeUserData[sizeof...(TUserData)];    // currently active user data
			bool activeUserDataReq[sizeof...(TUserData)]; // requested state


			// To correctly handle user data when tree is being changed,
			// driver should always call releaseNode before changing its children;
			// then it can perform arbitrary changes (in any order);
			// finally it should call validateTree on new tree root(s).

			// Calls split on node and all ancestors with valid user data top-down
			// marking them invalid.
			void releaseNode(Node *node) {
				assert(node->clusterType != ClusterType::BASE);
				InnerList<Node, &Node::tmpListNode> nodes;
				while (node && node->userDataValid) {
					nodes.push(node);
					node = node->parent;
				}
				while (Node *node = nodes.pop()) {
					node->userDataSplit();
					node->userDataValid = false;
				}
			}

			void releaseJustNode(Node *node) {
				assert(!node->parent || !node->parent->userDataValid);
				if (node->userDataValid) {
					node->userDataSplit();
					node->userDataValid = false;
				}
			}

			// Calls join on all descendants with invalid user data bottom-up
			// marking them valid.
			void validateTree(Node *root) {
				for (Node *node : root->postorder([](Node *node){return !node->userDataValid;})) {
					assert(node->clusterType != ClusterType::BASE);
					node->userDataJoin();
					node->userDataValid = true;
				}
			}

			// validateTree can be called either at the end of link/cut/(expose) operation,
			// or at the beginning of the following expose/search/... operation
			// allowing composite updates.


			// -- tree structure manipulation --

			// Methods newNode and deleteNode can be redefined (masked) in driver
			// to add further functionality,
			// the original methods should be called from new ones.
			// All nodes created by the base class are deleted in rollabck method.

			size_t nodesAllocated = 0;
			InnerList<Node, &Node::tmpListNode> freeNodes;

			Node *newNode() {
				Node *node = freeNodes.pop();
				if (!node) {
					node = new Node;
					node->index = nodesAllocated++;
				}
				return node;
			}

			void deleteNode(Node *node) { // TODO: add deallocation in TopTree destructor
				assert(!node->parent && !node->children[0] && !node->children[1]);
				freeNodes.push(node);
			}

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

			void markNodeAsNonRoot(Node *node) {
				vertexToNode[node->boundary[0]] = nullptr;
				vertexToNode[node->boundary[1]] = nullptr;
			}

		private:
			NodeInPath joinNodesInPath(NodeInPath np1, NodeInPath np2);
	};


	// --- TopTree methods' implementation ---

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
			releaseJustNode(npM.node);
			rollbackList.push(npM.node);

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
			if (npM.node) {
				npM.node->tmpMark = false;
			}
		}
		exposedRoot = npR.node;

		assert(Integrity::treeConsistency(exposedRoot));
		validateTree(exposedRoot);
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
	template <int TUserDataIndex, bool onExposedPath, class TSelect>
	inline void TopTree<TUserData...>::search(TSelect select) {
		if (!exposedRoot || (exposedRoot->clusterType == BASE)) {
			return;
		}
		if (!onExposedPath) {
			rollback();
		}

		// We maintain path: npL = (npLL, npRR), vM, npR = (npRL, npRR)
		// Inv.: npL either equals npLR, or is temporary; symmetrically for npR, npRL
		NodeInPath npL, npLL, npLR, npR, npRL, npRR;
		Vertex vM;

		exposedRoot->tmpMark = true;
		npL.node = npLR.node = exposedRoot;
		vM = npLR.node->boundary[0];

		while (npLR.node->clusterType != BASE) {

			// npL, npR are not set at this point; npRL is empty; npLR should be split

			releaseJustNode(npLR.node);
			rollbackList.push(npLR.node);

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
				validateTree(npRoot.node); // XXX not all user data need updating here

				reverse = select(
					npRoot.node->clusterType,
					std::get<TUserDataIndex>(npRoot.node->userData),
					std::get<TUserDataIndex>(npRoot.node->children[0]->userData),
					std::get<TUserDataIndex>(npRoot.node->children[1]->userData),
					npRoot.node->boundary[0],
					npRoot.node->boundary[1],
					npRoot.node->getInnerVertex());
				reverse ^= npRoot.node->children[0] != npL.node;

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
				npLL.onPath = false;
				npRR = joinNodesInPath(npRR, npLL);
				npLL = {};
				npLR.onPath = true;
			}
		}

		npLL.onPath = false;
		npRR.onPath = false;
		npL = joinNodesInPath(npLL, npLR);
		NodeInPath npRoot = joinNodesInPath(npL, npRR);

		exposedRoot = npRoot.node;
		assert(Integrity::treeConsistency(exposedRoot));
		validateTree(exposedRoot);
	}

	template <class... TUserData>
	bool TopTree<TUserData...>::rollback() {
		if (rollbackList.front()) {
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

				// destroy temporary tree if previous root was found
				if (node->tmpMark) {
					node->tmpMark = false;
					for (Node *node : exposedRoot->postorder()) {
						node->detach();
						deleteNode(node);
					}
					exposedRoot = node;
					assert(Integrity::treeConsistency(exposedRoot));
				}
			}

			validateTree(exposedRoot);
			return true;
		}
		return false;
	}

	// --- TopTree::Node methods' implementation ---


	template <class... TUserData>
	std::tuple<typename TopTree<TUserData...>::NodeInPath, typename TopTree<TUserData...>::NodeInPath> TopTree<TUserData...>::Node::getChildren(Vertex firstVertex, bool rev) {
		assert(clusterType != BASE);
		rev ^= (boundary[0] != firstVertex);
		NodeInPath ret[2] = {
			{ .node = children[0], .onPath = true },
			{ .node = children[1], .onPath = (clusterType == COMPRESS) }};
		return std::tuple(ret[rev], ret[!rev]);
	}

	// TODO: split and join only if enabled for user data

	template <class... TUserData>
	template <std::size_t... Is>
	void TopTree<TUserData...>::Node::userDataSplitSeq(std::index_sequence<Is...>) {
		Vertex innerVertex = this->getInnerVertex();
		(void)std::initializer_list<int>{
			(std::get<Is>(userData).split(
				clusterType,
				std::get<Is>(this->userData),
				std::get<Is>(children[0]->userData),
				std::get<Is>(children[1]->userData),
				boundary[0],
				boundary[1],
				innerVertex
			), 0)...
		};
	}


	template <class... TUserData>
	template <std::size_t... Is>
	void TopTree<TUserData...>::Node::userDataJoinSeq(std::index_sequence<Is...>) {
		Vertex innerVertex = this->getInnerVertex();
		(void)std::initializer_list<int>{
			(std::get<Is>(userData).join(
				clusterType,
				std::get<Is>(this->userData),
				std::get<Is>(children[0]->userData),
				std::get<Is>(children[1]->userData),
				boundary[0],
				boundary[1],
				innerVertex
			), 0)...
		};
	}


	template <class... TUserData>
	typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::attachChildren(ClusterType type, Node *child1, Node *child2) {
		attachChild(0, child1);
		attachChild(1, child2);
		setBoundary(type);
		return this;
	}

	template <class... TUserData>
	typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::attachChild(int childIndex, Node* child) {
		this->children[childIndex] = child;
		child->parent = this;
		return this;
	}

	template <class... TUserData>
	typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::detach() {
		Node *p = this->parent;
		if (p) {
			assert(p->children[ p->children[0] != this ] == this);
			p->children[ p->children[0] != this ] = nullptr;
			this->parent = nullptr;
		}
		return p;
	}

	template <class... TUserData>
	typename TopTree<TUserData...>::Node *TopTree<TUserData...>::Node::detachChildren() {
		for (size_t i : {0, 1}) {
			if (this->children[i]) {
				this->children[i]->detach();
			}
		}
		return this;
	}

	template <class... TUserData>
	void TopTree<TUserData...>::Node::setBoundary(ClusterType type) { // for COMPRESS or RAKE node
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

		assert(Integrity::nodeChildrenBoundary(this));
	}

	template <class... TUserData>
	void TopTree<TUserData...>::Node::setBoundary(Vertex v1, Vertex v2) { // for BASE node
		assert((this->children[0] == nullptr) && (this->children[1] == nullptr));
		this->clusterType = BASE;
		this->boundary[0] = v1;
		this->boundary[1] = v2;
		this->userDataValid = true;
	}

}

using TopTreeClusterType = TopTreeInternals::ClusterType;
using TopTreeVertex = TopTreeInternals::Vertex;
using TopTreeInternals::TopTree;

#undef assert
#endif
