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

namespace TopTreeInternals {

	template <class... TUserData>
	class TopTree {
		using Integrity = TopTreeIntegrity<TUserData...>; friend Integrity;

		public:
			virtual void expose(Vertex u, Vertex v) {rollback(); /* implicit worst-case impl. */ };
			virtual void exposeTree(Vertex v) {rollback(); /* implicit worst-case impl. */ }; // set tree with v as main
			virtual void link(Vertex u, Vertex v, TUserData... userData) = 0;
			virtual std::tuple<TUserData...> cut(Vertex u, Vertex v) = 0;

			virtual void searchBegin() {rollback(); /* implicit worst-case impl. */ };
			bool searchStep(bool chooseSecond) { /* impl. */ };

			Vertex newVertex() {
				vertexToNode.emplace_back(nullptr);
				return vertexToNode.size()-1;
			};

			int getVerticesCnt() {
				return vertexToNode.size();
			}

			bool sameComponent(Vertex u, Vertex v) {
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
			struct Node : SubtreeTraversability<Node> {  // node in top tree
				size_t index; // < nodesAllocated, to be used in vectors with extension data
				Vertex boundary[2];
				Node *children[2] = {nullptr, nullptr};
				Node *parent = nullptr;
				// constraints:
					// boundary[0] should be in children[0]->boundary;
					// if rake node, boundary should equal children[0]->boundary (incl. order);
				ClusterType clusterType;
				std::tuple<TUserData...> userData;
				bool userDataChanged[sizeof...(TUserData)];
				bool userDataValid = false;
				// Node *temporalListOfNodes; // several such lists for different purposes
				// ...

				Vertex getInnerVertex() { // or maybe store directly?
					assert(Integrity::nodeChildrenBoundary(this));
					return children[1]->boundary[ children[1]->boundary[0] == boundary[1] ];
				}

				// TODO: split and join only if enabled for user data

				template <std::size_t... Is>
				void userDataSplitSeq(std::index_sequence<Is...>) {
					Vertex innerVertex = this->getInnerVertex();
					(void)std::initializer_list<int>{
						(std::get<Is>(userData).split(
							clusterType,
							std::get<Is>(parent->userData),
							std::get<Is>(children[0]->userData),
							std::get<Is>(children[1]->userData),
							boundary[0],
							boundary[1],
							innerVertex
						), 0)...
					};
				}

				void userDataSplit() {
					userDataSplitSeq(std::index_sequence_for<TUserData...>{});
				}

				template <std::size_t... Is>
				void userDataJoinSeq(std::index_sequence<Is...>) {
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

				void userDataJoin() {
					userDataJoinSeq(std::index_sequence_for<TUserData...>{});
				}

			};

			Node *treeRoot(Node *node) {
				assert(node);
				while (node->parent) node = node->parent;
				return node;
			}

			// To be called after implicit expose or search;
			// can be extended by same actions for amortization in non-worst-case version.
			// Returns whether anything has been unrolled.
			virtual bool rollback() { /* implicit worst-case impl. */ };

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
				std::list<Node*> nodes; // TODO: rewrite without std::list
				while (node && node->userDataValid) {
					nodes.emplace_front(node);
					node = node->parent;
				}
				for (auto const& node : nodes) {
					node->userDataSplit();
					node->userDataValid = false;
				}
			}

			// Calls join on all descendants with invalid user data bottom-up
			// marking them valid.
			void validateTree(Node *root) {
				assert(Integrity::treeConsistency(root));

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

			// TODO: reuse old nodes, not to increase nodesAllocated much
			size_t nodesAllocated = 0;

			Node *newNode() {
				Node *node = new Node;
				node->index = nodesAllocated++;
				return node;
			}

			void deleteNode(Node *node) { // TODO: add deallocation in TopTree destructor
				delete node;
			}

			Node *detachSubtree(Node *subtree) {
				Node *p = subtree->parent;
				if (p) {
					if (p.children[0] && p.children[1]) {
						vertexToNode[p.getInnerVertex()] = nullptr;
					}
					releaseNode(p);
					p->children[ p->children[0] != subtree ] = nullptr;
					subtree->parent = nullptr;
				}
				return p;
			}

			Node *attachSubtree(Node *parent, int childIndex, Node* subtree) {
				parent->children[childIndex] = subtree;
				subtree->parent = parent;
				return parent;
			}

			void setNodeBoundary(Node *node, ClusterType type, Vertex v1, Vertex v2) {
				node->clusterType = type;
				node->boundary[0] = v1;
				node->boundary[1] = v2;
				if (type != BASE) {
					vertexToNode[node->getInnerVertex()] = node;
				}
				if (type == BASE) {
					node->userDataValid = true;
				}
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

	};

}

using TopTreeClusterType = TopTreeInternals::ClusterType;
using TopTreeVertex = TopTreeInternals::Vertex;
using TopTreeInternals::TopTree;

#undef assert
#endif
