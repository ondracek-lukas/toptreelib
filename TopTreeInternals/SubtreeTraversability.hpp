// TopTreeLibrary  Copyright (C) 2022  Lukáš Ondráček <ondracek.lukas@gmail.com>, use under MIT license

/* The SubtreeTraversability class
 * creates iterators for (conditional) preorder and postorder subtree traversal.
 *
 * It should be derived by a tree node class (TNode) containing
 *   TNode *parent,
 *   TNode *children[2].
 *
 * It contains the following methods returning an iterable object:
 *   preorder(cond = true),
 *   postorder(cond = true),
 *   preorderFB(cond),
 *   postorderFB(cond).
 * The cond parameter should be a lambda function deciding
 * whether to visit the current node and its subtree or not:
 *   bool cond(TNode *node).
 * FB (false border) versions visit also nodes not satisfying cond
 * but still not subtrees of their children.
 */

#ifndef TOP_TREE_SUBTREE_TRAVERSABILITY_HPP
#define TOP_TREE_SUBTREE_TRAVERSABILITY_HPP

namespace TopTreeInternals {

	template <typename TNode>
	class SubtreeTraversability {

		private:
			template <bool postorder, bool visitFalseBorder, typename TCond>
			class SubtreeTraversal {
				private:
					TNode *root;
					TCond cond;
				public:
					SubtreeTraversal(TNode *root, TCond cond) : root(root), cond(cond) {};
					class Iterator {
						friend SubtreeTraversal;
						private:
							TNode *root;
							TNode *node;
							enum Dir { LEFT = 0, RIGHT, UP } dir; // direction where to go next
							TCond cond;
							explicit Iterator(TNode *root, TNode *node, TCond cond) : root(root), node(node), cond(cond) {
								if (node) {
									bool condRes = (!visitFalseBorder || node->children[0] || node->children[1]) && cond(node);
									if ((node->children[0] || node->children[1]) && condRes) {
										dir = node->children[0] ? LEFT : RIGHT;
										if (postorder) ++*this;
									} else {
										dir = UP;
										if (!visitFalseBorder && !condRes) {
											this->node = nullptr;
										}
									}
								}
							}
						public:
							Iterator& operator++()
							{
								while (node) {
									bool condRes;
									if (dir != UP) {
										node = node->children[dir];
										condRes = (!visitFalseBorder || node->children[0] || node->children[1]) && cond(node);
										if ((node->children[0] || node->children[1]) && condRes) {
											dir = node->children[0] ? LEFT : RIGHT;
										} else {
											dir = UP;
										}
										if (!postorder && (visitFalseBorder || condRes)) break;
									} else {
										dir =
											!node->parent ||
											((node->parent->children[0] == node) && (node->parent->children[1]))
											? RIGHT : UP;
										node = node->parent;
										condRes = true;
									}
									if (postorder && (dir == UP) && (visitFalseBorder || condRes)) break;
								}
								return *this;
							}
							bool operator==(const Iterator& other) const { return node == other.node; }
							bool operator!=(const Iterator& other) const { return !(*this == other); }
							TNode *operator*() const { return node; }
					};
					Iterator begin() {return Iterator(root, root, cond); };
					Iterator end() { return Iterator(root, nullptr, cond); };
			};

			static constexpr auto SubtreeTraversalNoCond = [](TNode*){return true;};

		public:
			template <typename TCond = decltype(SubtreeTraversalNoCond)>
			SubtreeTraversal<false, false, TCond> preorder(TCond cond = SubtreeTraversalNoCond) {
				return SubtreeTraversal<false, false, TCond> (static_cast<TNode*>(this), cond);
			}

			template <typename TCond>
			SubtreeTraversal<false, true, TCond> preorderFB(TCond cond) {
				return SubtreeTraversal<false, true, TCond> (static_cast<TNode*>(this), cond);
			}

			template <typename TCond = decltype(SubtreeTraversalNoCond)>
			SubtreeTraversal<true, false, TCond> postorder(TCond cond = SubtreeTraversalNoCond) {
				return SubtreeTraversal<true, false, TCond> (static_cast<TNode*>(this), cond);
			}

			template <typename TCond>
			SubtreeTraversal<true, true, TCond> postorderFB(TCond cond) {
				return SubtreeTraversal<true, true, TCond> (static_cast<TNode*>(this), cond);
			}

	};
}

#endif
