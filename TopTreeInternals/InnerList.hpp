// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* Singly-linked list of TNode objects
 * with next pointers stored in TNode's member variable *MNext.
 */

#ifndef TOP_TREE_INNER_LIST
#define TOP_TREE_INNER_LIST

#include <utility>

#ifndef assert
#define assert(cond)
#endif

namespace TopTreeInternals {

	template <typename TNode, TNode *TNode::*MNext>
	class InnerList {
		private:
			TNode *first = nullptr;
			TNode *last = nullptr;

			void insert(TNode **listPtr, TNode *node) {
				assert(node && !(node->*MNext));
				if (!*listPtr) last = node;
				node->*MNext = *listPtr;
				*listPtr = node;
			}
			void insert(TNode **listPtr, InnerList<TNode, MNext> &&list) {
				if (list.last) {
					if (!*listPtr) last = list.last;
					list.last->*MNext = *listPtr;
					*listPtr = list.first;
					list.first = nullptr;
					list.last  = nullptr;
				}
			}

			TNode *remove(TNode **listPtr, TNode *prev) {
				TNode *node = *listPtr;
				if (node) {
					*listPtr = node->*MNext;
					node->*MNext = nullptr;
					if (!*listPtr) last = prev;
				}
				return node;
			}

		public:
			InnerList() = default;
			InnerList(InnerList &&old) : first(old.first), last(old.last) {
				old.first = nullptr;
				old.last  = nullptr;
			}
			InnerList &operator=(InnerList &&old) {
				while (pop());
				first = old.first;
				last  = old.last;
				old.first = nullptr;
				old.last = nullptr;
				return *this;
			}

			TNode *front() {
				return first;
			}
			TNode *back() {
				return last;
			}

			TNode *next(TNode *node) {
				assert(node);
				return node->*MNext;
			}

			void insertAfter(TNode *after, TNode *node) {
				insert(after ? &(after->*MNext) : &first, node);
			}
			void insertAfter(TNode *after, InnerList<TNode, MNext> &&list) {
				insert(after ? &(after->*MNext) : &first, std::move(list));
			}

			TNode *removeAfter(TNode *after) {
				assert(after);
				return remove(&(after->*MNext), after);
			}

			void push(TNode *node) {
				insert(&first, node);
			}
			void push(InnerList<TNode, MNext> &&list) {
				insert(&first, std::move(list));
			}

			TNode *pop() {
				return remove(&first, nullptr);
			}

			void append(TNode *node) {
				insert(last ? &(last->*MNext) : &first, node);
			}
			void append(InnerList<TNode, MNext> &&list) {
				insert(last ? &(last->*MNext) : &first, std::move(list));
			}

			~InnerList() {
				while (pop());
			}
	};
}

#endif
