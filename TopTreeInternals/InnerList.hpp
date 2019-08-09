// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

/* Singly-linked list of TNode objects
 * with next pointers stored in TNode's member variable *MNext.
 */

#ifndef TOP_TREE_INNER_LIST
#define TOP_TREE_INNER_LIST

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
			TNode *front() {
				return first;
			}

			TNode *next(TNode *node) {
				assert(node);
				return node->*MNext;
			}

			void insertAfter(TNode *after, TNode *node) {
				assert(after);
				insert(&(after->*MNext), node);
			}

			TNode *removeAfter(TNode *after) {
				assert(after);
				return remove(&(after->*MNext), after);
			}

			void push(TNode *node) {
				insert(&first, node);
			}

			TNode *pop() {
				return remove(&first, nullptr);
			}

			void append(TNode *node) {
				insert(last ? &(last->*MNext) : &first, node);
			}

			~InnerList() {
				while (pop());
			}
	};
}

#endif
