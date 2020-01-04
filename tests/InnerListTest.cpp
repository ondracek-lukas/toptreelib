// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#include <cassert>
#include "../TopTreeInternals/InnerList.hpp"
#include <vector>
#include <list>
#include <cstdio>
#include <utility>

struct A {
	int value = 0;
	A *innerListNode = nullptr;
};


void assertList(TopTreeInternals::InnerList<A, &A::innerListNode> &list, std::list<int> &&expected) {
	for (A *a = list.front(); a; a = list.next(a)) {
		assert(a->value == expected.front());
		expected.pop_front();
	}
	assert(expected.size() == 0);
}

int main() {
	int i;
	std::vector<A> vect;

	{ // vect init
		vect.resize(10);
		i = 0;
		for (A &a : vect) {
			a.value = i++;
		}
	}

	TopTreeInternals::InnerList<A, &A::innerListNode> list;
	// list = {}

	{ // test: push, pop
		for (A &a : vect) {
			list.push(&a);
		}

		// list = { 9, 8, ..., 0 }

		i = 10;
		while (A *a = list.pop()) {
			--i;
			assert(a->value == i);
		}
		assert(i == 0);
	}

	// list = {}

	{ // test: destructor
		TopTreeInternals::InnerList<A, &A::innerListNode> list2;
		for (A &a : vect) {
			list2.push(&a);
		}
	}
	for (A &a : vect) {
		assert(!a.innerListNode);
	}

	{ // test: append, front, next
		for (A &a : vect) {
			list.append(&a);
		}

		// list = { 0, 1, ..., 9 }

		i = 0;
		for (A *a = list.front(); a; a = list.next(a)) {
			assert(a->value == i);
			i++;
		}
		assert(i == 10);
	}

	// list = { 0, 1, ..., 9 }

	{ // test: removeAfter, insertAfter
		i = list.removeAfter(&vect[3])->value;
		assert(i == 4);
		i = list.removeAfter(&vect[8])->value;
		assert(i == 9);
		list.insertAfter(&vect[8], &vect[4]);
		list.insertAfter(&vect[3], &vect[9]);

		assertList(list, { 0, 1, 2, 3, 9, 5, 6, 7, 8, 4 });
	}

	// list = { 0, 1, 2, 3, 9, 5, 6, 7, 8, 4 }

	{ // test: push list
		TopTreeInternals::InnerList<A, &A::innerListNode> list2;
		for (i = 0; i < 5; i++) {
			list2.push(list.pop());
		}
		list.push(std::move(list2));
		assert(!list2.front() && !list2.back());
		assertList(list, { 9, 3, 2, 1, 0, 5, 6, 7, 8, 4 });
	}

	{ // test: insertAfter list
		TopTreeInternals::InnerList<A, &A::innerListNode> list2;
		for (i = 0; i < 5; i++) {
			list2.push(list.pop());
		}
		A *a = list.front();
		a = list.next(a);
		a = list.next(a);
		list.insertAfter(a, std::move(list2));
		assert(!list2.front() && !list2.back());
		assertList(list, { 5, 6, 7, 0, 1, 2, 3, 9, 8, 4 });
	}


	{ // test: append list
		TopTreeInternals::InnerList<A, &A::innerListNode> list2;
		for (i = 0; i < 6; i++) {
			list2.append(list.pop());
		}
		list.append(std::move(list2));
		assert(!list2.front() && !list2.back());
		assertList(list, { 3, 9, 8, 4, 5, 6, 7, 0, 1, 2 });
	}

	{ // test: move constructor and move assignment
		TopTreeInternals::InnerList<A, &A::innerListNode> list2 = std::move(list);
		assert(!list.front() && !list.back());
		A a, b;
		list.push(&a);
		list.push(&b);
		assert(b.innerListNode);
		list = std::move(list2);
		assert(!list2.front() && !list2.back());
		assertList(list, { 3, 9, 8, 4, 5, 6, 7, 0, 1, 2 });
		assert(!b.innerListNode);
	}



	return 0;
}
