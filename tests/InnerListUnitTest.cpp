// TopTreeLibrary  Copyright (C) 2019  Lukáš Ondráček <ondracek@ktiml.mff.cuni.cz>, use under GNU GPLv3

#include <cassert>
#include "../TopTreeInternals/InnerList.hpp"
#include <vector>
#include <list>
#include <cstdio>

struct A {
	int value;
	A *innerListNode;
};


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

		std::list<int> exp = { 0, 1, 2, 3, 9, 5, 6, 7, 8, 4 };
		for (A *a = list.front(); a; a = list.next(a)) {
			assert(a->value == exp.front());
			exp.pop_front();
		}
		assert(exp.size() == 0);
	}

	// list = { 0, 1, 2, 3, 9, 5, 6, 7, 8, 4 }



	return 0;
}
