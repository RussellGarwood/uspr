#ifndef INCLUDE_USPR_NEIGHBORS
#define INCLUDE_USPR_NEIGHBORS

// INCLUDES
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>
#include <time.h>

#include "utree.h"

using namespace std;

// FUNCTIONS

list<utree *> get_neighbors(utree *T, set<string> *known_trees = NULL);
void get_neighbors(utree *T, unode *prev, unode *current, list<utree *> &neighbors, set<string> *known_trees = NULL);
void get_neighbors(utree *T, unode *x, unode *y, unode *prev, unode *current, list<utree *> &neighbors, set<string> *known_trees = NULL);
void add_neighbor(utree *T, unode *x, unode *y, unode *w, unode *z, list<utree *> &neighbors, set<string> *known_trees = NULL);


list<utree *> get_neighbors(utree *T, set<string> *known_trees) {
	list<utree *> neighbors = list<utree *>();
	unode *root = T->get_node(T->get_smallest_leaf());
	get_neighbors(T, NULL, root, neighbors, known_trees);
	return neighbors;
}

// enumerate the source edges
void get_neighbors(utree *T, unode *prev, unode *current, list<utree *> &neighbors, set<string> *known_trees) {
	// continue enumerating choices of the first edge
	list<unode *> c_neighbors = current->get_neighbors();
	for (unode *next : c_neighbors) {
		if (next != prev) {
			get_neighbors(T, current, next, neighbors, known_trees);
		}
	}
	// try moving both sides of the edge (prev, current)
	if (prev != NULL) {
		get_neighbors(T, prev, current, prev, current, neighbors, known_trees);
		get_neighbors(T, current, prev, current, prev, neighbors, known_trees);
	}
}

// enumerate the target edges
void get_neighbors(utree *T, unode *x, unode *y, unode *prev, unode *current, list<utree *> &neighbors, set<string> *known_trees) {
	// continue enumerating choices of the second edge
	// copy the neighbor list as it may change
	list<unode *> c_neighbors = current->get_neighbors();
	for (unode *next : c_neighbors) {
		if (next != prev) {
			get_neighbors(T, x, y, current, next, neighbors, known_trees);
		}
	}
	// test the spr move (T, x, y, prev, current)
	if (prev != NULL) {
		add_neighbor(T, x, y, prev, current, neighbors, known_trees);
	}
}

void add_neighbor(utree *T, unode *x, unode *y, unode *w, unode *z, list<utree *> &neighbors, set<string> *known_trees) {
	// check for duplicate SPR moves
	if (x == y ||
			y == w ||
			y == z ) {
		return;
	}
	// TODO: other duplicates? probably NNIs, like with rooted SPR?

	// node info so the uspr can be reversed
	unode *yprime = NULL;
	unode *y1 = NULL;
	unode *y2 = NULL;
	// apply the spr
	cout << endl;
	cout << "T: " << T->str(true) << endl;
	cout << "\tx: " << x->get_label() << endl;
	cout << "\ty: " << y->get_label() << endl;
	cout << "\tw: " << w->get_label() << endl;
	cout << "\tz: " << z->get_label() << endl;
	T->uspr(x, y, w, z, &yprime, &y1, &y2);
	// normalize the tree
	T->normalize_order();
	// print the tree
	cout << "neighbor: " << T->str() << endl;
	// revert the SPR
	T->uspr(x, yprime, y1, y2);
	T->normalize_order();
	cout << "T: " << T->str() << endl;
	cout << endl;
	return;
}
#endif
