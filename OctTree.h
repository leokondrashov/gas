#ifndef GAS_OCTTREE_H
#define GAS_OCTTREE_H

#include <SFML/System/Vector3.hpp>
#include <cstdio>
#include "Molecule.h"

class OctTreeNode {
private:
	OctTreeNode *children[8];
	double centerx, centery, centerz;
	double sizex, sizey, sizez;
	Molecule *molecule;
	bool isLeaf;
	
public:
	OctTreeNode(double cx, double cy, double cz, double sx, double sy, double sz);
	~OctTreeNode();
	void AddMolecule(Molecule *molecule);
	OctTreeNode *getChild(int index) const;
	Molecule *getMolecule() const;
	
	bool isNear(Molecule *molecule) const;
	void printNode(FILE *file);
};

class OctTree {
private:
	OctTreeNode *root;
	double length, width, height;
	
public:
	OctTree(double length, double width, double height);
	~OctTree();
	void AddMolecule(Molecule *molecule);
	void clear();
	
	const OctTreeNode *gerRoot();
	void dump();
	void plotGraph();
};

#endif //GAS_OCTTREE_H
