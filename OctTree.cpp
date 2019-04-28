#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include "OctTree.h"

OctTreeNode::OctTreeNode(double cx, double cy, double cz, double sx, double sy, double sz) {
	centerx = cx;
	centery = cy;
	centerz = cz;
	
	sizex = sx;
	sizey = sy;
	sizez = sz;
	
	for (int i = 0; i < 8; i++) {
		children[i] = nullptr;
	}
	
	molecule = nullptr;
	isLeaf = true;
}

OctTreeNode::~OctTreeNode() {
	for (int i = 0; i < 8; i++) {
		if (children[i] != nullptr)
			delete children[i];
	}
}

void OctTreeNode::AddMolecule(Molecule *molecule) {
	assert(molecule);
	
	if (this->molecule == nullptr && isLeaf) {
		this->molecule = molecule;
		return;
	}
	
	if (this->molecule != nullptr) {
		const sf::Vector3<double> *coordinates = this->molecule->getCoordinates();
		int index = 0 | (coordinates->x > centerx) << 2 | (coordinates->y > centery) << 1 | (coordinates->z > centerz);
		double newCenterX = centerx + ((index & (1 << 2)) ? (sizex / 2) : (-sizex / 2));
		double newCenterY = centery + ((index & (1 << 1)) ? (sizey / 2) : (-sizey / 2));
		double newCenterZ = centerz + ((index & 1) ? (sizez / 2) : (-sizez / 2));
		children[index] = new OctTreeNode(newCenterX, newCenterY, newCenterZ, sizex / 2, sizey / 2, sizez / 2);
		children[index]->AddMolecule(this->molecule);
		isLeaf = false;
		
		this->molecule = nullptr;
	}
	
	const sf::Vector3<double> *coordinates = molecule->getCoordinates();
	int index = 0 | (coordinates->x > centerx) << 2 | (coordinates->y > centery) << 1 | (coordinates->z > centerz);
	
	if (children[index] == nullptr) {
		double newCenterX = centerx + ((index & (1 << 2)) ? (sizex / 2) : (- sizex / 2));
		double newCenterY = centery + ((index & (1 << 1)) ? (sizey / 2) : (- sizey / 2));
		double newCenterZ = centerz + ((index & 1) ? (sizez / 2) : (- sizez / 2));
		children[index] = new OctTreeNode(newCenterX, newCenterY, newCenterZ, sizex / 2, sizey / 2, sizez / 2);
	}
	
	children[index]->AddMolecule(molecule);
	
}

OctTreeNode *OctTreeNode::getChild(int index) const {
	return children[index];
}

Molecule *OctTreeNode::getMolecule() const {
	return molecule;
}

bool OctTreeNode::isNear(Molecule *molecule) const {
	const sf::Vector3<double> *coordinates = molecule->getCoordinates();
	double radius = molecule->getRadius();
	return (coordinates->x > centerx - sizex - 2 * radius) && (coordinates->x < centerx + sizex + 2 * radius) &&
			(coordinates->y > centery - sizey - 2 * radius) && (coordinates->y < centery + sizey + 2 * radius) &&
			(coordinates->z > centerz - sizez - 2 * radius) && (coordinates->z < centerz + sizez + 2 * radius);
}

void OctTreeNode::printNode(FILE *file) {
	
	fprintf(file, "\tn%p [shape=record, label=\"{%g,%g,%g,%g,%g,%g|{", this, centerx, centery, centerz, sizex, sizey, sizez);
	fprintf(file, "<c0%p> 0", this);
	for (int i = 1; i < 8; i++) {
		fprintf(file, "|<c%d%p> %d", i, this, i);
	}
	fprintf(file, "}}\"];\n");
//	Fprintf("\tn%p [shape=record, label=\"{%s|{<l%p> left|<r%p> right}}\"];\n", node, node->val, node, node);
//	if (node->left)
//		Fprintf("\tn%p:l%p -> n%p:n;\n", node, node, node->left);
//	if (node->right)
//		Fprintf("\tn%p:r%p -> n%p:n;\n", node, node, node->right);
//	if (node->parent)
//		Fprintf("\tn%p:n -> n%p;\n", node, node->parent);
	
	for (int i = 0; i < 8; i++) {
		if (this->children[i]) {
			fprintf(file, "\tn%p:c%d%p -> n%p:n;\n", this, i, this, this->children[i]);
			getChild(i)->printNode(file);
		}
	}
}

OctTree::OctTree(double length, double width, double height) {
	this->length = length;
	this->width = width;
	this->height = height;
	
	root = nullptr;
}

OctTree::~OctTree() {
	delete root;
}

void OctTree::clear() {
	delete root;
	root = nullptr;
}

void OctTree::AddMolecule(Molecule *molecule) {
	assert(molecule);
	
	if (root == nullptr) {
		root = new OctTreeNode(length / 2, width / 2, height / 2, length / 2, width / 2, height / 2);
	}
	
	root->AddMolecule(molecule);
}

const OctTreeNode *OctTree::gerRoot() {
	return root;
}

void OctTree::dump() {
	plotGraph();
}

void OctTree::plotGraph() {
	assert(this);
	
	FILE *dumpFile = fopen("tree.dv", "wb");
	fprintf(dumpFile, "digraph tree {\n");
	
	root->printNode(dumpFile);
	
	fprintf(dumpFile, "}\n");
	fclose(dumpFile);
	
	system("dotty tree.dv");
	
	remove("tree.dv");
}
