#pragma once

#include "bbox.h"
#include <stack>

// maximum leaf size
static const uint LEAF_SIZE = 8;

template<class GeomIdx>
struct GeomBlob
{
	BBox3d  bbox;
	Vec3d   point; // representative point, usually the box midpoint
	GeomIdx id;
};

template<class GeomIdx>
struct AABVHNode
{
	BBox3d                          bbox;
	AABVHNode                       *left;
	AABVHNode                       *right;
	ShortVec<uint, LEAF_SIZE>       blobids;
	inline bool isLeaf() const { return left == NULL; }
};

template<class GeomIdx>
class AABVH
{
public:
	AABVH(const std::vector< GeomBlob<GeomIdx> > &geoms) :	root(NULL), blobs(geoms), tmpids(geoms.size())
	{
		ENSURE(blobs.size() > 0);

		for (uint k = 0; k < tmpids.size(); k++)
			tmpids[k] = k;

		root = constructTree(0, tmpids.size(), 2);
	}
	~AABVH() {}

private:
	// process range of tmpids including begin, excluding end last_dim provides a hint by saying which dimension a
	// split was last made along
	AABVHNode<GeomIdx>* constructTree(uint begin, uint end, uint last_dim)
	{
		ENSURE(end - begin > 0); // don't tell me to build a tree from nothing
								 // base case
		if (end - begin <= LEAF_SIZE) {
			AABVHNode<GeomIdx> *node = node_pool.alloc();
			node->left = NULL;
			node->blobids.resize(end - begin);
			for (uint k = 0; k < end - begin; k++) {
				uint blobid = node->blobids[k] = tmpids[begin + k];
				node->bbox = convex(node->bbox, blobs[blobid].bbox);
			}
			return node;
		}
		// otherwise, let's try to split this geometry up

		uint dim = (last_dim + 1) % 3;
		uint mid = (begin + end) / 2;
		quickSelect(mid, begin, end, dim);

		// now recurse
		AABVHNode<GeomIdx> *node = node_pool.alloc();
		node->left = constructTree(begin, mid, dim);
		node->right = constructTree(mid, end, dim);
		node->bbox = convex(node->left->bbox, node->right->bbox);
		return node;
	}

	// precondition: begin <= select < end
	void quickSelect(uint select, uint begin, uint end, uint dim)
	{
		// NOTE: values equal to the pivot may appear on either side of the split
		if (end - 1 == select)     return;

		// p(ivot)i(ndex) and p(ivot)v(alue)
		uint pi = randMod(end - begin) + begin;
		double pv = blobs[tmpids[pi]].point[dim];

		int front = begin;
		int back = end - 1;
		while (front < back) {
			if (blobs[tmpids[front]].point[dim] < pv) {
				front++;
			}
			else if (blobs[tmpids[back]].point[dim] > pv) {
				back--;
			}
			else {
				std::swap(tmpids[front], tmpids[back]);
				front++;
				back--;
			}
		}
		if (front == back && blobs[tmpids[front]].point[dim] <= pv) {
			front++;
		}

		if (select < uint(front)) {
			quickSelect(select, begin, front, dim);
		}
		else {
			quickSelect(select, front, end, dim);
		}
	}

public:
	AABVHNode<GeomIdx>                  *root;
	std::vector< GeomBlob<GeomIdx> >    blobs;

private:
	IterPool< AABVHNode<GeomIdx> >      node_pool;
	std::vector<uint>                   tmpids; // used during construction
};