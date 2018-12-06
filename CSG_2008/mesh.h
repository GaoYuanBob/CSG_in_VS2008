#pragma once

// SIMPLE USAGE:
//  In order to get standard template inclusion behavior/usage just include this file.
//  Then the entire template code will be included in the compilation unit.
//  ADVANCED USAGE:
//  Only include "mesh.decl.h" where-ever you would normally include a header file.
//  This will avoid including the implementation code in the current compilation unit.
//  Then, create a seperate cpp file which includes "mesh.h" and explicitly instantiates the template with the desired template parameters.
//  By following this scheme, you can prevent re-compiling the entire template implementation in every usage compilation unit and every
//  time those compilation units are recompiled during development.

// Declaration
#include "mesh.decl.h"

// Implementation
//////////////////////////////////////////////////////////////////////////////////////////////////////// mesh.tpp

#include <algorithm>
#include "unsafeRayTriIsct.h"
#include <cfloat>
#include <cmath>
#include <sstream>

#include "unionFind.h"

#include <map>
#include <set>

#include "memPool.h"
//////////////////////////////////////////////////////////////////////////////////////////////////// mesh.topoCache.tpp

#include "iterPool.h"

#include <vector>

/*
*  Allows for topological algorithms to manipulate a more familiar pointer data structure based on a simplicial complex.
*  This structure can be regenerated from the more basic vertex/triangle arrays using / createTopoCache() /
*  Once manipulations have been satisfactorily performed, the underlying vertex/triangle arrays can be cleaned up for
*  further use by topologically insensitive algorithms by /  commitTopoCache()  /
*/

#define INVALID_ID uint(-1)

struct TopoVert {
	uint                    ref;        // index to actual data
	void*                   data;       // algorithm specific handle

	ShortVec<Tptr, 8>       tris;       // triangles this vertex is incident on
	ShortVec<Eptr, 8>       edges;      // edges this vertex is incident on
};

struct TopoEdge {
	void*                   data;       // algorithm specific handle

	Vptr                    verts[2];   // endpoint vertices
	ShortVec<Tptr, 2>       tris;       // incident triangles
};

struct TopoTri {
	uint                    ref;        // index to actual data
	void*                   data;       // algorithm specific handle

	Vptr                    verts[3];   // vertices of this triangle
	Eptr                    edges[3];   // edges of this triangle   opposite to the given vertex
};

template<class VertData, class TriData>
struct Mesh<VertData, TriData>::TopoCache {
	IterPool<TopoVert>    verts;
	IterPool<TopoEdge>    edges;
	IterPool<TopoTri>     tris;

	Mesh *mesh;
	TopoCache(Mesh *owner);
	virtual ~TopoCache() {}

	// until commit() is called, the Mesh::verts and Mesh::tris arrays will still contain garbage entries
	void commit();
	bool isValid();
	void print();

	// helpers to create bits and pieces
	inline Vptr newVert();
	inline Eptr newEdge();
	inline Tptr newTri();

	// helpers to release bits and pieces
	inline void freeVert(Vptr);
	inline void freeEdge(Eptr);
	inline void freeTri(Tptr);

	// helper to delete geometry in a structured way
	inline void deleteTri(Tptr);

	// helper to flip triangle orientation
	inline void flipTri(Tptr);

private:
	void init();
};


template<class VertData, class TriData> inline
Vptr Mesh<VertData, TriData>::TopoCache::newVert()
{
	uint ref = mesh->verts.size();
	mesh->verts.push_back(VertData());
	Vptr v = verts.alloc(); // cache.verts
	v->ref = ref;
	return v;
}
template<class VertData, class TriData> inline
Eptr Mesh<VertData, TriData>::TopoCache::newEdge()
{
	Eptr e = edges.alloc(); // cache.edges
	return e;
}
template<class VertData, class TriData> inline
Tptr Mesh<VertData, TriData>::TopoCache::newTri()
{
	uint ref = mesh->tris.size();
	mesh->tris.push_back(Tri());
	Tptr t = tris.alloc(); // cache.tris
	t->ref = ref;
	return t;
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeVert(Vptr v)
{
	verts.free(v);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeEdge(Eptr e)
{
	edges.free(e);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeTri(Tptr t)
{
	tris.free(t);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::deleteTri(Tptr tri)
{
	// first, unhook the triangle from its faces
	for (uint k = 0; k<3; k++)
	{
		Vptr v = tri->verts[k];
		v->tris.erase(tri);
		Eptr e = tri->edges[k];
		e->tris.erase(tri);
	}
	// now, let's check for any edges which no longer border triangles
	for (uint k = 0; k<3; k++) {
		Eptr e = tri->edges[k];
		if (e->tris.size() == 0)
		{
			// delete edge unhook from vertices
			Vptr v0 = e->verts[0];
			v0->edges.erase(e);
			Vptr v1 = e->verts[1];
			v1->edges.erase(e);
			freeEdge(e);
		}
	}
	// now, let's check for any vertices which no longer border triangles
	for (uint k = 0; k<3; k++) {
		Vptr v = tri->verts[k];
		if (v->tris.size() == 0) {
			freeVert(v);
		}
	}
	// finally, release the triangle
	freeTri(tri);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::flipTri(Tptr t)
{
	std::swap(t->verts[0], t->verts[1]);
	std::swap(t->edges[0], t->edges[1]);
	std::swap(mesh->tris[t->ref].v[0], mesh->tris[t->ref].v[1]);
}

template<class VertData, class TriData>
Mesh<VertData, TriData>::TopoCache::TopoCache(Mesh *owner) : mesh(owner)
{
	init();
}

// support structure for cache construction
struct TopoEdgePrototype {
	uint vid;
	ShortVec<Tptr, 2> tris;
	TopoEdgePrototype() {}
	TopoEdgePrototype(uint v) : vid(v) {}
};

inline TopoEdgePrototype& getTopoEdgePrototype(uint a, uint b, std::vector< ShortVec<TopoEdgePrototype, 8> > &prototypes)
{
	uint N = prototypes[a].size();
	for (uint i = 0; i<N; i++) {
		if (prototypes[a][i].vid == b)
			return prototypes[a][i];
	}
	prototypes[a].push_back(TopoEdgePrototype(b));
	return prototypes[a][N];
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::init()
{
	// first lay out vertices
	std::vector<Vptr> temp_verts(mesh->verts.size()); // need temp. reference
	for (uint i = 0; i<mesh->verts.size(); i++) {
		Vptr vert = verts.alloc(); // cache.verts.alloc()
		vert->ref = i;
		temp_verts[i] = vert;
	}

	// We need to still do the following
	//  * Generate TopoTris
	//  * Generate TopoEdges
	// ---- Hook up references between
	//  * Triangles and Vertices
	//  * Triangles and Edges
	//  * Vertices and Edges

	// We handle two of these items in a pass over the triangles,
	//  * Generate TopoTris
	//  * Hook up Triangles and Vertices
	// building a structure to handle the edges as we go:
	std::vector< ShortVec<TopoEdgePrototype, 8> > edgeacc(mesh->verts.size());
	for (uint i = 0; i<mesh->tris.size(); i++) {
		Tptr tri = tris.alloc(); // cache.tris.alloc()
		tri->ref = i;
		const Tri &ref_tri = mesh->tris[i];

		// triangles <--> verts
		uint vids[3];
		for (uint k = 0; k<3; k++) {
			uint vid = vids[k] = ref_tri.v[k];
			tri->verts[k] = temp_verts[vid];
			temp_verts[vid]->tris.push_back(tri);
		}
		// then, put these in arbitrary but globally consistent order
		if (vids[0] > vids[1])   std::swap(vids[0], vids[1]);
		if (vids[1] > vids[2])   std::swap(vids[1], vids[2]);
		if (vids[0] > vids[1])   std::swap(vids[0], vids[1]);
		// and accrue in structure
		getTopoEdgePrototype(vids[0], vids[1], edgeacc).tris.push_back(tri);
		getTopoEdgePrototype(vids[0], vids[2], edgeacc).tris.push_back(tri);
		getTopoEdgePrototype(vids[1], vids[2], edgeacc).tris.push_back(tri);
	}

	// Now, we can unpack the edge accumulation to
	//  * Generate TopoEdges
	//  * Hook up Triangles and Edges
	//  * Hook up Vertices and Edges
	for (uint vid0 = 0; vid0 < edgeacc.size(); vid0++)
	{
		// for (TopoEdgePrototype &proto : edgeacc[vid0])
		for (uint ii = 0; ii < edgeacc[vid0].size(); ii++)
		{
			TopoEdgePrototype &proto = edgeacc[vid0][ii];
			uint vid1 = proto.vid;
			Vptr v0 = temp_verts[vid0];
			Vptr v1 = temp_verts[vid1];

			Eptr edge = edges.alloc(); // cache.edges.alloc()
									   // edges <--> verts
			edge->verts[0] = v0;
			v0->edges.push_back(edge);
			edge->verts[1] = v1;
			v1->edges.push_back(edge);
			// edges <--> tris
			// for (Tptr tri : proto.tris)
			for (uint j = 0; j < proto.tris.size(); j++)
			{
				Tptr tri = proto.tris[j];
				edge->tris.push_back(tri);
				for (uint k = 0; k<3; k++) {
					if (v0 != tri->verts[k] && v1 != tri->verts[k]) {
						tri->edges[k] = edge;
						break;
					}
				}
			}
		}
	}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::commit()
{
	// record which vertices are live
	std::vector<bool> live_verts(mesh->verts.size(), false);
	//verts.for_each(	[&](Vptr vert)		// 4
	//				{ // cache.verts
	//					live_verts[vert->ref] = true;
	//				}
	//);
	for (IterPool<TopoVert>::Block *block = verts.block_list; block != NULL; block = block->next)
	{
		live_verts[((Vptr)(block))->ref] = true;
	}

	// record which triangles are live, and record connectivity
	std::vector<bool> live_tris(mesh->tris.size(), false);
	//tris.for_each(	[&](Tptr tri)		// 5
	//				{ // cache.tris
	//					live_tris[tri->ref] = true;
	//					for (uint k = 0; k<3; k++)
	//						mesh->tris[tri->ref].v[k] = tri->verts[k]->ref;
	//				}
	//);
	for (IterPool<TopoTri>::Block *block = tris.block_list; block != NULL; block = block->next)
	{
		live_tris[((Tptr)(block))->ref] = true;
		for (uint k = 0; k < 3; k++)
			mesh->tris[((Tptr)(block))->ref].v[k] = ((Tptr)(block))->verts[k]->ref;
	}

	// compact the vertices and build a remapping function
	std::vector<uint> vmap(mesh->verts.size());
	uint write = 0;
	for (uint read = 0; read < mesh->verts.size(); read++) {
		if (live_verts[read]) {
			vmap[read] = write;
			mesh->verts[write] = mesh->verts[read];
			write++;
		}
		else
			vmap[read] = INVALID_ID;
	}
	mesh->verts.resize(write);

	// rewrite the vertex reference ids
	//verts.for_each(	[&](Vptr vert)		// 6
	//					{ // cache.verts
	//						vert->ref = vmap[vert->ref];
	//					}
	//);
	for (IterPool<TopoVert>::Block *block = verts.block_list; block != NULL; block = block->next) {
		((Vptr)(block))->ref = vmap[((Vptr)(block))->ref];
	}

	std::vector<uint> tmap(mesh->tris.size());
	write = 0;
	for (uint read = 0; read < mesh->tris.size(); read++) {
		if (live_tris[read]) {
			tmap[read] = write;
			mesh->tris[write] = mesh->tris[read];
			for (uint k = 0; k<3; k++)
				mesh->tris[write].v[k] = vmap[mesh->tris[write].v[k]];
			write++;
		}
		else
			tmap[read] = INVALID_ID;
	}
	mesh->tris.resize(write);

	// rewrite the triangle reference ids
	//tris.for_each(	[&](Tptr tri)		// 7
	//				{ // cache.tris
	//					tri->ref = tmap[tri->ref];
	//				}
	//);
	for (IterPool<TopoTri>::Block *block = tris.block_list; block != NULL; block = block->next)
	{
		((Tptr)(block))->ref = tmap[((Tptr)(block))->ref];
	}
}

// support functions for validity check
template<class T, class Container> inline
bool count(const Container &contain, const T &val) {
	uint c = 0;
	// for (const T &t : contain)
	for (uint i = 0; i < contain.size(); i++)
		if (contain[i] == val)
			c++;
	return c;
}

template<class T> inline
bool count2(const T arr[], const T &val) {
	return ((arr[0] == val) ? 1 : 0) + ((arr[1] == val) ? 1 : 0);
}

template<class T> inline
bool count3(const T arr[], const T &val) {
	return ((arr[0] == val) ? 1 : 0) + ((arr[1] == val) ? 1 : 0) + ((arr[2] == val) ? 1 : 0);
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::TopoCache::isValid()
{
	std::set<Vptr> vaddr;
	std::set<Eptr> eaddr;
	std::set<Tptr> taddr;
	//verts.for_each([&vaddr](Vptr v) { vaddr.insert(v); });		// 8
	for (IterPool<TopoVert>::Block *block = verts.block_list; block != NULL; block = block->next)
	{
		vaddr.insert(((Vptr)(block)));
	}

	//edges.for_each([&eaddr](Eptr e) { eaddr.insert(e); });		// 9
	for (IterPool<TopoEdge>::Block *block = edges.block_list; block != NULL; block = block->next)
	{
		eaddr.insert(((Eptr)(block)));
	}

	//tris.for_each([&taddr](Tptr t) { taddr.insert(t); });			// 10
	for (IterPool<TopoTri>::Block *block = tris.block_list; block != NULL; block = block->next)
	{
		taddr.insert(((Tptr)(block)));
	}

	// check verts
	//verts.for_each([&](Vptr v)		// 11
	//				{
	//					ENSURE(v->ref < mesh->verts.size());
	//					// make sure each edge pointer goes somewhere and that
	//					// the pointed-to site also points back correctly
	//					for (Eptr e : v->edges) {
	//						ENSURE(eaddr.count(e) > 0); // pointer is good
	//						ENSURE(count2(e->verts, v) == 1); // back-pointer is good
	//					}
	//					for (Tptr t : v->tris) {
	//						ENSURE(taddr.count(t) > 0);
	//						ENSURE(count3(t->verts, v) == 1);
	//					}
	//				}
	//);
	for (IterPool<TopoVert>::Block *block = verts.block_list; block != NULL; block = block->next)
	{
		ENSURE(((Vptr)(block))->ref < mesh->verts.size());
		// make sure each edge pointer goes somewhere and that the pointed-to site also points back correctly
		// for (Eptr e : ((Vptr)(block))->edges)
		for (uint i = 0; i < ((Vptr)(block))->edges.size(); i++)
		{
			Eptr e = ((Vptr)(block))->edges[i];
			ENSURE(eaddr.count(e) > 0); // pointer is good
			ENSURE(count2(e->verts, ((Vptr)(block))) == 1); // back-pointer is good
		}
		// for (Tptr t : ((Vptr)(block))->tris)
		for (uint i = 0; i < ((Vptr)(block))->tris.size(); i++)
		{
			Tptr t = ((Vptr)(block))->tris[i];
			ENSURE(taddr.count(t) > 0);
			ENSURE(count3(t->verts, ((Vptr)(block))) == 1);
		}
	}

	// check edges
	//edges.for_each([&](Eptr e)		// 12
	//				{
	//					// check for non-degeneracy
	//					ENSURE(e->verts[0] != e->verts[1]);
	//					for (uint k = 0; k<2; k++) {
	//						Vptr v = e->verts[k];
	//						ENSURE(vaddr.count(v) > 0);
	//						ENSURE(count(v->edges, e) == 1);
	//					}
	//					for (Tptr t : e->tris) {
	//						ENSURE(taddr.count(t) > 0);
	//						ENSURE(count3(t->edges, e) == 1);
	//					}
	//				}
	//);

	for (IterPool<TopoEdge>::Block *block = edges.block_list; block != NULL; block = block->next)
	{
		// func((T*)(block));
		// check for non-degeneracy
		ENSURE(block->verts[0] != block->verts[1]);
		for (uint k = 0; k < 2; k++) {
			Vptr v = block->verts[k];
			ENSURE(vaddr.count(v) > 0);
			ENSURE(count(v->edges, block) == 1);
		}
		// for (Tptr t : block->tris)
		for (uint i = 0; i < block->tris.size(); i++)
		{
			Tptr t = block->tris[i];
			ENSURE(taddr.count(t) > 0);
			ENSURE(count3(t->edges, block) == 1);
		}
	}

	// check triangles
	//tris.for_each(	[&](Tptr t)			// 13
	//				{
	//					// check for non-degeneracy
	//					ENSURE(t->verts[0] != t->verts[1] && t->verts[1] != t->verts[2]	&& t->verts[0] != t->verts[2]);
	//					for (uint k = 0; k<3; k++) {
	//						Vptr v = t->verts[k];
	//						ENSURE(vaddr.count(v) > 0);
	//						ENSURE(count(v->tris, t) == 1);

	//						Eptr e = t->edges[k];
	//						ENSURE(eaddr.count(e) == 1);
	//						ENSURE(count(e->tris, t) == 1);

	//						// also need to ensure that the edges are opposite the vertices as expected
	//						Vptr v0 = e->verts[0];
	//						Vptr v1 = e->verts[1];
	//						ENSURE((v0 == t->verts[(k + 1) % 3] && v1 == t->verts[(k + 2) % 3])	|| (v0 == t->verts[(k + 2) % 3] && v1 == t->verts[(k + 1) % 3]));
	//					}
	//				}
	//);

	for (IterPool<TopoTri>::Block *block = tris.block_list; block != NULL; block = block->next)
	{
		// func((T*)(block));
		// check for non-degeneracy
		ENSURE(block->verts[0] != block->verts[1] && block->verts[1] != block->verts[2] && block->verts[0] != block->verts[2]);
		for (uint k = 0; k < 3; k++) {
			Vptr v = block->verts[k];
			ENSURE(vaddr.count(v) > 0);
			ENSURE(count(v->tris, block) == 1);

			Eptr e = block->edges[k];
			ENSURE(eaddr.count(e) == 1);
			ENSURE(count(e->tris, block) == 1);

			// also need to ensure that the edges are opposite the vertices as expected
			Vptr v0 = e->verts[0];
			Vptr v1 = e->verts[1];
			ENSURE((v0 == block->verts[(k + 1) % 3] && v1 == block->verts[(k + 2) % 3])
				|| (v0 == block->verts[(k + 2) % 3] && v1 == block->verts[(k + 1) % 3]));
		}
	}

	return true;
}

std::ostream& operator<<(std::ostream &out, const TopoVert& vert)
{
	out << "ref(" << vert.ref << ") "
		<< "e(" << vert.edges.size() << "):";
	// for (Eptr e : vert.edges)
	for (uint i = 0; i < vert.edges.size(); i++)
		out << vert.edges[i] << ";";
	out << " " << "t(" << vert.tris.size() << "):";
	// for (Tptr t : vert.tris)
	for (uint i = 0; i < vert.tris.size(); i++)
		out << vert.tris[i] << ";";
	return out;
}

std::ostream& operator<<(std::ostream &out, const TopoEdge& edge)
{
	out << "v(2):" << edge.verts[0] << "(" << edge.verts[0]->ref << ");" << edge.verts[1] << "(" << edge.verts[1]->ref << ");";
	out << " "	<< "t(" << edge.tris.size() << "):";
	// for (Tptr t : edge.tris)
	for (uint i = 0; i < edge.tris.size(); i++)
		out << edge.tris[i] << ";";
	return out;
}

std::ostream& operator<<(std::ostream &out, const TopoTri& tri)
{
	out << "ref(" << tri.ref << ") ";
	out << "v(3):" << tri.verts[0] << "(" << tri.verts[0]->ref << ");"
		<< tri.verts[1] << "(" << tri.verts[1]->ref << ");"
		<< tri.verts[2] << "(" << tri.verts[2]->ref << ");";
	out << " ";
	out << "e(3):" << tri.edges[0] << ";" << tri.edges[1] << ";" << tri.edges[2] << ";";
	return out;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::print()
{
	using std::cout;
	using std::endl;

	cout << "dumping remeshing cache for debug..." << endl;
	cout << "TRIS" << endl;
	int tri_count = 0;
	//tris.for_each(	[&](Tptr t)			// 14
	//				{
	//					cout << " " << t << ": " << *t << endl;
	//					tri_count++;
	//				}
	//);

	for (IterPool<TopoTri>::Block *block = tris.block_list; block != NULL; block = block->next)
	{
		//func((T*)(block));
		{
			cout << " " << block << ": " << *block << endl;
			tri_count++;
		}
	}


	cout << "There were " << tri_count << " TRIS" << endl;
	cout << "EDGES" << endl;
	int edge_count = 0;
	//edges.for_each(	[&](Eptr e)			// 15
	//				{
	//					cout << " " << e << ": " << endl;
	//					cout << "  v " << e->verts[0] << "; "
	//						<< e->verts[1] << endl;
	//					cout << "  t (" << e->tris.size() << ")" << endl;
	//					for (Tptr t : e->tris)
	//						cout << "    " << t << endl;
	//					edge_count++;
	//				}
	//);

	for (IterPool<TopoEdge>::Block *block = edges.block_list; block != NULL; block = block->next)
	{
		//func((T*)(block));
		cout << " " << block << ": " << endl;
		cout << "  v " << block->verts[0] << "; " << block->verts[1] << endl;
		cout << "  t (" << block->tris.size() << ")" << endl;
		// for (Tptr t : block->tris)
		for (uint i = 0; i < block->tris.size(); i++)
			cout << "    " << block->tris[i] << endl;
		edge_count++;
	}

	cout << "There were " << edge_count << " EDGES" << endl;
	cout << "VERTS" << endl;
	int vert_count = 0;
	//verts.for_each(	[&](Vptr v)				 // 16
	//				{
	//					cout << " " << v << ": ref(" << v->ref << ")" << endl;
	//					cout << "  e (" << v->edges.size() << ")" << endl;
	//					for (Eptr e : v->edges)
	//						cout << "    " << e << endl;
	//					cout << "  t (" << v->tris.size() << ")" << endl;
	//					for (Tptr t : v->tris)
	//						cout << "    " << t << endl;
	//					vert_count++;
	//				}
	//);

	for (IterPool<TopoVert>::Block *block = verts.block_list; block != NULL; block = block->next)
	{
		//func((T*)(block));
		cout << " " << block << ": ref(" << block->ref << ")" << endl;
		cout << "  e (" << block->edges.size() << ")" << endl;
		// for (Eptr e : block->edges)
		for (uint i = 0; i < block->edges.size(); i++)
			cout << "    " << block->edges[i] << endl;
		cout << "  t (" << block->tris.size() << ")" << endl;
		// for (Tptr t : block->tris)
		for (uint i = 0; i < block->tris.size(); i++)
			cout << "    " << block->tris[i] << endl;
		vert_count++;
	}

	cout << "There were " << vert_count << " VERTS" << endl;
}

// constructors
template<class VertData, class TriData>
Mesh<VertData, TriData>::Mesh(Mesh &cp) : tris(cp.tris), verts(cp.verts) {}

template<class VertData, class TriData>
Mesh<VertData, TriData>::Mesh(const RawMesh<VertData, TriData> &raw) : tris(raw.triangles.size()), verts(raw.vertices)
{
	// fill out the triangles
	for (uint i = 0; i < raw.triangles.size(); i++) {
		tris[i].data = raw.triangles[i];
		tris[i].a = raw.triangles[i].a;
		tris[i].b = raw.triangles[i].b;
		tris[i].c = raw.triangles[i].c;
	}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::operator=(Mesh &src)
{
	tris = src.tris;
	verts = src.verts;
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::valid() const
{
	for (uint i = 0; i < verts.size(); i++) {
		if (!std::isfinite(verts[i].pos.x) || !std::isfinite(verts[i].pos.y) || !std::isfinite(verts[i].pos.z)) {
			std::ostringstream message;
			message << "vertex #" << i << " has non-finite coordinates: " << verts[i].pos;
			CORK_ERROR(message.str());
			return false;
		}
	}

	for (uint i = 0; i < tris.size(); i++) {
		if (tris[i].a >= verts.size() || tris[i].b >= verts.size() || tris[i].c >= verts.size())
		{
			std::ostringstream message;
			message << "triangle #" << i << " should have indices in the range 0 to " << (verts.size() - 1)
				<< ", but it has invalid indices: "	<< tris[i].a << ", " << tris[i].b << ", " << tris[i].c;
			CORK_ERROR(message.str());
			return false;
		}
	}
	return true;
}

template<class VertData, class TriData>
RawMesh<VertData, TriData> Mesh<VertData, TriData>::raw() const
{
	RawMesh<VertData, TriData> result;
	result.vertices = verts;
	result.triangles.resize(tris.size());
	for (uint i = 0; i < tris.size(); i++) {
		result.triangles[i] = tris[i].data;
		result.triangles[i].a = tris[i].a;
		result.triangles[i].b = tris[i].b;
		result.triangles[i].c = tris[i].c;
	}
	return result;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::disjointUnion(const Mesh &cp)
{
	uint oldVsize = verts.size();
	uint oldTsize = tris.size();
	uint cpVsize = cp.verts.size();
	uint cpTsize = cp.tris.size();
	uint newVsize = oldVsize + cpVsize;
	uint newTsize = oldTsize + cpTsize;

	std::vector<int> v_remap(cpVsize); // oh this is obvious...
	verts.resize(newVsize);
	tris.resize(newTsize);

	for (uint i = 0; i < cpVsize; i++)
		verts[oldVsize + i] = cp.verts[i];

	for (uint i = 0; i < cpTsize; i++) {
		Tri &tri = tris[oldTsize + i];
		tri = cp.tris[i];
		tri.a += oldVsize;
		tri.b += oldVsize;
		tri.c += oldVsize;
	}
}

// Picking.
// Dumb Implementation just passes over all triangles w/o any precomputed acceleration structure
template<class VertData, class TriData>
typename Mesh<VertData, TriData>::Isct
Mesh<VertData, TriData>::pick(Ray3d ray)
{
	Isct result;
	result.ray = ray;
	result.exists = false;

	double mint = DBL_MAX;

	// pass all triangles over ray
	for (uint i = 0; i < tris.size(); i++) {
		const Tri  &tri = tris[i];

		uint   a = tri.a;
		uint   b = tri.b;
		uint   c = tri.c;
		Vec3d va = verts[a].pos;
		Vec3d vb = verts[b].pos;
		Vec3d vc = verts[c].pos;
		// normalize vertex order (to prevent leaks)
		if (a > b) { std::swap(a, b); std::swap(va, vb); }
		if (b > c) { std::swap(b, c); std::swap(vb, vc); }
		if (a > b) { std::swap(a, b); std::swap(va, vb); }

		double t;
		Vec3d  bary;
		if (isct_ray_triangle(ray, va, vb, vc, &t, &bary)) {
			if (t > 0 && t < mint) {
				result.exists = true;
				mint = t;
				result.tri_id = i;
				result.isct = ray.p + t * ray.r;
				result.bary = bary;
			}
		}
	}
	return result;
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::isClosed()
{
	EGraphCache<int> chains = createEGraphCache<int>();
	//chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry)		// 17
	//				{
	//					entry.data = 0;
	//				}
	//);

	for (uint i = 0; i < chains.skeleton.size(); i++)
	{
		for (uint j = 0; j < chains.skeleton[i].size(); j++)
			chains.skeleton[i][j].data = 0;
	}

	// count up how many times each edge is encountered in one orientation vs. the other
	// for (Tri &tri : tris)
	for (uint i = 0; i < tris.size(); i++)
	{
		Tri& tri = tris[i];
		chains(tri.a, tri.b).data++;
		chains(tri.b, tri.a).data--;

		chains(tri.b, tri.c).data++;
		chains(tri.c, tri.b).data--;

		chains(tri.c, tri.a).data++;
		chains(tri.a, tri.c).data--;
	}
	// now go through and see if any of these are non-zero
	bool closed = true;
	//chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry)		// 18
	//				{
	//					if (entry.data != 0)
	//						closed = false;
	//				}
	//);

	for (uint i = 0; i < chains.skeleton.size(); i++)
	{
		for (uint j = 0; j < chains.skeleton[i].size(); j++)
			if (chains.skeleton[i][j].data != 0)
				closed = false;
	}

	return closed;
}

static inline bool contains(const ShortVec<uint, 8> &list, uint item)
{
	// for (uint k : list)
	for (uint i = 0; i < list.size(); i++)
		if (list[i] == item)
			return true;
	return false;
}

template<class VertData, class TriData>
typename Mesh<VertData, TriData>::NeighborCache
Mesh<VertData, TriData>::createNeighborCache()
{
	NeighborCache result;
	result.skeleton.resize(verts.size());

	for (uint tid = 0; tid < tris.size(); tid++)
	{
		const Tri &tri = tris[tid];

		result(tri.a, tri.b).tids.push_back(tid);
		result(tri.b, tri.a).tids.push_back(tid);

		result(tri.a, tri.c).tids.push_back(tid);
		result(tri.c, tri.a).tids.push_back(tid);

		result(tri.b, tri.c).tids.push_back(tid);
		result(tri.c, tri.b).tids.push_back(tid);
	}

	return result;
}

// This f_unction signature is an amazing disaster...
#ifdef _WIN32
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData, TriData>::EGraphCache<Edata>
#else
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData, TriData>::template EGraphCache<Edata>
#endif
Mesh<VertData, TriData>::createEGraphCache()
{
	EGraphCache<Edata> result;
	result.skeleton.resize(verts.size());

	for (uint tid = 0; tid < tris.size(); tid++)
	{
		const Tri &tri = tris[tid];

		result(tri.a, tri.b).tids.push_back(tid);
		result(tri.b, tri.a).tids.push_back(tid);

		result(tri.a, tri.c).tids.push_back(tid);
		result(tri.c, tri.a).tids.push_back(tid);

		result(tri.b, tri.c).tids.push_back(tid);
		result(tri.c, tri.b).tids.push_back(tid);
	}

	return result;
}

template<class VertData, class TriData>
std::vector<uint> Mesh<VertData, TriData>::getComponentIds()
{
	UnionFind uf(verts.size());
	// for (const Tri &tri : tris)
	for (uint i = 0; i < tris.size(); i++)
	{
		const Tri& tri = tris[i];
		uf.unionIds(tri.a, tri.b);
		uf.unionIds(tri.a, tri.c);
	}
	return uf.dump();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////// mesh.remesh.tpp

enum EdgeRemeshOperation {
	EDGE_SPLIT,
	EDGE_COLLAPSE,
	EDGE_NOTHING
};

struct RemeshEdgeAuxiliary {
	double                  score;      // priority for queue
	EdgeRemeshOperation     op;         // operation to perform
};

template<class VertData, class TriData>
struct Mesh<VertData, TriData>::RemeshScratchpad {
	TopoCache                               cache;
	MemPool<RemeshEdgeAuxiliary>            edge_data;
	std::set< std::pair<double, Eptr> >     queue;

	RemeshScratchpad(Mesh<VertData, TriData> *mesh) : cache(mesh) {}
};


template<class VertData, class TriData>
void Mesh<VertData, TriData>::remesh()
{
	//std::cout << "remesh initial tri count: " << tris.size() << std::endl;
	if (verts.size() == 0)   return; // pathology guard
							 
	RemeshScratchpad scratchpad(this);// create a scratchpad and set it up

	// compute which vertices are boundary or not
	// for (VertData &v : verts)
	for (uint i = 0; i < verts.size(); i++)
	{
		verts[i].manifold = true;
	}
	//scratchpad.cache.edges.for_each([this, &scratchpad](Eptr edge) {		// 19
	//	if (edge->tris.size() != 2) {
	//		verts[edge->verts[0]->ref].manifold = false;
	//		verts[edge->verts[1]->ref].manifold = false;
	//	}
	//});

	for (IterPool<TopoEdge>::Block *block = scratchpad.cache.edges.block_list; block != NULL; block = block->next)
	{
		//func((T*)(block));
		if (block->tris.size() != 2) {
			verts[block->verts[0]->ref].manifold = false;
			verts[block->verts[1]->ref].manifold = false;
		}
	}

	// Then, we set up the priority queue (allocating auxiliary data as we go)
	//scratchpad.cache.edges.for_each([this, &scratchpad](Eptr edge)			// 20
	//								{
	//									edge->data = scratchpad.edge_data.alloc();
	//									scoreAndEnqueue(scratchpad.queue, edge);
	//								}
	//);

	for (IterPool<TopoEdge>::Block *block = scratchpad.cache.edges.block_list; block != NULL; block = block->next)
	{
		//func((T*)(block));
		block->data = scratchpad.edge_data.alloc();
		scoreAndEnqueue(scratchpad.queue, block);
	}

	int cutoff = 10000;
	// now let's go into a loop pulling work off of the queue
	while (scratchpad.queue.size() > 0) {
		// pop the end of the queue
		set< std::pair<double, Eptr> >::iterator it_top = scratchpad.queue.end();

		it_top--;
		Eptr top_edge = it_top->second;
		scratchpad.queue.erase(it_top);

		// Which operation should be performed to this edge?
		EdgeRemeshOperation op = reinterpret_cast<RemeshEdgeAuxiliary*>(top_edge->data)->op;
		if (op == EDGE_SPLIT) {
			//std::cout << "edge split" << std::endl;
			edgeSplit(scratchpad, top_edge);
		}
		else if (op == EDGE_COLLAPSE) {
			//std::cout << "edge collapse" << std::endl;
			edgeCollapse(scratchpad, top_edge, true);
		} // else do nothing (should have score <= 0.0 then...)
		cutoff--;
		if (cutoff < 0)
			break;
	}

	// finally, we commit all of the results of this remeshing operation causing the data storage to defrag and clean out dead items
	scratchpad.cache.commit();

	//std::cout << "remesh final tri count: " << tris.size() << std::endl;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::scoreAndEnqueue(std::set< std::pair<double, Eptr> > &queue, Eptr edge)
{
	double score = computeEdgeScore(edge);
	if (score > 0.0) // only enqueue edges with actual work to do
		queue.insert(std::make_pair(score, edge));
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::dequeue(std::set< std::pair<double, Eptr> > &queue, Eptr edge)
{
	double score = reinterpret_cast<RemeshEdgeAuxiliary*>(edge->data)->score;
	if ((queue.find(std::make_pair(score, edge))) != queue.end())
		queue.erase(it);
}

template<class VertData, class TriData>
double Mesh<VertData, TriData>::computeEdgeScore(Eptr edge) {
	double edge_length;
	double min_angle = 360.0; // clearly overkill
	double max_angle = 0.0;

	edgeNeighborhood_GY(edge, edge_length, min_angle, max_angle);

	RemeshEdgeAuxiliary *edge_aux = reinterpret_cast<RemeshEdgeAuxiliary*>(edge->data);

	// extract violation quantities:
	//      if any of these are positive,
	//      then an operation should be performed to this edge
	double excess_length = edge_length - remesh_options.maxEdgeLength;
	double deficient_length = remesh_options.minEdgeLength - edge_length;
	double excess_angle = max_angle - remesh_options.maxAngle;
	double deficient_angle = remesh_options.minAngle - max_angle;
	if (excess_length <= 0.0 && deficient_length <= 0.0 && excess_angle <= 0.0 && deficient_angle <= 0.0)
	{
		// no op on edge
		edge_aux->score = -1.0;
		edge_aux->op = EDGE_NOTHING;
		return edge_aux->score;
	}

	// rescale violation quantities:
	//      make quantities more commensurate for use as priorities 
	excess_length /= remesh_options.maxEdgeLength;
	deficient_length /= remesh_options.minEdgeLength;
	excess_angle /= (180.0 - remesh_options.maxAngle);
	deficient_angle /= remesh_options.minAngle;
	edge_aux->score = -1.0;
	if (excess_length > edge_aux->score) {
		edge_aux->score = excess_length;
		edge_aux->op = EDGE_SPLIT;
	}
	if (deficient_length > edge_aux->score) {
		edge_aux->score = deficient_length;
		edge_aux->op = EDGE_COLLAPSE;
	}
	if (excess_angle > edge_aux->score) {
		edge_aux->score = excess_angle;
		edge_aux->op = EDGE_SPLIT;
	}
	if (deficient_angle > edge_aux->score) {
		edge_aux->score = deficient_angle;
		edge_aux->op = EDGE_COLLAPSE;
	}

	return edge_aux->score;
}

// erase the first occurrence of a value from a short vector,
// !!! KNOWING that it MUST BE PRESENT
template<class T, uint LEN> inline
void remove(ShortVec<T, LEN> &vec, T val)
{
	uint last_i = vec.size() - 1;
	for (uint i = 0; i<last_i; i++) {
		if (vec[i] == val) {
			std::swap(vec[i], vec[last_i]);
			break;
		}
	}
	vec.resize(last_i);
}

// erase the first occurrence of a value, IF THE VALUE OCCURS
template<class T, uint LEN> inline
void maybe_erase(ShortVec<T, LEN> &vec, T val)
{
	uint N = vec.size();
	for (uint i = 0; i<N; i++) {
		if (vec[i] == val) {
			std::swap(vec[i], vec[N - 1]);
			vec.resize(N - 1);
			break;
		}
	}
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::populateTriFromTopoTri(Tptr tri) {
	Tri &tri_ref = tris[tri->ref];
	for (uint k = 0; k<3; k++)
		tri_ref.v[k] = tri->verts[k]->ref;
}

template<class VertData, class TriData>
Eptr Mesh<VertData, TriData>::allocateRemeshEdge(
	RemeshScratchpad &scratchpad
	) {
	Eptr        e = scratchpad.cache.newEdge();
	e->data = scratchpad.edge_data.alloc();
	return      e;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::deallocateRemeshEdge(RemeshScratchpad &scratchpad, Eptr e)
{
	dequeue(scratchpad.queue, e);
	scratchpad.cache.freeEdge(e);
}

/*
Strategy:
1. Construct new geometry in parallel to old geometry w/o tampering with old geometry in any way.
2. Fill out data for new geometry
3. Determine whether or not to commit or abort the operation
4. Connect new geometry to old geometry simultaneously to old geometry.
5. Delete the old geometry appropriately
*/

struct EdgeWedge {
	Eptr e0, e1;

	EdgeWedge() : e0(NULL), e1(NULL) {}
	bool full() const { return e0 && e1; }
	Eptr one()  const { return (e0) ? e0 : e1; }
};

struct TriWedge {
	Tptr t0, t1;

	TriWedge() : t0(NULL), t1(NULL) {}
	bool full() const { return t0 && t1; }
	Tptr one()  const { return (t0) ? t0 : t1; }
};

struct VptrRemap {
	Vptr v0, v1, v_merged;
	VptrRemap(Vptr v0_, Vptr v1_, Vptr v_merged_) : v0(v0_), v1(v1_), v_merged(v_merged_) {}
	Vptr operator[](Vptr key) {
		if (key == v0 || key == v1)          return v_merged;
		else                                return key;
	}
};

template<class T>
struct PtrRemap {
	std::map<T*, T*> remap;

	void set(T* key, T* val) {
		remap[key] = val;
	}

	T* operator[](T* key)
	{
		std::map<T*, T*>::iterator it = remap.find(key);
		if (it == remap.end())
			return key;
		else
			return it->second;
	}
};

template<class VertData, class TriData>
void Mesh<VertData, TriData>::edgeCollapse(RemeshScratchpad &scratchpad, Eptr e_collapse, bool collapsing_tetrahedra_disappear)
{
	/*
	*
	*  Two cases: both with and without a triangle filling the triangular arrangement (i.e. wedge) of edges
	*
	*                     *
	*                   /   \
	*                 / # # # \
	*               /  # # # #  \
	*             / # # # # # # # \
	*           /  # # # # # # # #  \
	*         / # # # # # # # # # # # \
	*  vid0 *---------------------------* vid1
	*                eid_collapse
	*
	*                     *
	*                   /   \
	*                 /       \
	*               /           \
	*             /               \
	*           /                   \
	*         /                       \
	*  vid0 *---------------------------* vid1
	*                eid_collapse
	*
	*
	*  Additionally, we need to worry about how tetrahedral structures can collapse.
	*  If the collapse will identify the two non-collapsing faces, then we call this structure a triangle wedge
	*
	*                     *
	*                    /|\
	*                   / | \
	*                  /  |  \
	*                 / # | # \
	*                / #  |  # \
	*               /   # | # # \
	*              / # #  *  # # \
	*             /  #  /   \  #  \
	*            / #  /       \  # \
	*           /   /           \   \
	*          /  /               \  \
	*         / /                   \ \
	*        //                       \\
	*  vid0 *---------------------------* vid1
	*                 eid_collapse
	*
	In summary, there are three features we should be interested in:
	--  hollow edge wedges that will be collapsed
	--  filled edge wedges that will be collapsed
	--  triangle wedges that will be collapsed

	We can identify all edge wedges (hollow or filled) via:
	JOIN FILTER(EDGES(A,X), X!=B) WITH
	FILTER(EDGES(B,Y), Y!=A)
	WHERE X=Y
	We can identify all triangle wedges via:
	JOIN FILTER(TRIANGLES(A,X,Y), X!=B && Y!=B) WITH
	FILTER(TRIANGLES(B,Z,W), Z!=A && W!=A)
	WHERE (X,Y)=(Z,W) or (X,Y)=(W,Z)
	*/

	/* Identify all the following components
	* and separate them according to type:
	*  -   e_collapse:     the collapsing edge
	*  -   tris_collapse:  the collapsing triangles
	*  -   v0, v1:         the merging vertices
	*  -   edge_wedges:    the merging edges
	*  -   tri_wedges:     the merging (or disappearing) triangles
	*  -   edges_moving:   the persisting, but moving edges
	*  -   tris_moving:    the persisting, but moving triangles
	*/
	ShortVec<Tptr, 2>       tris_collapse = e_collapse->tris;

	Vptr                    v0 = e_collapse->verts[0];
	Vptr                    v1 = e_collapse->verts[1];

	ShortVec<EdgeWedge, 2>  edge_wedges;
	ShortVec<Eptr, 10>      edges_moving;

	ShortVec<TriWedge, 2>   tri_wedges;
	ShortVec<Tptr, 14>      tris_moving;

	// BUILD all the edge wedges
	std::map<Vptr, EdgeWedge>   edge_w_map;
	// for (Eptr e : v0->edges)
	for (uint ii = 0; ii < v0->edges.size(); ii++)
	{
		Eptr e = v0->edges[ii];
		if (e == e_collapse)
			continue;
		Vptr ev0 = e->verts[0];
		Vptr ev1 = e->verts[1];
		Vptr key = (ev0 != v0) ? ev0 : ev1;
		edge_w_map[key].e0 = e;
	}
	// for (Eptr e : v1->edges)
	for (uint ii = 0; ii < v1->edges.size(); i++)
	{
		Eptr e = v1->edges[ii];
		if (e == e_collapse) continue;
		Vptr ev0 = e->verts[0];
		Vptr ev1 = e->verts[1];
		Vptr key = (ev0 != v1) ? ev0 : ev1;
		edge_w_map[key].e1 = e;
	}
	for (uint ii = 0; ii < edge_w_map.size())
	{
		const 
		if (pair.second.full())
			edge_wedges.push_back(pair.second);
		else
			edges_moving.push_back(pair.second.one());
	}

	// BUILD all the triangle wedges
	std::map<Eptr, TriWedge> tri_w_map;
	// for (Tptr t : v0->tris)
	for (uint ii = 0; ii < v0->tris.size(); ii++)
	{
		Tptr t = v0->tris[ii];
		if (t->edges[0] == e_collapse || t->edges[1] == e_collapse || t->edges[2] == e_collapse)
			continue;

		Eptr key = NULL;
		for (uint k = 0; k<3; k++)
			if (t->verts[k] == v0)
				key = t->edges[k];
		ENSURE(key != NULL);
		tri_w_map[key].t0 = t;
	}

	//for (Tptr t : v1->tris)
	for (uint ii = 0; ii < v1->tris.size(); ii++)
	{
		Tptr t = v1->tris[ii];

		if (t->edges[0] == e_collapse || t->edges[1] == e_collapse || t->edges[2] == e_collapse)
			continue;

		Eptr key = NULL;
		for (uint k = 0; k<3; k++)
			if (t->verts[k] == v1)
				key = t->edges[k];
		ENSURE(key != NULL);
		tri_w_map[key].t1 = t;
	}
	for (std::map<Eptr, TriWedge>::iterator Iter_Tri = tri_w_map.begin(); Iter_Tri != tri_w_map.end(); Iter_Tri++)
	{
		if (Iter_Tri->second.full())
			tri_wedges.push_back(Iter_Tri->second);
		else
			tris_moving.push_back(Iter_Tri->second.one());
	}

	/* Next, we need to create the new geometry which will supplant
	* the just enumerated parts:
	*  -   v_merged:       the merged vertex
	*  -   edges_merged:   the merged edges (parallel to edge_wedges)
	*  -   tris_merged:    the merged triangles (parallel to tri_wedges)
	*  -   edges_moved:    the moved edges (parallel to edges_moving)
	*  -   tris_moved:     the moved triangles (parallel to tris_moving)
	*/
	/* As part of this, we will build a record of the remapping to aid us when we start to connect all of this geometry
	*/

	Vptr                    v_merged = scratchpad.cache.newVert();
	ShortVec<Eptr, 2>       edges_merged(edge_wedges.size());
	ShortVec<Eptr, 10>      edges_moved(edges_moving.size());
	ShortVec<Tptr, 2>       tris_merged(tri_wedges.size());
	ShortVec<Tptr, 14>      tris_moved(tris_moving.size());

	VptrRemap               vptr_remap(v0, v1, v_merged);
	PtrRemap<TopoEdge>    eptr_remap;
	PtrRemap<TopoTri>     tptr_remap;

	// Create new edges and enter re-mappings
	eptr_remap.set(e_collapse, NULL); // mark this edge as dying
	for (uint i = 0; i<edge_wedges.size(); i++) {
		edges_merged[i] = allocateRemeshEdge(scratchpad);
		eptr_remap.set(edge_wedges[i].e0, edges_merged[i]);
		eptr_remap.set(edge_wedges[i].e1, edges_merged[i]);
	}
	for (uint i = 0; i<edges_moving.size(); i++) {
		edges_moved[i] = allocateRemeshEdge(scratchpad);
		eptr_remap.set(edges_moving[i], edges_moved[i]);
	}

	// Create new triangles and enter re-mappings
	// for (Tptr t : tris_collapse)
	for (uint i = 0; i < tris_collapse.size(); i++)// mark this triangle as dying	
		tptr_remap.set(tris_collapse[i], NULL);

	for (uint i = 0; i<tri_wedges.size(); i++) {
		if (collapsing_tetrahedra_disappear) {
			tris_merged[i] = NULL;
		}
		else {
			tris_merged[i] = scratchpad.cache.newTri();
		}
		tptr_remap.set(tri_wedges[i].t0, tris_merged[i]);
		tptr_remap.set(tri_wedges[i].t1, tris_merged[i]);
	}
	for (uint i = 0; i<tris_moving.size(); i++) {
		tris_moved[i] = scratchpad.cache.newTri();
		tptr_remap.set(tris_moving[i], tris_moved[i]);
	}

	/* Now we want to connect up all of this new geometry to each other
	* and to the existing geometry,
	*  BUT!!!
	* we also want to be careful not to point any existing geometry
	* at the new pieces.  Before doing that, we want to be able to
	* compute data for the new geometry and confirm or deny that we
	* want to commit this operation. (e.g. check for unwanted collisions)
	*
	* We adopt the following strategy:
	*  -   first point the new triangles at edges and vertices
	*  -   next, point the new edges at vertices; and
	*          point the new edges at new triangles.
	*  -   finally, point the new vertex at
	*  -       the new edges and new triangles
	*
	* If an edge cannot find any valid triangles it is incident to, then it must be marked for deletion, etc.
	* If the merged vertex cannot find any valid tris/edges it is incident to, then it must be marked for deletion too!
	*/

	// First, the triangles
	if (!collapsing_tetrahedra_disappear) {
		for (uint i = 0; i<tri_wedges.size(); i++) {
			Tptr t0 = tri_wedges[i].t0;
			Tptr t_new = tris_merged[i];

			for (uint k = 0; k<3; k++) {
				t_new->verts[k] = vptr_remap[t0->verts[k]];
				t_new->edges[k] = eptr_remap[t0->edges[k]];
			}
			populateTriFromTopoTri(t_new);
		}
	}
	for (uint i = 0; i<tris_moving.size(); i++) {
		Tptr                t_old = tris_moving[i];
		Tptr                t_new = tris_moved[i];

		for (uint k = 0; k<3; k++) {
			t_new->verts[k] = vptr_remap[t_old->verts[k]];
			t_new->edges[k] = eptr_remap[t_old->edges[k]];
		}
		populateTriFromTopoTri(t_new);
	}

	// Next, the edges
	for (uint i = 0; i<edge_wedges.size(); i++) {
		Eptr                e0 = edge_wedges[i].e0;
		Eptr                e1 = edge_wedges[i].e1;
		Eptr                e_new = edges_merged[i];

		// plug in all the valid triangles...
		// for (Tptr t : e0->tris)
		for (uint j = 0; j < e0->tris.size(); j++)
		{
			Tptr t = e0->tris[j];
			t = tptr_remap[t];
			if (t)
				e_new->tris.push_back(t);
		}
		// for (Tptr t : e1->tris)
		for (uint j = 0; j < e1->tris.size(); j++)
		{
			Tptr t = e1->tris[j];
			t = tptr_remap[t];
			if (t)
				e_new->tris.push_back(t);
		}

		if (e_new->tris.size() == 0) { // if there are no parent triangles left	then we need to kill this edge
			eptr_remap.set(e0, NULL);
			eptr_remap.set(e1, NULL);
			deallocateRemeshEdge(scratchpad, e_new);
			edges_merged[i] = NULL;
		}
		else { // otherwise, let's go ahead and finish hooking up this edge
			for (uint k = 0; k<2; k++)
				e_new->verts[k] = vptr_remap[e0->verts[k]];
		}
	}
	for (uint i = 0; i<edges_moving.size(); i++) {
		Eptr                e_old = edges_moving[i];
		Eptr                e_new = edges_moved[i];

		// note: should never have any dead/null triangles
		// for (Tptr t : e_old->tris)
		for (uint j = 0; j < e_old->tris.size(); j++)
		{
			Tptr t = e_old->tris[j];
			t = tptr_remap[t];
			ENSURE(t);
			e_new->tris.push_back(t);
		}
		for (uint k = 0; k<2; k++)
			e_new->verts[k] = vptr_remap[e_old->verts[k]];
	}

	// Finally, the vertex
	{
		// Should do this directly, not via the remap translation.
		// Working via re-maps will lead to duplicates of merged geometry.
		// However, we can exploit the fact that this vertex is unique, and that we already have lists of all the incident geometry

		// for (Tptr t : tris_merged)
		for (uint j = 0; j < tris_merged.size(); j++)
			if (tris_merged[j])
				v_merged->tris.push_back(tris_merged[j]);
		// for (Tptr t : tris_moved) // cannot be dead
		for (uint j = 0; j < tris_moved.size(); j++)
			v_merged->tris.push_back(tris_moved[j]);

		if (v_merged->tris.size() == 0) {
			scratchpad.cache.freeVert(v_merged);
			v_merged = NULL;
		}
		else
		{
			// for (Eptr e : edges_merged)
			for (uint j = 0; j < edges_merged.size(); j++)
				if (edges_merged[j])
					v_merged->edges.push_back(edges_merged[j]);

			// for (Eptr e : edges_moved) // cannot be dead
			for (uint j = 0; j < edges_moved.size(); j++)
				v_merged->edges.push_back(edges_moved[j]);
			// it's impossible to have triangles incident w/o edges too.
			ENSURE(v_merged->edges.size() > 0);
		}
	}

	// OK, here we get to finally compute data for all the new geometry
	// Once we've done that, we can also check to see whether we actually want to commit this operation or not.

	// merge vertices' data
	if (v_merged) { // REMEMBER: vertex could be deleted by now
		VertData            &data_new = verts[v_merged->ref];
		const VertData      &data0 = verts[v0->ref];
		const VertData      &data1 = verts[v1->ref];
		data_new.merge(data0, data1);
	}
	// merge triangles' data
	for (uint i = 0; i<tri_wedges.size(); i++) {
		Tptr                t_new = tris_merged[i];
		// REMEMBER: these triangles could be deleted
		if (!t_new)          continue;

		Tptr t0 = tri_wedges[i].t0;
		Tptr t1 = tri_wedges[i].t1;

		merge_tris(t_new->ref, t0->ref, t1->ref);
	}

	// update moved triangles' data
	for (uint i = 0; i<tris_moving.size(); i++) {
		// NOTE: moved triangles cannot be deleted
		Tptr t_new = tris_moved[i];
		Tri  &tri_new = tris[t_new->ref];

		Tptr t_old = tris_moving[i];
		Tri  &tri_old = tris[t_old->ref];

		move_tri(tri_new, tri_old);
	}

	// Find and Store all of the existing, unchanged edges that border this operation.  We will need these references
	// when we account for changes to edge operation priorities later
	ShortVec<Eptr, 16>      borderEdges;
	// for (const TriWedge &e_wedge : tri_wedges)
	for (uint i = 0; i < tri_wedges.size(); i++)
	{
		Tptr t0 = tri_wedges[i].t0;
		for (uint k = 0; k<3; k++) {
			if (t0->verts[k] == v0) {
				borderEdges.push_back(t0->edges[k]);
				break;
			}
		}
	}
	// for (Tptr t : tris_moved)
	for (uint i = 0; i < tris_moved.size(); i++)
	{
		for (uint k = 0; k<3; k++) {
			if (tris_moved[i]->verts[k] == v_merged) {
				borderEdges.push_back(tris_moved[i]->edges[k]);
				break;
			}
		}
	}

	/* Now that we've got the go ahead, let's hook in the new geometry to
	*  the existing geometry!
	* We can do this in the following order:
	*  -   take all the new edges, and add them to their existing endpoint's
	*      edge list.  (all new edges must have exactly one such endpoint)
	*  -   take all the new triangles, and add them to their two existing
	*      endpoints' and one existing edge's triangle lists.
	*/
	// first, consolidate arrays
	ShortVec<Eptr, 12>      new_edges;
	ShortVec<Tptr, 16>      new_tris;
	{
		//for (Eptr e : edges_merged)
		for (uint i = 0; i < edges_merged.size(); i++)
			if (edges_merged[i])
				new_edges.push_back(edges_merged[i]);

		// for (Eptr e : edges_moved) // cannot be deleted
		for (uint i = 0; i < edges_moved.size(); i++)
			new_edges.push_back(edges_moved[i]);

		// for (Tptr t : tris_merged)
		for (uint i = 0; i < tris_merged.size(); i++)
			if (tris_merged[i])
				new_tris.push_back(tris_merged[i]);

		// for (Tptr t : tris_moved) // cannot be deleted
		for (uint i = 0; i < tris_moved.size(); i++)
			new_tris.push_back(tris_moved[i]);
	}
	// hook up edges to existing geometry
	// for (Eptr edge : new_edges)
	for (uint i = 0; i < new_edges.size(); i++)
	{
		Eptr edge = new_edges[i];
		Vptr ev0 = edge->verts[0];
		Vptr ev1 = edge->verts[1];
		Vptr v_old = (ev0 != v_merged) ? ev0 : ev1;
		v_old->edges.push_back(edge);
	}
	// hook up triangles to existing geometry
	// for (Tptr tri : new_tris)
	for (uint i = 0; i < new_tris.size(); i++)
	{
		Tptr tri = new_tris[i];
		for (uint k = 0; k<3; k++)
		{
			if (tri->verts[k] != v_merged)       continue;
			Eptr e = tri->edges[k];
			Vptr tv0 = e->verts[0];
			Vptr tv1 = e->verts[1];
			e->tris.push_back(tri);
			tv0->tris.push_back(tri);
			tv1->tris.push_back(tri);
			break;
		}
	}

	/* Now, let's kill all the old geometry.  We can use the checklist
	* we built at the begining:
	*  -   e_collapse:     the collapsing edge
	*  -   tris_collapse:  the collapsing triangles
	*  -   v0, v1:         the merging vertices
	*  -   edge_wedges:    the merging edges
	*  -   tri_wedges:     the merging (or disappearing) triangles
	*  -   edges_moving:   the persisting, but moving edges
	*  -   tris_moving:    the persisting, but moving triangles
	*
	* We need to be careful to free this geometry in top down order;
	* starting with the triangles and moving towards the vertices.
	* If we furthermore guarantee that any singular edges or vertices
	* created by a triangle deletion are also deleted, then we can focus all of our attention on just deleting triangles
	*/
	ShortVec<Tptr, 16>      dead_tris;
	ShortVec<Eptr, 16>      dead_edges;
	ShortVec<Vptr, 2>       dead_verts;

	// assemble the list of triangles to kill
	// for (Tptr t : tris_collapse)
	for (uint i = 0; i < tris_collapse.size(); i++)
		dead_tris.push_back(tris_collapse[i]);

	for (uint i = 0; i < tri_wedges.size(); i++)
	{
		dead_tris.push_back(tri_wedges[i].t0);
		dead_tris.push_back(tri_wedges[i].t1);
	}

	// for (Tptr t : tris_moving)
	for (uint i = 0; i < tris_moving.size(); i++)
		dead_tris.push_back(tris_moving[i]);

	// process the list of triangles
	// for (Tptr tri : dead_tris)
	for (uint i = 0; i < dead_tris.size(); i++)
	{
		Tptr tri = dead_tris[i];
		// Let's unhook this triangle from its faces first
		for (uint k = 0; k<3; k++)
		{
			Vptr v = tri->verts[k];
			remove(v->tris, tri);
			if (v->tris.size() == 0)
				dead_verts.push_back(v);

			Eptr e = tri->edges[k];
			remove(e->tris, tri);
			if (e->tris.size() == 0)
				dead_edges.push_back(e);
		}
		// now that we're disconnected, go ahead and jettison the triangle
		scratchpad.cache.freeTri(tri);
	}

	// now, we can process the list of edges
	// for (Eptr edge : dead_edges)
	for (uint i = 0; i < dead_edges.size(); i++)
	{
		Eptr edge = dead_edges[i];
		// Let's unhook this edge from its vertices
		for (uint k = 0; k<2; k++) {
			Vptr v = edge->verts[k];
			remove(v->edges, edge);
			// the triangle removal was enough to determine which vertices should die.
			// re-adding them here would lead to duplicates
		}
		// and then jetisson the edge
		deallocateRemeshEdge(scratchpad, edge);
		// If this edge is in the border edge list, then we need to remove it right away!
		maybe_erase(borderEdges, edge);
	}

	// Finally, polish off by getting rid of any vertices that talked too much
	// for (Vptr vert : dead_verts)
	for (uint i = 0; i < dead_verts.size(); i++)
	{
		Vptr vert = dead_verts[i];
		scratchpad.cache.freeVert(vert);
		if (vert == v_merged)
			v_merged = NULL;
	}

	// We pause a moment here to update the manifoldness of any vertices for which it might have changed
	// ONLY do if the merged vertex is still alive...
	if (v_merged)
	{
		verts[v_merged->ref].manifold = true;
		// for (Eptr e : v_merged->edges)
		for (uint i = 0; i < v_merged->edges.size(); i++)
		{
			Eptr e = v_merged->edges[i];
			if (e->tris.size() != 2)
				verts[v_merged->ref].manifold = false;

			// process neighboring point
			Vptr v = e->verts[0];
			if (v == v_merged)
				v = e->verts[1];
			verts[v->ref].manifold = true;
			// for (Eptr ee : v->edges)
			for (uint j = 0; j < v->edges.size(); j++)
			{
				if (v->edges[j]->tris.size() != 2)
				{
					verts[v->ref].manifold = false;
					break;
				}
			}
		}
	}

	// Before we're completely done, we will go through and adjust priorities for edges which might have been effected by this op.
	// Only explicitly dequeue pre-existing edges we did not delete!
	// for (Eptr e : edges_merged)
	for (uint j = 0; j < edges_merged.size(); j++)
	{
		if (edges_merged[j])
		{ // might be deleted
			scoreAndEnqueue(scratchpad.queue, edges_merged[j]);
		}
	}
	// for (Eptr e : edges_moved)
	for(uint i = 0; i < edges_moved.size(); i++)
	{ // def. not deleted
		scoreAndEnqueue(scratchpad.queue, edges_moved[i]);
	}
	// for (Eptr e : borderEdges)
	for (uint i = 0; i < borderEdges.size(); i++)
	{
		// border edges could have been deleted...
		dequeue(scratchpad.queue, borderEdges[i]);
		scoreAndEnqueue(scratchpad.queue, borderEdges[i]);
	}
	// that should more or less complete an edge collapse
}

/*
Strategy:
1. Construct new geometry in parallel to old geometry w/o
tampering with old geometry in any way.
2. Fill out data for new geometry
3. Determine whether or not to commit or abort the operation
4. Connect new geometry to old geometry simultaneously to old geometry.
5. Delete the old geometry appropriately
*/

template<class VertData, class TriData>
void Mesh<VertData, TriData>::edgeSplit(RemeshScratchpad &scratchpad, Eptr e_split)
{
	/*
	Clean-up picture please...

	*                     *  vs_opp
	*                   /   \
	*                 /       \
	*       e0s     /           \     e1s
	*             /               \
	*           /      ts_orig      \
	*         /                       \
	*    v0 *---------------------------* v1
	*                  e_split
	*                     _
	*                     |
	*                     v
	*
	*                     *   vs_opp
	*                   / | \
	*                 /   |   \
	*        e0s    /     | <----------- es_mid
	*             /       |       \
	*           /         |         \
	*         /  t0s_new  |  t1s_new  \
	*    v0 *-------------*-------------* v1
	*           e0_new  v_new  e1_new
	*

	BEGIN
	new vertex
	-- invoke interpolation callback
	FOR EACH TRIANGLE SUB-PROBLEM:
	new next triangle (from split)
	new prev triangle (from split)
	-- invoke triangle split callback
	delete triangle
	FIXUP all of the Tri-Edges
	*/

	/* Identify all the following components
	* and separate them according to type:
	*  -   e_split:        the edge to be split
	*  -   ts_orig:        the triangles to be split
	*  -   v0, v1:         the two endpoint vertices of the edge
	*  (the following is helper data; do not delete these vertices)
	*  -   vs_opp:         the vertices opposite eid_split for each triangle
	*/

	ShortVec<Tptr, 2>       ts_orig = e_split->tris;

	Vptr                    v0 = e_split->verts[0];
	Vptr                    v1 = e_split->verts[1];

	ShortVec<Vptr, 2> vs_opp(ts_orig.size());
	for (uint i = 0; i<ts_orig.size(); i++) {
		Tptr t_orig = ts_orig[i];
		for (uint k = 0; k<3; k++) {
			if (t_orig->edges[k] == e_split) {
				vs_opp[i] = t_orig->verts[k];
				break;
			}
		}
	}

	/* Next, we need to create the new geometry which will supplant
	* the just enumerated parts:
	*  -   v_new:				the vertex introduced by the split
	*  -   e0_new, e1_new:		the two pieces of e_split
	*  -   t0s_new, t1s_new:    the two pieces of each triangle in tids_orig
	*  -   es_mid:				the new edges splitting each triangle
	*/
	Vptr                    v_new = scratchpad.cache.newVert();
	Eptr                    e0_new = allocateRemeshEdge(scratchpad);
	Eptr                    e1_new = allocateRemeshEdge(scratchpad);
	ShortVec<Tptr, 2>       t0s_new(ts_orig.size());
	ShortVec<Tptr, 2>       t1s_new(ts_orig.size());
	ShortVec<Eptr, 2>       es_mid(ts_orig.size());

	for (uint i = 0; i<ts_orig.size(); i++) {
		t0s_new[i] = scratchpad.cache.newTri();
		t1s_new[i] = scratchpad.cache.newTri();
		es_mid[i] = allocateRemeshEdge(scratchpad);
	}

	/* Now we want to connect up all of this new geometry to each other and to the existing geometry,
	*  BUT!!!
	* we also want to be careful not to point any existing geometry at the new pieces.
	Before doing that, we want to be able to compute data for the new geometry and confirm or deny that we want to commit this operation.
	(e.g. check for unwanted collisions)
	*
	* We adopt the following strategy:
	*  -   first point the new triangles at edges and vertices
	*  -   next, point the new edges at vertices; and
	*          point the new edges at new triangles.
	*  -   finally, point the new vertex at the new edges and new triangles
	*/

	// hook up t0s_new and t1s_new also go ahead and hook up es_mid
	for (uint i = 0; i<ts_orig.size(); i++) 
	{
		Tptr t_orig = ts_orig[i];
		Tptr t0 = t0s_new[i];
		Tptr t1 = t1s_new[i];
		Eptr e_mid = es_mid[i];

		// replace every edge and vertex appropriately for the two variants
		for (uint k = 0; k<3; k++) {
			Vptr            v_orig = t_orig->verts[k];
			Eptr            e_orig = t_orig->edges[k];
			if (v_orig == v0) {
				t0->verts[k] = v_orig;
				t0->edges[k] = e_mid;
				t1->verts[k] = v_new;
				t1->edges[k] = e_orig;
			}
			else if (v_orig == v1) {
				t0->verts[k] = v_new;
				t0->edges[k] = e_orig;
				t1->verts[k] = v_orig;
				t1->edges[k] = e_mid;
			}
			else {
				t0->verts[k] = v_orig;
				t0->edges[k] = e0_new;
				t1->verts[k] = v_orig;
				t1->edges[k] = e1_new;
			}
		}
		populateTriFromTopoTri(t0);
		populateTriFromTopoTri(t1);
		// set up the mid edge from the split
		e_mid->verts[0] = v_new;
		e_mid->verts[1] = vs_opp[i];
		e_mid->tris.resize(2);
		e_mid->tris[0] = t0;
		e_mid->tris[1] = t1;
	}
	// hook up e0_new and e1_new
	e0_new->verts[0] = v0;
	e0_new->verts[1] = v_new;
	e1_new->verts[0] = v_new;
	e1_new->verts[1] = v1;
	for (uint i = 0; i<ts_orig.size(); i++) {
		e0_new->tris.push_back(t0s_new[i]);
		e1_new->tris.push_back(t1s_new[i]);
	}
	// hook up v_new
	v_new->edges.push_back(e0_new);
	v_new->edges.push_back(e1_new);
	for (uint i = 0; i<ts_orig.size(); i++) {
		v_new->edges.push_back(es_mid[i]);
		v_new->tris.push_back(t0s_new[i]);
		v_new->tris.push_back(t1s_new[i]);
	}


	// OK, here we get to finally compute data for all the new geometry
	// Once we've done that, we can also check to see whether we actually want to commit this operation or not.

	// interpolate data onto the new vertex
	{
		VertData            &data_new = verts[v_new->ref];
		const VertData      &data0 = verts[v0->ref];
		const VertData      &data1 = verts[v1->ref];
		data_new.interpolate(data0, data1);
		data_new.manifold = (ts_orig.size() == 2);
	}

	// split triangles' data
	for (uint i = 0; i<ts_orig.size(); i++) {
		Tptr t_orig = ts_orig[i];

		Tptr t0 = t0s_new[i];
		Tptr t1 = t1s_new[i];

		split_tris(t0->ref, t1->ref, t_orig->ref);
	}

	// TODO: COMMIT OPTION will be left as a stub for now!
	if (false) {
		// DESTROY THE GEOMETRY WE JUST CREATED AND RETURN
	}

	// record border edges for later priority updates
	ShortVec<Eptr, 4>       borderEdges;
	// for (Tptr t : ts_orig)
	for (uint i = 0; i < ts_orig.size(); i ++)
	{ // TODO: evacuate all of this to the end...
		Tptr t = ts_orig[i];
		for (uint k = 0; k<3; k++) {
			if (t->edges[k] == e_split)  continue;
			borderEdges.push_back(t->edges[k]);
		}
	}

	/* Now that we've got the go ahead, let's hook in the new geometry to
	*  the existing geometry!
	* We can do this in the following order:
	*  -   take all the new edges, and add them to their existing endpoint's
	*      edge list.  (all new edges must have exactly one such endpoint)
	*  -   take all the new triangles, and add them to their two existing
	*      endpoints' and one existing edge's triangle lists.
	*/

	// add new edges to v0 and v1
	v0->edges.push_back(e0_new);
	v1->edges.push_back(e1_new);

	// now, let's tackle the other edges and triangles in tandem...
	for (uint i = 0; i<ts_orig.size(); i++) {
		Tptr t_orig = ts_orig[i];
		Vptr v_opp = vs_opp[i];

		// add mid edge and two tris to v_opp
		v_opp->edges.push_back(es_mid[i]);
		v_opp->tris.push_back(t0s_new[i]);
		v_opp->tris.push_back(t1s_new[i]);
		// add resp. tris to v0 and v1
		v0->tris.push_back(t0s_new[i]);
		v1->tris.push_back(t1s_new[i]);
		// find the two non-split edges and add resp. tris
		for (uint k = 0; k<3; k++) {
			if (t_orig->verts[k] == v0) {
				Eptr e1 = t_orig->edges[k];
				e1->tris.push_back(t1s_new[i]);
			}
			else if (t_orig->verts[k] == v1) {
				Eptr e0 = t_orig->edges[k];
				e0->tris.push_back(t0s_new[i]);
			}
		}
	}

	/* Now, let's kill all the old geometry.  This consists of:
	*  -   e_split:        the edge to be split
	*  -   ts_orig:        the triangles to be split
	*
	* Luckily, in triangle splitting we know exactly which things must be deleted.  A split cannot make any geometry newly singular.
	*/

	// kill triangles
	// for (Tptr t : ts_orig)
	for (uint i = 0; i < ts_orig.size(); i++)
	{
		Tptr t = ts_orig[i];
		// First, unhook this triangle from its faces
		for (uint k = 0; k<3; k++) {
			Vptr v = t->verts[k];
			remove(v->tris, t);

			Eptr e = t->edges[k];
			remove(e->tris, t);
		}
		// now that we're disconnected, jettison the triangle
		scratchpad.cache.freeTri(t);
	}

	// now, kill the edge that we split
	remove(v0->edges, e_split);
	remove(v1->edges, e_split);
	deallocateRemeshEdge(scratchpad, e_split);

	// recompute edge scores for all edges whose scores might be effected
	// Don't need to dequeue newly created edges...
	scoreAndEnqueue(scratchpad.queue, e0_new);
	scoreAndEnqueue(scratchpad.queue, e1_new);
	// for (Eptr e : es_mid)
	for (uint i = 0; i < es_mid.size(); i++)
	{
		scoreAndEnqueue(scratchpad.queue, es_mid[i]);
	}
	// for (Eptr e : borderEdges)
	for (uint i = 0; i < borderEdges.size(); i++)
	{
		dequeue(scratchpad.queue, borderEdges[i]);
		scoreAndEnqueue(scratchpad.queue, borderEdges[i]);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////// mesh.isct.tpp
#include "bbox.h"
#include "quantization.h"
#include "empty3d.h"

#include "aabvh.h"

#define REAL double
extern "C" {
#include "triangle.h"
}

struct GenericVertType;
struct IsctVertType;
struct OrigVertType;
struct GenericEdgeType;
struct IsctEdgeType;
struct OrigEdgeType;
struct SplitEdgeType;
struct GenericTriType;

struct GluePointMarker;

typedef GenericVertType*    GVptr;
typedef IsctVertType*           IVptr;
typedef OrigVertType*	        OVptr;
typedef GenericEdgeType*    GEptr;
typedef IsctEdgeType*           IEptr;
typedef OrigEdgeType*	        OEptr;
typedef SplitEdgeType*	        SEptr;
typedef GenericTriType*     GTptr;

typedef GluePointMarker*    GluePt;

struct GenericVertType
{
	virtual ~GenericVertType() {}
	Vptr                    concrete;
	Vec3d                   coord;

	bool                    boundary;
	uint                    idx; // temporary for triangulation marshalling

	ShortVec<GEptr, 2>       edges;
};
struct IsctVertType : public GenericVertType
{
	GluePt                  glue_marker;
};
struct OrigVertType : public GenericVertType {};

struct GenericEdgeType
{
	virtual ~GenericEdgeType() {}
	Eptr                    concrete;

	bool                    boundary;
	uint                    idx; // temporary for triangulation marshalling

	GVptr                   ends[2];
	ShortVec<IVptr, 1>      interior;
};
struct IsctEdgeType : public GenericEdgeType
{
public:
	// use to detect duplicate instances within a triangle
	Tptr                    other_tri_key;
};

struct OrigEdgeType : public GenericEdgeType {};
struct SplitEdgeType : public GenericEdgeType {};

struct GenericTriType
{
	Tptr concrete;
	GVptr verts[3];
};

struct GluePointMarker
{
	// list of all the vertices to be glued...
	ShortVec<IVptr, 3>      copies;
	bool                    split_type; // splits are introduced manually, not via intersection and therefore use only e pointer
	bool                    edge_tri_type; // true if edge-tri intersection false if tri-tri-tri
	Eptr                    e;
	Tptr                    t[3];
};

template<uint LEN> inline
IEptr find_edge(ShortVec<IEptr, LEN> &vec, Tptr key)
{
	// for (IEptr ie : vec)
	for (uint i = 0; i < vec.size(); i++)
	{
		IEptr ie = vec[i];
		if (ie->other_tri_key == key)
			return ie;
	}
	return NULL;
}

inline Vptr commonVert(Tptr t0, Tptr t1)
{
	for (uint i = 0; i<3; i++) {
		for (uint j = 0; j<3; j++) {
			if (t0->verts[i] == t1->verts[j])
				return t0->verts[i];
		}
	}
	return NULL;
}

inline bool hasCommonVert(Tptr t0, Tptr t1)
{
	return (t0->verts[0] == t1->verts[0] ||
		t0->verts[0] == t1->verts[1] ||
		t0->verts[0] == t1->verts[2] ||
		t0->verts[1] == t1->verts[0] ||
		t0->verts[1] == t1->verts[1] ||
		t0->verts[1] == t1->verts[2] ||
		t0->verts[2] == t1->verts[0] ||
		t0->verts[2] == t1->verts[1] ||
		t0->verts[2] == t1->verts[2]);
}

inline bool hasCommonVert(Eptr e, Tptr t)
{
	return (e->verts[0] == t->verts[0] ||
		e->verts[0] == t->verts[1] ||
		e->verts[0] == t->verts[2] ||
		e->verts[1] == t->verts[0] ||
		e->verts[1] == t->verts[1] ||
		e->verts[1] == t->verts[2]);
}

inline void disconnectGE(GEptr ge)
{
	ge->ends[0]->edges.erase(ge);
	ge->ends[1]->edges.erase(ge);
	// for (IVptr iv : ge->interior)
	for (uint i = 0; i < ge->interior.size(); i++)
		ge->interior[i]->edges.erase(ge);
}

// should deal with via pointers
template<class VertData, class TriData>
class Mesh<VertData, TriData>::TriangleProblem
{
public:
	TriangleProblem() {}
	~TriangleProblem() {}

	inline void init(IsctProblem *iprob, Tptr t) {
		the_tri = t;
		// extract original edges/verts
		for (uint k = 0; k<3; k++)
			overts[k] = iprob->newOrigVert(the_tri->verts[k]);
		for (uint k = 0; k<3; k++) {
			oedges[k] = iprob->newOrigEdge(the_tri->edges[k],
				overts[(k + 1) % 3],
				overts[(k + 2) % 3]);
		}
	}

private: // may actually not add edge, but instead just hook up endpoint
	inline void addEdge(IsctProblem *iprob, IVptr iv, Tptr tri_key)
	{
		IEptr ie = find_edge(iedges, tri_key);
		if (ie) { // if the edge is already present
			ie->ends[1] = iv;
			iv->edges.push_back(ie);
		}
		else { // if the edge is being added
			ie = iprob->newIsctEdge(iv, tri_key);
			iedges.push_back(ie);
		}
	}
	void addBoundaryHelper(Eptr edge, IVptr iv)
	{
		iv->boundary = true;
		iverts.push_back(iv);
		// hook up point to boundary edge interior!
		for (uint k = 0; k<3; k++) {
			OEptr   oe = oedges[k];
			if (oe->concrete == edge) {
				oe->interior.push_back(iv);
				iv->edges.push_back(oe);
				break;
			}
		}
	}
public:
	// specify reference glue point and edge piercing this triangle.
	IVptr addInteriorEndpoint(IsctProblem *iprob, Eptr edge, GluePt glue)
	{
		IVptr iv = iprob->newIsctVert(edge, the_tri, glue);
		iv->boundary = false;
		iverts.push_back(iv);
		// for (Tptr tri_key : edge->tris)
		for (uint i = 0; i < edge->tris.size(); i++)
		{
			addEdge(iprob, iv, edge->tris[i]);
		}
		return iv;
	}
	// specify the other triangle cutting this one, the edge cut, and the resulting point of intersection
	void addBoundaryEndpoint(IsctProblem *iprob, Tptr tri_key, Eptr edge, IVptr iv)
	{
		iv = iprob->copyIsctVert(iv);
		addBoundaryHelper(edge, iv);
		// handle edge extending into interior
		addEdge(iprob, iv, tri_key);
	}
	IVptr addBoundaryEndpoint(IsctProblem *iprob, Tptr tri_key, Eptr edge, Vec3d coord, GluePt glue)
	{
		IVptr iv = iprob->newSplitIsctVert(coord, glue);
		addBoundaryHelper(edge, iv);
		// handle edge extending into interior
		addEdge(iprob, iv, tri_key);
		return iv;
	}
	// Should only happen for manually inserted split points on edges, not for points computed via intersection...
	IVptr addBoundaryPointAlone(IsctProblem *iprob, Eptr edge, Vec3d coord, GluePt glue)
	{
		IVptr iv = iprob->newSplitIsctVert(coord, glue);
		addBoundaryHelper(edge, iv);
		return iv;
	}
	void addInteriorPoint(IsctProblem *iprob, Tptr t0, Tptr t1, GluePt glue)
	{
		// note this generates wasted re-computation of coordinates 3X
		IVptr iv = iprob->newIsctVert(the_tri, t0, t1, glue);
		iv->boundary = false;
		iverts.push_back(iv);
		// find the 2 interior edges
		// for (IEptr ie : iedges)
		for (uint i = 0; i < iedges.size(); i++)
		{
			IEptr ie = iedges[i];
			if (ie->other_tri_key == t0 || ie->other_tri_key == t1)
			{
				ie->interior.push_back(iv);
				iv->edges.push_back(ie);
			}
		}
	}

	// run after we've accumulated all the elements
	void consolidate(IsctProblem *iprob)
	{
		// identify all intersection edges missing endpoints and check to see if we can assign an original vertex
		// as the appropriate endpoint.
		// for (IEptr ie : iedges)
		for(uint i = 0; i < iedges.size(); i++)
		{
			IEptr ie = iedges[i];
			if (ie->ends[1] == NULL)
			{
				// try to figure out which vertex must be the endpoint...
				Vptr vert = commonVert(the_tri, ie->other_tri_key);
				if (!vert)
				{
					std::cout << "the  edge is " << ie->ends[0] << ",  " << ie->ends[1] << std::endl;
					IVptr iv = dynamic_cast<IVptr>(ie->ends[0]);
					std::cout << "   " << iv->glue_marker->edge_tri_type << std::endl;
					std::cout << "the   tri is " << the_tri << ": " << *the_tri << std::endl;
					std::cout << "other tri is " << ie->other_tri_key << ": " << *(ie->other_tri_key) << std::endl;
					std::cout << "coordinates for triangles" << std::endl;
					std::cout << "the tri" << std::endl;
					for (uint k = 0; k<3; k++)
						std::cout << iprob->vPos(the_tri->verts[k]) << std::endl;
					for (uint k = 0; k<3; k++)
						std::cout << iprob->vPos(ie->other_tri_key->verts[k]) << std::endl;
					std::cout << "degen count:" << Empty3d::degeneracy_count << std::endl;
					std::cout << "exact count: " << Empty3d::exact_count << std::endl;
				}
				ENSURE(vert); // bad if we can't find a common vertex
							  // then, find the corresponding OVptr, and connect
				for (uint k = 0; k<3; k++) {
					if (overts[k]->concrete == vert) {
						ie->ends[1] = overts[k];
						overts[k]->edges.push_back(ie);
						break;
					}
				}
			}
		}

		ENSURE(isValid());
	}

	bool isValid() const {
		ENSURE(the_tri);
		return true;
	}

	void subdivide(IsctProblem *iprob)
	{
		// collect all the points, and create more points as necessary
		ShortVec<GVptr, 7> points;
		for (uint k = 0; k<3; k++) {
			points.push_back(overts[k]);
		}
		// for (IVptr iv : iverts)
		for (uint i = 0; i < iverts.size(); i++)
		{
			points.push_back(iverts[i]);
		}
		for (uint i = 0; i<points.size(); i++)
			points[i]->idx = i;

		// split edges and marshall data for safety, we zero out references to pre-subdivided edges, which may have been destroyed
		ShortVec<GEptr, 8> edges;
		for (uint k = 0; k<3; k++)
		{
			subdivideEdge(iprob, oedges[k], edges);
			oedges[k] = NULL;
		}
		// for (IEptr &ie : iedges)
		for (uint i = 0; i < iedges.size(); i++)
		{
			subdivideEdge(iprob, iedges[i], edges);
			iedges[i] = NULL;
		}
		for (uint i = 0; i<edges.size(); i++)
			edges[i]->idx = i;

		// find 2 dimensions to project onto
		// get normal
		Vec3d normal = cross(overts[1]->coord - overts[0]->coord, overts[2]->coord - overts[0]->coord);
		uint normdim = maxDim(abs(normal));
		uint dim0 = (normdim + 1) % 3;
		uint dim1 = (normdim + 2) % 3;
		double sign_flip = (normal.v[normdim] < 0.0) ? -1.0 : 1.0;

		struct triangulateio in, out;

		/* Define input points. */
		in.numberofpoints = points.size();
		in.numberofpointattributes = 0;
		in.pointlist = new REAL[in.numberofpoints * 2];
		in.pointattributelist = NULL;
		in.pointmarkerlist = new int[in.numberofpoints];
		for (int k = 0; k<in.numberofpoints; k++) {
			in.pointlist[k * 2 + 0] = points[k]->coord.v[dim0];
			in.pointlist[k * 2 + 1] = points[k]->coord.v[dim1] * sign_flip;
			in.pointmarkerlist[k] = (points[k]->boundary) ? 1 : 0;
		}

		/* Define the input segments */
		in.numberofsegments = edges.size();
		in.numberofholes = 0;// yes, zero
		in.numberofregions = 0;// not using regions
		in.segmentlist = new int[in.numberofsegments * 2];
		in.segmentmarkerlist = new int[in.numberofsegments];
		for (int k = 0; k<in.numberofsegments; k++) {
			in.segmentlist[k * 2 + 0] = edges[k]->ends[0]->idx;
			in.segmentlist[k * 2 + 1] = edges[k]->ends[1]->idx;
			in.segmentmarkerlist[k] = (edges[k]->boundary) ? 1 : 0;
		}

		// to be safe... declare 0 triangle attributes on input
		in.numberoftriangles = 0;
		in.numberoftriangleattributes = 0;

		/* set for flags.... */
		out.pointlist = NULL;
		out.pointattributelist = NULL; // not necessary if using -N or 0 attr
		out.pointmarkerlist = NULL;
		out.trianglelist = NULL; // not necessary if using -E
									//out.triangleattributelist = null; // not necessary if using -E or 0 attr
									//out.trianglearealist = // only needed with -r and -a
									//out.neighborlist = null; // only neccesary if -n is used
		out.segmentlist = NULL; // NEED THIS; output segments go here
		out.segmentmarkerlist = NULL; // NEED THIS for OUTPUT SEGMENTS
										 //out.edgelist = null; // only necessary if -e is used
										 //out.edgemarkerlist = null; // only necessary if -e is used

										 // solve the triangulation problem
		char *params = (char*)("pzQYY");
		triangulate(params, &in, &out, NULL);

		if (out.numberofpoints != in.numberofpoints) {
			std::cout << "out.numberofpoints: " << out.numberofpoints << std::endl;
			std::cout << "points.size(): " << points.size() << std::endl;
			std::cout << "dumping out the points' coordinates" << std::endl;
			for (uint k = 0; k<points.size(); k++) {
				GVptr gv = points[k];
				std::cout << "  " << gv->coord << "  " << gv->idx << std::endl;
			}

			std::cout << "dumping out the segments" << std::endl;
			for (int k = 0; k<in.numberofsegments; k++)
				std::cout << "  " << in.segmentlist[k * 2 + 0]
				<< "; " << in.segmentlist[k * 2 + 1]
				<< " (" << in.segmentmarkerlist[k]
				<< ") " << std::endl;

			std::cout << "dumping out the solved for triangles now..." << std::endl;
			for (int k = 0; k<out.numberoftriangles; k++) {
				std::cout << "  "
					<< out.trianglelist[(k * 3) + 0] << "; "
					<< out.trianglelist[(k * 3) + 1] << "; "
					<< out.trianglelist[(k * 3) + 2] << std::endl;
			}
		}
		ENSURE(out.numberofpoints == in.numberofpoints);

		gtris.resize(out.numberoftriangles);
		for (int k = 0; k<out.numberoftriangles; k++) {
			GVptr       gv0 = points[out.trianglelist[(k * 3) + 0]];
			GVptr       gv1 = points[out.trianglelist[(k * 3) + 1]];
			GVptr       gv2 = points[out.trianglelist[(k * 3) + 2]];
			gtris[k] = iprob->newGenericTri(gv0, gv1, gv2);
		}

		// clean up after triangulate...
		// in free
		free(in.pointlist);
		free(in.pointmarkerlist);
		free(in.segmentlist);
		free(in.segmentmarkerlist);
		// out free
		free(out.pointlist);
		free(out.pointmarkerlist);
		free(out.trianglelist);
		free(out.segmentlist);
		free(out.segmentmarkerlist);
	}

private:
	void subdivideEdge(IsctProblem *iprob, GEptr ge, ShortVec<GEptr, 8> &edges)
	{
		if (ge->interior.size() == 0) {
			edges.push_back(ge);
		}
		else if (ge->interior.size() == 1) { // common case
			SEptr se0 = iprob->newSplitEdge(ge->ends[0], ge->interior[0], ge->boundary);
			SEptr se1 = iprob->newSplitEdge(ge->interior[0], ge->ends[1], ge->boundary);

			edges.push_back(se0);
			edges.push_back(se1);

			// get rid of old edge
			iprob->releaseEdge(ge);
		}
		else
		{ // sorting is the uncommon case determine the primary dimension and direction of the edge
			Vec3d       dir = ge->ends[1]->coord - ge->ends[0]->coord;
			uint        dim = (fabs(dir.x) > fabs(dir.y)) ? ((fabs(dir.x) > fabs(dir.z)) ? 0 : 2) : ((fabs(dir.y) > fabs(dir.z)) ? 1 : 2);
			double      sign = (dir.v[dim] > 0.0) ? 1.0 : -1.0;

			// pack the interior vertices into a vector for sorting
			std::vector< std::pair<double, IVptr> > verts;
			// for (IVptr iv : ge->interior)
			for (uint j = 0; j < ge->interior.size(); j++)
			{
				IVptr iv = ge->interior[j];
				// if the sort is ascending, then we're good...
				verts.push_back(std::make_pair(sign * iv->coord.v[dim], iv));
			}
			// ... and sort the vector
			std::sort(verts.begin(), verts.end());
			// then, write the verts into a new container with the endpoints
			std::vector<GVptr>  allv(verts.size() + 2);
			allv[0] = ge->ends[0];
			allv[allv.size() - 1] = ge->ends[1];
			for (uint k = 0; k<verts.size(); k++)
				allv[k + 1] = verts[k].second;

			// now create and accumulate new split edges
			for (uint i = 1; i<allv.size(); i++)
			{
				SEptr se = iprob->newSplitEdge(allv[i - 1], allv[i], ge->boundary);
				edges.push_back(se);
			}
			// get rid of old edge
			iprob->releaseEdge(ge);
		}
	}

public: // data
	ShortVec<IVptr, 4>      iverts;
	ShortVec<IEptr, 2>      iedges;
	// original triangle elements
	OVptr                   overts[3];
	OEptr                   oedges[3];

	ShortVec<GTptr, 8>      gtris;

	Tptr                    the_tri;
};


template<class VertData, class TriData>
class Mesh<VertData, TriData>::IsctProblem : public TopoCache
{
public:
	IsctProblem(Mesh *owner) : TopoCache(owner)
	{
		// initialize all the triangles to NOT have an associated tprob
		// TopoCache::tris.for_each([](Tptr t) {			// 21
		//	t->data = NULL;
		//});

		for (IterPool<TopoTri>::Block* block = TopoCache::tris.block_list; block != NULL; block = block->next)
		{
			//func((T*)(block));
			((TopoTri*)block)->data = NULL;
		}

		// Callibrate the quantization unit...
		double maxMag = 0.0;
		// for (VertData &v : TopoCache::mesh->verts)
		for (uint i = 0; i < TopoCache::mesh->verts.size(); i++)
		{
			maxMag = std::max(maxMag, max(abs(TopoCache::mesh->verts[i].pos)));
		}
		Quantization::callibrate(maxMag);

		// and use vertex auxiliary data to store quantized vertex coordinates
		uint N = TopoCache::mesh->verts.size();
		quantized_coords.resize(N);
		uint write = 0;
		//TopoCache::verts.for_each(	[&](Vptr v)				// 22
		//							{
		//							#ifdef _WIN32
		//								Vec3d raw = mesh->verts[v->ref].pos;
		//							#else
		//								Vec3d raw = TopoCache::mesh->verts[v->ref].pos;
		//							#endif
		//								quantized_coords[write].x = Quantization::quantize(raw.x);
		//								quantized_coords[write].y = Quantization::quantize(raw.y);
		//								quantized_coords[write].z = Quantization::quantize(raw.z);
		//								v->data = &(quantized_coords[write]);
		//								write++;
		//							}
		//);

		for (IterPool<TopoVert>::Block *block = TopoCache::verts.block_list; block != NULL; block = block->next)
		{
			//func((T*)(block));
#ifdef _WIN32
			Vec3d raw = mesh->verts[((Vptr)block)->ref].pos;
#else
			Vec3d raw = TopoCache::mesh->verts[((Vptr)block)->ref].pos;
#endif
			quantized_coords[write].x = Quantization::quantize(raw.x);
			quantized_coords[write].y = Quantization::quantize(raw.y);
			quantized_coords[write].z = Quantization::quantize(raw.z);
			((Vptr)block)->data = &(quantized_coords[write]);
			write++;
		}

	}

	virtual ~IsctProblem() {}

	// access auxiliary quantized coordinates
	inline Vec3d vPos(Vptr v) const {
		return *(reinterpret_cast<Vec3d*>(v->data));
	}

	Tprob getTprob(Tptr t) {
		Tprob prob = reinterpret_cast<Tprob>(t->data);
		if (!prob) {
			t->data = prob = tprobs.alloc();
			prob->init(this, t);
		}
		return prob;
	}
	GluePt newGluePt() {
		GluePt glue = glue_pts.alloc();
		glue->split_type = false;
		return glue;
	}

	inline IVptr newIsctVert(Eptr e, Tptr t, GluePt glue) {
		IVptr       iv = ivpool.alloc();
		iv->concrete = NULL;
		iv->coord = computeCoords(e, t);
		iv->glue_marker = glue;
		glue->copies.push_back(iv);
		return      iv;
	}
	inline IVptr newIsctVert(Tptr t0, Tptr t1, Tptr t2, GluePt glue) {
		IVptr       iv = ivpool.alloc();
		iv->concrete = NULL;
		iv->coord = computeCoords(t0, t1, t2);
		iv->glue_marker = glue;
		glue->copies.push_back(iv);
		return      iv;
	}
	inline IVptr newSplitIsctVert(Vec3d coords, GluePt glue) {
		IVptr       iv = ivpool.alloc();
		iv->concrete = NULL;
		iv->coord = coords;
		iv->glue_marker = glue;
		glue->copies.push_back(iv);
		return      iv;
	}
	inline IVptr copyIsctVert(IVptr orig) {
		IVptr       iv = ivpool.alloc();
		iv->concrete = NULL;
		iv->coord = orig->coord;
		iv->glue_marker = orig->glue_marker;
		orig->glue_marker->copies.push_back(iv);
		return      iv;
	}
	inline IEptr newIsctEdge(IVptr endpoint, Tptr tri_key) {
		IEptr       ie = iepool.alloc();
		ie->concrete = NULL;
		ie->boundary = false;
		ie->ends[0] = endpoint;
		endpoint->edges.push_back(ie);
		ie->ends[1] = NULL; // other end null
		ie->other_tri_key = tri_key;
		return      ie;
	}

	inline OVptr newOrigVert(Vptr v) {
		OVptr       ov = ovpool.alloc();
		ov->concrete = v;
		ov->coord = vPos(v);
		ov->boundary = true;
		return      ov;
	}
	inline OEptr newOrigEdge(Eptr e, OVptr v0, OVptr v1) {
		OEptr       oe = oepool.alloc();
		oe->concrete = e;
		oe->boundary = true;
		oe->ends[0] = v0;
		oe->ends[1] = v1;
		v0->edges.push_back(oe);
		v1->edges.push_back(oe);
		return      oe;
	}
	inline SEptr newSplitEdge(GVptr v0, GVptr v1, bool boundary) {
		SEptr       se = sepool.alloc();
		se->concrete = NULL;
		se->boundary = boundary;
		se->ends[0] = v0;
		se->ends[1] = v1;
		v0->edges.push_back(se);
		v1->edges.push_back(se);
		return      se;
	}

	inline GTptr newGenericTri(GVptr v0, GVptr v1, GVptr v2) {
		GTptr       gt = gtpool.alloc();
		gt->verts[0] = v0;
		gt->verts[1] = v1;
		gt->verts[2] = v2;
		gt->concrete = NULL;
		return      gt;
	}

	inline void releaseEdge(GEptr ge) {
		disconnectGE(ge);
		IEptr       ie = dynamic_cast<IEptr>(ge);
		if (ie) {
			iepool.free(ie);
		}
		else {
			OEptr   oe = dynamic_cast<OEptr>(ge);
			ENSURE(oe);
			oepool.free(oe);
		}
	}

	inline void killIsctVert(IVptr iv) 
	{
		iv->glue_marker->copies.erase(iv);
		if (iv->glue_marker->copies.size() == 0)
			glue_pts.free(iv->glue_marker);

		// for (GEptr ge : iv->edges)
		for (uint i = 0; i < iv->edges.size(); i++)
		{
			GEptr ge = iv->edges[i];
			// disconnect
			ge->interior.erase(iv);
			if (ge->ends[0] == iv)   ge->ends[0] = NULL;
			if (ge->ends[1] == iv)   ge->ends[1] = NULL;
		}

		ivpool.free(iv);
	}

	inline void killIsctEdge(IEptr ie) {
		// an endpoint may be an original vertex
		if (ie->ends[1])
			ie->ends[1]->edges.erase(ie);
		iepool.free(ie);
	}

	inline void killOrigVert(OVptr ov) {
		ovpool.free(ov);
	}

	inline void killOrigEdge(OEptr oe) {
		oepool.free(oe);
	}

	bool hasIntersections(); // test for iscts, exit if one is found
	void findIntersections();
	void resolveAllIntersections();
private:
	// if we encounter ambiguous degeneracies, then this routine returns false, indicating that the computation aborted.
	bool tryToFindIntersections();
	// In that case, we can perturb the positions of points
	void perturbPositions();
	// in order to give things another try, discard partial work
	void reset();
public:

	void dumpIsctPoints(std::vector<Vec3d> *points);
	void dumpIsctEdges(std::vector< std::pair<Vec3d, Vec3d> > *edges);

protected: // DATA
	IterPool<GluePointMarker>   glue_pts;
	IterPool<TriangleProblem>   tprobs;

	IterPool<IsctVertType>      ivpool;
	IterPool<OrigVertType>      ovpool;
	IterPool<IsctEdgeType>      iepool;
	IterPool<OrigEdgeType>      oepool;
	IterPool<SplitEdgeType>     sepool;
	IterPool<GenericTriType>    gtpool;
private:
	std::vector<Vec3d>          quantized_coords;
private:

	// 源码中共有两种调用
	inline void bvh_edge_tri_1(const int& _degeneracy_count)
	{
		std::vector< GeomBlob<Eptr> > edge_geoms;

		for (IterPool<TopoEdge>::Block *block = TopoCache::edges.block_list; block != NULL; block = block->next)
		{
			// func((T*)(block));
			edge_geoms.push_back(edge_blob((Eptr)block));	// 有时候需要强制类型转换，有时候不需要
		}

		AABVH<Eptr> edgeBVH(edge_geoms);

		// use the acceleration structure
		bool aborted = false;

		for (IterPool<TopoTri>::Block *block = TopoCache::tris.block_list; block != NULL; block = block->next)
		{
			// compute BBox
			BBox3d bbox = buildBox(((Tptr)(block)));
			if (!aborted)
			{
				// do a recursive search and invoke the action at each piece of geometry
				std::stack< AABVHNode<Eptr>* >  nodes;
				nodes.push(edgeBVH.root);

				while (!nodes.empty())
				{
					AABVHNode<Eptr> *node = nodes.top();
					nodes.pop();

					// check bounding box isct
					if (!hasIsct(node->bbox, bbox))  continue;

					// otherwise...
					if (node->isLeaf())
					{
						// for (uint bid : node->blobids)
						for (uint ii = 0; ii < node->blobids.size(); ii++)
						{
							uint bid = node->blobids[ii];
							if (hasIsct(bbox, edgeBVH.blobs[bid].bbox))
							{
								// action(edgeBVH.blobs[bid].id);
								if (checkIsct(edgeBVH.blobs[bid].id, ((Tptr)(block))))
								{
									GluePt glue = newGluePt();
									glue->edge_tri_type = true;
									glue->e = edgeBVH.blobs[bid].id;
									glue->t[0] = ((Tptr)(block));
									// first add point and edges to the pierced triangle
									IVptr iv = getTprob(((Tptr)(block)))->addInteriorEndpoint(this, edgeBVH.blobs[bid].id, glue);
									// for (Tptr tri : edgeBVH.blobs[bid].id->tris)
									for (uint jj = 0; jj < edgeBVH.blobs[bid].id->tris.size(); jj++)
									{
										getTprob(edgeBVH.blobs[bid].id->tris[jj])->addBoundaryEndpoint(this, ((Tptr)(block)), edgeBVH.blobs[bid].id, iv);
									}
								}
								if (_degeneracy_count > 0)
									aborted = true;
							}
						}
					}
					else {
						nodes.push(node->left);
						nodes.push(node->right);
					}
				}
			}
		}
	}

	inline void bvh_edge_tri_2(bool& foundIsct, const int& degeneracy_count)
	{
		std::vector< GeomBlob<Eptr> > edge_geoms;

		for (IterPool<TopoEdge>::Block *block = TopoCache::edges.block_list; block != NULL; block = block->next)
		{
			// func((T*)(block));
			edge_geoms.push_back(edge_blob((Eptr)block));	// 有时候需要强制类型转换，有时候不需要
		}

		AABVH<Eptr> edgeBVH(edge_geoms);

		// use the acceleration structure
		bool aborted = false;

		for (IterPool<TopoTri>::Block *block = TopoCache::tris.block_list; block != NULL; block = block->next)
		{
			// compute BBox
			BBox3d bbox = buildBox(((Tptr)(block)));
			if (!aborted)
			{
				// do a recursive search and invoke the action at each piece of geometry
				std::stack< AABVHNode<Eptr>* >  nodes;
				nodes.push(edgeBVH.root);

				while (!nodes.empty())
				{
					AABVHNode<Eptr> *node = nodes.top();
					nodes.pop();

					// check bounding box isct
					if (!hasIsct(node->bbox, bbox))  continue;

					// otherwise...
					if (node->isLeaf())
					{
						// for (uint bid : node->blobids)
						for (uint ii = 0; ii < node->blobids.size(); ii++)
						{
							uint bid = node->blobids[ii];
							if (hasIsct(bbox, edgeBVH.blobs[bid].bbox))
							{
								// action(blobs[bid].id);
								if (checkIsct(edgeBVH.blobs[bid].id, ((Tptr)(block))))
								{
									foundIsct = true;
									aborted = true;
								}
								if (Empty3d::degeneracy_count > 0)
									aborted = true;
							}
						}
					}
					else {
						nodes.push(node->left);
						nodes.push(node->right);
					}
				}
			}
		}
	}

	inline GeomBlob<Eptr> edge_blob(Eptr e)
	{
		GeomBlob<Eptr>  blob;
		blob.bbox = buildBox(e);
		blob.point = (blob.bbox.minp + blob.bbox.maxp) / 2.0;
		blob.id = e;
		return blob;
	}

	inline BBox3d buildBox(Eptr e) const
	{
		Vec3d p0 = vPos(e->verts[0]);
		Vec3d p1 = vPos(e->verts[1]);
		return BBox3d(min(p0, p1), max(p0, p1));
	}

	inline BBox3d buildBox(Tptr t) const
	{
		Vec3d p0 = vPos(t->verts[0]);
		Vec3d p1 = vPos(t->verts[1]);
		Vec3d p2 = vPos(t->verts[2]);
		return BBox3d(min(p0, min(p1, p2)), max(p0, max(p1, p2)));
	}

	inline void marshallArithmeticInput(Empty3d::TriIn &input, Tptr t) const
	{
		input.p[0] = vPos(t->verts[0]);
		input.p[1] = vPos(t->verts[1]);
		input.p[2] = vPos(t->verts[2]);
	}

	inline void marshallArithmeticInput(Empty3d::EdgeIn &input, Eptr e) const
	{
		input.p[0] = vPos(e->verts[0]);
		input.p[1] = vPos(e->verts[1]);
	}

	inline void marshallArithmeticInput(Empty3d::TriEdgeIn &input, Eptr e, Tptr t) const
	{
		marshallArithmeticInput(input.edge, e);
		marshallArithmeticInput(input.tri, t);
	}

	inline void marshallArithmeticInput(Empty3d::TriTriTriIn &input, Tptr t0, Tptr t1, Tptr t2) const
	{
		marshallArithmeticInput(input.tri[0], t0);
		marshallArithmeticInput(input.tri[1], t1);
		marshallArithmeticInput(input.tri[2], t2);
	}

	bool checkIsct(Eptr e, Tptr t) const;
	bool checkIsct(Tptr t0, Tptr t1, Tptr t2) const;

	Vec3d computeCoords(Eptr e, Tptr t) const;
	Vec3d computeCoords(Tptr t0, Tptr t1, Tptr t2) const;

	void fillOutVertData(GluePt glue, VertData &data);
	void fillOutTriData(Tptr tri, Tptr parent);
private:
	class EdgeCache;

private: // functions here to get around a GCC bug...
	void createRealPtFromGluePt(GluePt glue);
	void createRealTriangles(Tprob tprob, EdgeCache &ecache);
};

struct TriTripleTemp
{
	Tptr t0, t1, t2;
	TriTripleTemp(Tptr tp0, Tptr tp1, Tptr tp2) : t0(tp0), t1(tp1), t2(tp2) {}
};

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::IsctProblem::tryToFindIntersections()
{
	Empty3d::degeneracy_count = 0;
	// Find all edge-triangle intersection points
	bvh_edge_tri_1(Empty3d::degeneracy_count);

	if (Empty3d::degeneracy_count > 0) {
		return false;   // restart / abort
	}

	// we're going to peek into the triangle problems in order to identify potential candidates for Tri-Tri-Tri intersections
	std::vector<TriTripleTemp> triples;
	//tprobs.for_each([&](Tprob tprob)				// 25
	//				{
	//					Tptr t0 = tprob->the_tri;
	//					// Scan pairs of existing edges to create candidate triples
	//					for (uint i = 0; i < tprob->iedges.size(); i++)
	//					{
	//						for (uint j = i + 1; j < tprob->iedges.size(); j++)
	//						{
	//							//func(vec[i], vec[j]);
	//							IEptr& ie1 = tprob->iedges[i];
	//							IEptr& ie2 = tprob->iedges[j];

	//							Tptr t1 = ie1->other_tri_key;
	//							Tptr t2 = ie2->other_tri_key;
	//							// This triple might be considered three times, one for each triangle it contains.
	//							// To prevent duplication, only proceed if this is the least triangle according to an arbitrary ordering
	//							if (t0 < t1 && t0 < t2)
	//							{
	//								Mesh<VertData, TriData>::Tprob prob1 = reinterpret_cast<Mesh<VertData, TriData>::Tprob>(t1->data);
	//								for (IEptr ie : prob1->iedges) {
	//									if (ie->other_tri_key == t2) {
	//										// ADD THE TRIPLE
	//										triples.push_back(TriTripleTemp(t0, t1, t2));
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//);

	for (IterPool<TriangleProblem>::Block *block = tprobs.block_list; block != NULL; block = block->next)
	{
		//  func((T*)(block));
		Tptr t0 = ((Tprob)block)->the_tri;
		// Scan pairs of existing edges to create candidate triples

		for (uint i = 0; i < ((Tprob)block)->iedges.size(); i++)
		{
			for (uint j = i + 1; j < ((Tprob)block)->iedges.size(); j++)
			{
				//func(vec[i], vec[j]);
				IEptr& ie1 = ((Tprob)block)->iedges[i];
				IEptr& ie2 = ((Tprob)block)->iedges[j];

				Tptr t1 = ie1->other_tri_key;
				Tptr t2 = ie2->other_tri_key;
				// This triple might be considered three times, one for each triangle it contains.
				// To prevent duplication, only proceed if this is the least triangle according to an arbitrary ordering
				if (t0 < t1 && t0 < t2)
				{
					Mesh<VertData, TriData>::Tprob prob1 = reinterpret_cast<Mesh<VertData, TriData>::Tprob>(t1->data);
					// for (IEptr ie : prob1->iedges)
					for (uint k = 0; k < prob1->iedges.size(); k++)
					{
						if (prob1->iedges[k]->other_tri_key == t2) {
							// ADD THE TRIPLE
							triples.push_back(TriTripleTemp(t0, t1, t2));
						}
					}
				}
			}
		}
	}


	// Now, we've collected a list of Tri-Tri-Tri intersection candidates.
	// Check to see if the intersections actually exist.
	// for (TriTripleTemp t : triples)
	for (uint i = 0; i < triples.size(); i++)
	{
		TriTripleTemp t = triples[i];
		if (!checkIsct(t.t0, t.t1, t.t2))    continue;

		// Abort if we encounter a degeneracy
		if (Empty3d::degeneracy_count > 0)   break;

		GluePt glue = newGluePt();
		glue->edge_tri_type = false;
		glue->t[0] = t.t0;
		glue->t[1] = t.t1;
		glue->t[2] = t.t2;
		getTprob(t.t0)->addInteriorPoint(this, t.t1, t.t2, glue);
		getTprob(t.t1)->addInteriorPoint(this, t.t0, t.t2, glue);
		getTprob(t.t2)->addInteriorPoint(this, t.t0, t.t1, glue);
	}
	if (Empty3d::degeneracy_count > 0) {
		return false;   // restart / abort
	}
	return true;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::perturbPositions()
{
	const double EPSILON = 1.0e-5; // perturbation epsilon
	// for (Vec3d &coord : quantized_coords)
	for (uint i = 0; i < quantized_coords.size(); i++)
	{
		Vec3d perturbation(Quantization::quantize(drand(-EPSILON, EPSILON)),
			Quantization::quantize(drand(-EPSILON, EPSILON)),
			Quantization::quantize(drand(-EPSILON, EPSILON)));
		quantized_coords[i] += perturbation;
	}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::reset()
{
	// the data pointer in the triangles points to tproblems that we're about to destroy, so zero out all those pointers first!
	//tprobs.for_each([](Tprob tprob) {				// 26
	//	Tptr t = tprob->the_tri;
	//	t->data = NULL;
	//});

	for (IterPool<TriangleProblem>::Block *block = tprobs.block_list; block != NULL; block = block->next)
	{
		// 用 ((Tprob)(block)) 替换 tprob
		Tptr t = ((Tprob)(block))->the_tri;
		t->data = NULL;
	}

	glue_pts.clear();
	tprobs.clear();

	ivpool.clear();
	ovpool.clear();
	iepool.clear();
	oepool.clear();
	sepool.clear();
	gtpool.clear();
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::findIntersections()
{
	int nTrys = 5;
	perturbPositions(); // always perturb for safety...
	while (nTrys > 0) {
		if (!tryToFindIntersections()) {
			reset();
			perturbPositions();
			nTrys--;
		}
		else
			break;
	}
	if (nTrys <= 0) {
		CORK_ERROR("Ran out of tries to perturb the mesh");
		exit(1);
	}
	// ok all points put together, all triangle problems assembled.
	// Some intersection edges may have original vertices as endpoints we consolidate the problems to check for cases like these.
	//tprobs.for_each([&](Tprob tprob)				// 27
	//				{
	//					tprob->consolidate(this);
	//				}
	//);

	for (IterPool<TriangleProblem>::Block *block = tprobs.block_list; block != NULL; block = block->next)
	{
		// 用 ((Tprob)(block)) 替换 tprob
		((Tprob)(block))->consolidate(this);
	}
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::IsctProblem::hasIntersections()
{
	bool foundIsct = false;
	Empty3d::degeneracy_count = 0;
	// Find some edge-triangle intersection point...
	//bvh_edge_tri_2(foundIsct, Empty3d::degeneracy_count);

	if (Empty3d::degeneracy_count > 0 || foundIsct)
	{
		std::cout << "This self-intersection might be spurious. Degeneracies were detected." << std::endl;
		return true;
	}
	else
		return false;
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::IsctProblem::checkIsct(Eptr e, Tptr t) const
{
	// simple bounding box cull; for acceleration, not correctness
	BBox3d      ebox = buildBox(e);
	BBox3d      tbox = buildBox(t);
	if (!hasIsct(ebox, tbox))
		return      false;

	// must check whether the edge and triangle share a vertex
	// if so, then trivially we know they intersect in exactly that vertex so we discard this case from consideration.
	if (hasCommonVert(e, t))
		return false;

	Empty3d::TriEdgeIn input;
	marshallArithmeticInput(input, e, t);
	bool empty = Empty3d::emptyExact(input);
	return !empty;
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::IsctProblem::checkIsct(Tptr t0, Tptr t1, Tptr t2) const
{
	// This f_unction should only be called if we've already
	// identified that the intersection edges  (t0,t1), (t0,t2), (t1,t2)  exist.
	// From this, we can conclude that each pair of triangles shares no more than a single vertex in common.
	// If each of these shared vertices is different from each other, then we could legitimately have a triple intersection point,
	// but if all three pairs share the same vertex in common, then the intersection of the three triangles must be that vertex.
	// So, we must check for such a single vertex in common amongst the three triangles
	Vptr common = commonVert(t0, t1);
	if (common) {
		for (uint i = 0; i<3; i++)
			if (common == t2->verts[i])
				return false;
	}
	Empty3d::TriTriTriIn input;
	marshallArithmeticInput(input, t0, t1, t2);
	bool empty = Empty3d::emptyExact(input);
	return !empty;
}

template<class VertData, class TriData>
Vec3d Mesh<VertData, TriData>::IsctProblem::computeCoords(Eptr e, Tptr t) const
{
	Empty3d::TriEdgeIn input;
	marshallArithmeticInput(input, e, t);
	Vec3d coords = Empty3d::coordsExact(input);
	return coords;
}

template<class VertData, class TriData>
Vec3d Mesh<VertData, TriData>::IsctProblem::computeCoords(Tptr t0, Tptr t1, Tptr t2) const
{
	Empty3d::TriTriTriIn input;
	marshallArithmeticInput(input, t0, t1, t2);
	Vec3d coords = Empty3d::coordsExact(input);
	return coords;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::fillOutVertData(GluePt glue, VertData &data)
{
	if (glue->split_type) { // manually inserted split point
		uint v0i = glue->e->verts[0]->ref;
		uint v1i = glue->e->verts[1]->ref;
		data.isctInterpolate(TopoCache::mesh->verts[v0i],
			TopoCache::mesh->verts[v1i]);
	}
	else
		if (glue->edge_tri_type)
		{
			IsctVertEdgeTriInput<VertData, TriData>      input;
			for (uint k = 0; k<2; k++) {
				uint    vid = glue->e->verts[k]->ref;
				input.e[k] = &(TopoCache::mesh->verts[vid]);
			}
			for (uint k = 0; k<3; k++) {
				uint    vid = glue->t[0]->verts[k]->ref;
				input.t[k] = &(TopoCache::mesh->verts[vid]);
			}
			data.isct(input);
		}
		else {
			IsctVertTriTriTriInput<VertData, TriData>    input;
			for (uint i = 0; i<3; i++) {
				for (uint j = 0; j<3; j++) {
					uint    vid = glue->t[i]->verts[j]->ref;
					input.t[i][j] = &(TopoCache::mesh->verts[vid]);
				}
			}
			data.isct(input);
		}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::fillOutTriData(Tptr piece, Tptr parent)
{
	TopoCache::mesh->subdivide_tri(piece->ref, parent->ref);
}

template<class VertData, class TriData>
class Mesh<VertData, TriData>::IsctProblem::EdgeCache
{
public:
	EdgeCache(IsctProblem *ip) : iprob(ip), edges(ip->mesh->verts.size()) {}

	Eptr operator()(Vptr v0, Vptr v1) {
		uint i = v0->ref;
		uint j = v1->ref;
		if (i > j) std::swap(i, j);

		uint N = edges[i].size();
		for (uint k = 0; k<N; k++)
			if (edges[i][k].vid == j)
				return edges[i][k].e;
		// if not existing, create it
		edges[i].push_back(EdgeEntry(j));
		Eptr e = edges[i][N].e = iprob->newEdge();
		e->verts[0] = v0;
		e->verts[1] = v1;
		v0->edges.push_back(e);
		v1->edges.push_back(e);

		return e;
	}

	// k = 0, 1, or 2
	Eptr getTriangleEdge(GTptr gt, uint k, Tptr big_tri)
	{
		GVptr   gv0 = gt->verts[(k + 1) % 3];
		GVptr   gv1 = gt->verts[(k + 2) % 3];
		Vptr    v0 = gv0->concrete;
		Vptr    v1 = gv1->concrete;
		// if neither of these are intersection points, then this is a pre-existing edge...
		Eptr    e = NULL;
		if (typeid(gv0) == typeid(OVptr) &&
			typeid(gv1) == typeid(OVptr)
			) {
			// search through edges of original triangle...
			for (uint c = 0; c<3; c++) {
				Vptr corner0 = big_tri->verts[(c + 1) % 3];
				Vptr corner1 = big_tri->verts[(c + 2) % 3];
				if ((corner0 == v0 && corner1 == v1) || (corner0 == v1 && corner1 == v0))
				{
					e = big_tri->edges[c];
				}
			}
			ENSURE(e); // Yell if we didn't find an edge
		}
		// otherwise, we need to check the cache to find this edge
		else
		{
			e = operator()(v0, v1);
		}
		return e;
	}

	Eptr maybeEdge(GEptr ge)
	{
		uint i = ge->ends[0]->concrete->ref;
		uint j = ge->ends[1]->concrete->ref;
		if (i > j) std::swap(i, j);

		uint N = edges[i].size();
		for (uint k = 0; k<N; k++)
			if (edges[i][k].vid == j)
				return edges[i][k].e;
		// if we can't find it
		return NULL;
	}

private:
	struct EdgeEntry {
		EdgeEntry(uint id) : vid(id) {}
		EdgeEntry() {}
		uint vid;
		Eptr e;
	};

	IsctProblem *iprob;
	std::vector< ShortVec<EdgeEntry, 8> >   edges;
};

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::createRealPtFromGluePt(GluePt glue) {
	ENSURE(glue->copies.size() > 0);
	Vptr        v = TopoCache::newVert();
	VertData    &data = TopoCache::mesh->verts[v->ref];
	data.pos = glue->copies[0]->coord;
	fillOutVertData(glue, data);
	// for (IVptr iv : glue->copies)
	for (uint i = 0; i < glue->copies.size(); i++)
		glue->copies[i]->concrete = v;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::createRealTriangles(Tprob tprob, EdgeCache &ecache)
{
	// for (GTptr gt : tprob->gtris)
	for (uint i = 0; i < tprob->gtris.size(); i++)
	{
		GTptr gt = tprob->gtris[i];
		Tptr t = TopoCache::newTri();
		gt->concrete = t;
		Tri &tri = TopoCache::mesh->tris[t->ref];
		for (uint k = 0; k<3; k++)
		{
			Vptr v = gt->verts[k]->concrete;
			t->verts[k] = v;
			v->tris.push_back(t);
			tri.v[k] = v->ref;

			Eptr e = ecache.getTriangleEdge(gt, k, tprob->the_tri);
			e->tris.push_back(t);
			t->edges[k] = e;
		}
		fillOutTriData(t, tprob->the_tri);
	}
	// Once all the pieces are hooked up, let's kill the old triangle!
	TopoCache::deleteTri(tprob->the_tri);
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::resolveAllIntersections()
{
	// solve a subdivision problem in each triangle
	//tprobs.for_each([&](Tprob tprob)					// 28
	//				{
	//					tprob->subdivide(this);
	//				}
	//);

	for (IterPool<TriangleProblem>::Block *block = tprobs.block_list; block != NULL; block = block->next)
	{
		// 用 ((Tprob)(block)) 替换 tprob
		((Tprob)(block))->subdivide(this);
	}

	// now we have diced up triangles inside each triangle problem

	// Let's go through the glue points and create a new concrete vertex object for each of these.
	//glue_pts.for_each([&](GluePt glue) {				// 29
	//	createRealPtFromGluePt(glue);
	//});

	for (IterPool<GluePointMarker>::Block *block = glue_pts.block_list; block != NULL; block = block->next)
	{
		//func((T*)(block));
		// ((GluePt)(block)) 替换 glue
		createRealPtFromGluePt(((GluePt)(block)));
	}

	EdgeCache ecache(this);

	// Now that we have concrete vertices plugged in, we can go through the diced triangle pieces and create concrete triangles for each of those.
	// Along the way, let's go ahead and hook up edges as appropriate
	//tprobs.for_each([&](Tprob tprob) {					// 30
	//	createRealTriangles(tprob, ecache);
	//});

	for (IterPool<TriangleProblem>::Block *block = tprobs.block_list; block != NULL; block = block->next)
	{
		createRealTriangles(((Tprob)(block)), ecache);
	}

	// mark all edges as normal by zero-ing out the data pointer
	//TopoCache::edges.for_each([](Eptr e) {				// 31
	//	e->data = 0;
	//});

	for (IterPool<TopoEdge>::Block *block = TopoCache::edges.block_list; block != NULL; block = block->next)
	{
		((Eptr)block)->data = 0;
	}

	// then iterate over the edges formed by intersections
	// (i.e. those edges without the boundary flag set in each triangle) and mark those by setting the data pointer
	//iepool.for_each([&](IEptr ie)							// 32
	//				{
	//					// every ie must be non-boundary
	//					Eptr e = ecache.maybeEdge(ie);
	//					ENSURE(e);
	//					e->data = (void*)1;
	//				}
	//);

	for (IterPool<IsctEdgeType>::Block *block = iepool.block_list; block != NULL; block = block->next)
	{
		// every ie must be non-boundary
		Eptr e = ecache.maybeEdge(((IEptr)(block)));
		ENSURE(e);
		e->data = (void*)1;
	}

	//sepool.for_each([&](SEptr se)							// 33
	//				{
	//					Eptr e = ecache.maybeEdge(se);
	//					ENSURE(e);
	//					e->data = (void*)1;
	//				}
	//);

	for (IterPool<SplitEdgeType>::Block *block = sepool.block_list; block != NULL; block = block->next)
	{
		Eptr e = ecache.maybeEdge(((SEptr)(block)));
		ENSURE(e);
		e->data = (void*)1;
	}

	// This basically takes care of everything EXCEPT one detail
	// *) The base mesh data structures still need to be compacted

	// This detail should be handled by the calling code...
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::dumpIsctPoints(std::vector<Vec3d> *points)
{
	points->resize(glue_pts.size());
	uint write = 0;
	//glue_pts.for_each(	[&](GluePt glue)				// 34
	//					{
	//						ENSURE(glue->copies.size() > 0);
	//						IVptr       iv = glue->copies[0];
	//						(*points)[write] = iv->coord;
	//						write++;
	//					}
	//);

	for (IterPool<GluePointMarker>::Block *block = glue_pts.block_list; block != NULL; block = block->next)
	{
		ENSURE(((GluePt)(block))->copies.size() > 0);
		IVptr iv = ((GluePt)(block))->copies[0];
		(*points)[write] = iv->coord;
		write++;
	}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::IsctProblem::dumpIsctEdges(std::vector< std::pair<Vec3d, Vec3d> > *edges)
{
	edges->clear();
	//tprobs.for_each([&](Tprob tprob)							// 35
	//				{
	//					for (IEptr ie : tprob->iedges) {
	//						GVptr gv0 = ie->ends[0];
	//						GVptr gv1 = ie->ends[1];
	//						edges->push_back(std::make_pair(gv0->coord, gv1->coord));
	//					}
	//				}
	//);

	for (IterPool<TriangleProblem>::Block *block = tprobs.block_list; block != NULL; block = block->next)
	{
		// for (IEptr ie : ((Tprob)(block))->iedges)
		for (uint i = 0; i < ((Tprob)(block))->iedges.size(); i++)
		{
			IEptr ie = ((Tprob)(block))->iedges[i];
			GVptr gv0 = ie->ends[0];
			GVptr gv1 = ie->ends[1];
			edges->push_back(std::make_pair(gv0->coord, gv1->coord));
		}
	}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::testingComputeStaticIsctPoints(std::vector<Vec3d> *points)
{
	IsctProblem iproblem(this);
	iproblem.findIntersections();
	iproblem.dumpIsctPoints(points);
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::testingComputeStaticIsct(std::vector<Vec3d> *points, std::vector< std::pair<Vec3d, Vec3d> > *edges)
{
	IsctProblem iproblem(this);

	iproblem.findIntersections();
	iproblem.dumpIsctPoints(points);
	iproblem.dumpIsctEdges(edges);
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::resolveIntersections()
{
	IsctProblem iproblem(this);

	iproblem.findIntersections();
	iproblem.resolveAllIntersections();
	iproblem.commit();
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::isSelfIntersecting()
{
	IsctProblem iproblem(this);
	return iproblem.hasIntersections();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////// mesh.bool.tpp
#include <queue>

template<class VertData, class TriData>
class Mesh<VertData, TriData>::BoolProblem
{
public:
	BoolProblem(Mesh *owner) : mesh(owner) {}
	virtual ~BoolProblem() {}

	// do things
	void doSetup(Mesh &rhs);

	// choose what to remove
	enum TriCode { KEEP_TRI, DELETE_TRI, FLIP_TRI };

	inline void doDeleteAndFlip_1()
	{
		TopoCache topocache(mesh);
		std::vector<Tptr> toDelete;

		for (IterPool<TopoTri>::Block *block = topocache.tris.block_list; block != NULL; block = block->next)
		{
			byte data = boolData(((Tptr)(block))->ref);

			TriCode code;
			if ((data & 2) == 2)     // part of op 0/1 INSIDE op 1/0
				code = BoolProblem::DELETE_TRI;
			else
				code = BoolProblem::KEEP_TRI;
			switch (code)
			{
			case DELETE_TRI:
				toDelete.push_back(((Tptr)(block)));
				break;
			case FLIP_TRI:
				topocache.flipTri(((Tptr)(block)));
				break;
			case KEEP_TRI:
			default:
				break;
			}
		}
		for (uint i = 0; i < toDelete.size(); i++)
		{
			topocache.deleteTri(toDelete[i]);
		}
		topocache.commit();
	}

	inline void doDeleteAndFlip_2()
	{
		TopoCache topocache(mesh);
		std::vector<Tptr> toDelete;

		for (IterPool<TopoTri>::Block *block = topocache.tris.block_list; block != NULL; block = block->next)
		{
			byte data = boolData(((Tptr)(block))->ref);

			TriCode code;

			if (data == 2 ||         // part of op 0 INSIDE op 1
				data == 1)           // part of op 1 OUTSIDE op 0
				code = BoolProblem::DELETE_TRI;
			else if (data == 3)      // part of op 1 INSIDE op 1
				code = BoolProblem::FLIP_TRI;
			else                    // part of op 0 OUTSIDE op 1
				code = BoolProblem::KEEP_TRI;

			switch (code)
			{
			case DELETE_TRI:
				toDelete.push_back(((Tptr)(block)));
				break;
			case FLIP_TRI:
				topocache.flipTri(((Tptr)(block)));
				break;
			case KEEP_TRI:
			default:
				break;
			}
		}
		// for (Tptr tptr : toDelete)
		for (uint i = 0; i < toDelete.size(); i++)
		{
			topocache.deleteTri(toDelete[i]);
		}
		topocache.commit();
	}

	inline void doDeleteAndFlip_3()
	{
		TopoCache topocache(mesh);
		std::vector<Tptr> toDelete;

		for (IterPool<TopoTri>::Block *block = topocache.tris.block_list; block != NULL; block = block->next)
		{
			byte data = boolData(((Tptr)(block))->ref);

			TriCode code;
			if ((data & 2) == 0)     // part of op 0/1 OUTSIDE op 1/0
				code = BoolProblem::DELETE_TRI;
			else                    // part of op 0/1 INSIDE op 1/0
				code = BoolProblem::KEEP_TRI;
			switch (code)
			{
			case DELETE_TRI:
				toDelete.push_back(((Tptr)(block)));
				break;
			case FLIP_TRI:
				topocache.flipTri(((Tptr)(block)));
				break;
			case KEEP_TRI:
			default:
				break;
			}
		}
		// for (Tptr tptr : toDelete)
		for (uint i = 0; i < toDelete.size(); i++)
		{
			topocache.deleteTri(toDelete[i]);
		}
		topocache.commit();
	}

	inline void doDeleteAndFlip_4()
	{
		TopoCache topocache(mesh);
		std::vector<Tptr> toDelete;

		for (IterPool<TopoTri>::Block *block = topocache.tris.block_list; block != NULL; block = block->next)
		{
			byte data = boolData(((Tptr)(block))->ref);

			TriCode code;
			if ((data & 2) == 0)     // part of op 0/1 OUTSIDE op 1/0
				code = BoolProblem::KEEP_TRI;
			else                    // part of op 0/1 INSIDE op 1/0
				code = BoolProblem::FLIP_TRI;
			switch (code)
			{
			case DELETE_TRI:
				toDelete.push_back(((Tptr)(block)));
				break;
			case FLIP_TRI:
				topocache.flipTri(((Tptr)(block)));
				break;
			case KEEP_TRI:
			default:
				break;
			}
		}
		//for (Tptr tptr : toDelete) {
		for (uint i = 0; i < toDelete.size(); i++)
		{
			topocache.deleteTri(toDelete[i]);
		}
		topocache.commit();
	}

private: // methods
	struct BoolEdata {
		bool is_isct;
	};

	inline byte& boolData(uint tri_id) {
		return mesh->tris[tri_id].data.bool_alg_data;
	}

	void populateECache()
	{
		ecache = mesh->createEGraphCache<BoolEdata>();

		// label some of the edges as intersection edges and others as not
		//ecache.for_each( [&](uint i, uint j, EGraphEntry<BoolEdata> &entry)			// 36
		//				{
		//					entry.data.is_isct = false;
		//					byte operand = boolData(entry.tids[0]);
		//					for (uint k = 1; k < entry.tids.size(); k++)
		//					{
		//						if (boolData(entry.tids[k]) != operand)
		//						{
		//							entry.data.is_isct = true;
		//							break;
		//						}
		//					}
		//				}
		//);

		for (uint i = 0; i < ecache.skeleton.size(); i++)
		{
			for (uint jj = 0; jj < ecache.skeleton[i].size(); jj++)
			{
				EGraphEntry<BoolEdata>& entry = ecache.skeleton[i][jj];
				//action(i, entry.vid, entry);
				entry.data.is_isct = false;
				byte operand = boolData(entry.tids[0]);
				for (uint k = 1; k < entry.tids.size(); k++)
				{
					if (boolData(entry.tids[k]) != operand)
					{
						entry.data.is_isct = true;
						break;
					}
				}
			}
		}

	}

	bool isInside(uint tid, byte operand) {
		// find the point to trace outward from...
		Vec3d p(0, 0, 0);
		p += mesh->verts[mesh->tris[tid].a].pos;
		p += mesh->verts[mesh->tris[tid].b].pos;
		p += mesh->verts[mesh->tris[tid].c].pos;
		p /= 3.0;
		// ok, we've got the point, now let's pick a direction
		Ray3d r;
		r.p = p;
		r.r = Vec3d(drand(0.5, 1.5), drand(0.5, 1.5), drand(0.5, 1.5));

		int winding = 0;
		// pass all triangles over ray
		// for (Tri &tri : mesh->tris)
		for (uint i = 0; i < mesh->tris.size(); i++)
		{
			Tri &tri = mesh->tris[i];
			// ignore triangles from the same operand surface
			if ((tri.data.bool_alg_data & 1) == operand)   continue;

			double flip = 1.0;
			uint   a = tri.a;
			uint   b = tri.b;
			uint   c = tri.c;
			Vec3d va = mesh->verts[a].pos;
			Vec3d vb = mesh->verts[b].pos;
			Vec3d vc = mesh->verts[c].pos;
			// normalize vertex order (to prevent leaks)
			if (a > b) { std::swap(a, b); std::swap(va, vb); flip = -flip; }
			if (b > c) { std::swap(b, c); std::swap(vb, vc); flip = -flip; }
			if (a > b) { std::swap(a, b); std::swap(va, vb); flip = -flip; }

			double t;
			Vec3d bary;
			if (isct_ray_triangle(r, va, vb, vc, &t, &bary)) {
				Vec3d normal = flip * cross(vb - va, vc - va);
				if (dot(normal, r.r) > 0.0) { // UNSAFE
					winding++;
				}
				else {
					winding--;
				}
			}
		}

		// now, we've got a winding number to work with...
		return winding > 0;
	}

private: // data
	Mesh                        *mesh;
	EGraphCache<BoolEdata>      ecache;
};


static inline double triArea(Vec3d a, Vec3d b, Vec3d c)
{
	return len(cross(b - a, c - a));
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::BoolProblem::doSetup(Mesh &rhs)
{
	// Label surfaces...
	//mesh->for_tris( [](TriData &tri, VertData&, VertData&, VertData&)
	//				{
	//					tri.bool_alg_data = 0;
	//				}
	//);

	mesh->mesh_for_tris();

	//rhs.for_tris( [](TriData &tri, VertData&, VertData&, VertData&){	tri.bool_alg_data = 1;} );

	rhs.rhs_for_tris();

	mesh->disjointUnion(rhs);
	mesh->resolveIntersections();

	populateECache();

	// form connected components;
	// we get one component for each connected component in one of the two input meshes.
	// These components are not necessarily uniformly inside or outside of the other operand mesh.
	UnionFind uf(mesh->tris.size());
	//for_ecache([&](uint, uint, bool, const ShortVec<uint, 2> &tids) {
	//	uint tid0 = tids[0];
	//	for (uint k = 1; k < tids.size(); k++)
	//		uf.unionIds(tid0, tids[k]);
	//});

	//ecache.for_each([&](uint i, uint j, EGraphEntry<BoolEdata> &entry)		// 37
	//				{
	//					if (entry.data.is_isct)
	//					{
	//						ShortVec<uint, 2> tid0s;
	//						ShortVec<uint, 2> tid1s;
	//						for (uint tid : entry.tids)
	//						{
	//							if (boolData(tid) & 1)
	//								tid1s.push_back(tid);
	//							else
	//								tid0s.push_back(tid);
	//						}
	//						uint tid1 = tid1s[0];
	//						for (uint k = 1; k < tid1s.size(); k++)
	//							uf.unionIds(tid1, tid1s[k]);

	//						uint tid0 = tid0s[0];
	//						for (uint k = 1; k < tid0s.size(); k++)
	//							uf.unionIds(tid0, tid0s[k]);
	//					}
	//					else
	//					{
	//						uint tid0 = entry.tids[0];
	//						for (uint k = 1; k < entry.tids.size(); k++)
	//							uf.unionIds(tid0, entry.tids[k]);
	//					}
	//				}
	//);

	for (uint i = 0; i < ecache.skeleton.size(); i++)
	{
		for (uint j = 0; j < ecache.skeleton[i].size(); j++)
		{
			EGraphEntry<BoolEdata> &entry = ecache.skeleton[i][j];
			// action(i, entry.vid, entry);
			if (entry.data.is_isct)
			{
				ShortVec<uint, 2> tid0s;
				ShortVec<uint, 2> tid1s;
				// for (uint tid : entry.tids)
				for (uint k = 0; k < entry.tids.size(); k++)
				{
					uint tid = entry.tids[k];
					if (boolData(tid) & 1)
						tid1s.push_back(tid);
					else
						tid0s.push_back(tid);
				}
				uint tid1 = tid1s[0];
				for (uint k = 1; k < tid1s.size(); k++)
					uf.unionIds(tid1, tid1s[k]);

				uint tid0 = tid0s[0];
				for (uint k = 1; k < tid0s.size(); k++)
					uf.unionIds(tid0, tid0s[k]);
			}
			else
			{
				uint tid0 = entry.tids[0];
				for (uint k = 1; k < entry.tids.size(); k++)
					uf.unionIds(tid0, entry.tids[k]);
			}
		}
	}

	// we re-organize the results of the union find as follows:
	std::vector<uint> uq_ids(mesh->tris.size(), uint(-1));
	std::vector< std::vector<uint> > components;
	for (uint i = 0; i < mesh->tris.size(); i++) {
		uint ufid = uf.find(i);
		if (uq_ids[ufid] == uint(-1)) { // unassigned
			uint N = components.size();
			components.push_back(std::vector<uint>());

			uq_ids[ufid] = uq_ids[i] = N;
			components[N].push_back(i);
		}
		else { // assigned already
			uq_ids[i] = uq_ids[ufid]; // propagate assignment
			components[uq_ids[i]].push_back(i);
		}
	}

	std::vector<bool> visited(mesh->tris.size(), false);

	// find the "best" triangle in each component, and ray cast to determine inside-ness vs. outside-ness
	for (uint i = 0; i < components.size(); i++)
	{
		std::vector<uint> &comp = components[i];
		// find max according to score
		uint best_tid = comp[0];
		double best_area = 0.0;
		// SEARCH
		// for (uint tid : comp)
		for (uint j = 0; j < comp.size(); j++)
		{
			uint tid = comp[j];
			Vec3d va = mesh->verts[mesh->tris[tid].a].pos;
			Vec3d vb = mesh->verts[mesh->tris[tid].b].pos;
			Vec3d vc = mesh->verts[mesh->tris[tid].c].pos;

			double area = triArea(va, vb, vc);
			if (area > best_area) {
				best_area = area;
				best_tid = tid;
			}
		}

		byte operand = boolData(best_tid);
		bool inside = isInside(best_tid, operand);

		// NOW PROPAGATE classification throughout the component. do a breadth first propagation
		std::queue<uint> work;

		// begin by tagging the first triangle
		boolData(best_tid) |= (inside) ? 2 : 0;
		visited[best_tid] = true;
		work.push(best_tid);

		while (!work.empty())
		{
			uint curr_tid = work.front();
			work.pop();

			for (uint k = 0; k < 3; k++)
			{
				uint a = mesh->tris[curr_tid].v[k];
				uint b = mesh->tris[curr_tid].v[(k + 1) % 3];
				EGraphEntry<BoolEdata> &entry = ecache(a, b);
				byte inside_sig = boolData(curr_tid) & 2;
				if (entry.data.is_isct)  inside_sig ^= 2;
				// for (uint tid : entry.tids)
				for (uint kk = 0; kk < entry.tids.size(); kk++)
				{
					uint tid = entry.tids[kk];
					if (visited[tid])
						continue;
					if ((boolData(tid) & 1) != operand)
						continue;

					boolData(tid) |= inside_sig;
					visited[tid] = true;
					work.push(tid);
				}
			}
		}
	}
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::boolUnion(Mesh &rhs)
{
	BoolProblem bprob(this);
	bprob.doSetup(rhs);
	bprob.doDeleteAndFlip_1();
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::boolDiff(Mesh &rhs)
{
	BoolProblem bprob(this);
	bprob.doSetup(rhs);
	bprob.doDeleteAndFlip_2();
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::boolIsct(Mesh &rhs)
{
	BoolProblem bprob(this);
	bprob.doSetup(rhs);
	bprob.doDeleteAndFlip_3();
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::boolXor(Mesh &rhs)
{
	BoolProblem bprob(this);
	bprob.doSetup(rhs);
	bprob.doDeleteAndFlip_4();
}