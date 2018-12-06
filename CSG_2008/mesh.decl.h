#pragma once

#include "rawMesh.h"
#include "prelude.h"
#include "vec.h"
#include "ray.h"
#include "shortVec.h"
#include "iterPool.h"

#include <vector>
#include <set>
#include <functional>

struct BoolVertexData {};
struct BoolTriangleData {
	byte bool_alg_data; // internal use by algorithm, please copy value when the triangle is subdivided
};

template<class VertData, class TriData>
struct IsctVertEdgeTriInput
{
	VertData*   e[2];
	VertData*   t[3];
};

template<class VertData, class TriData>
struct IsctVertTriTriTriInput
{
	VertData*   t[3][3];
};

template<class VertData, class TriData>
struct SubdivideTriInput
{
	TriData*    pt;
	VertData*   pv[3];
	VertData*   v[3];
};

// in order to perform remeshing, VertData and TriData must support
struct RemeshVertexData {
	bool manifold; // whether this point is manifold.
				   // useful for modifying interpolation behavior
				   // specify how to compute new data for a vertex in the event of either an edge collapse (via merge) or edge split (via interpolate)

};

struct RemeshOptions
{
	double maxEdgeLength;
	double minEdgeLength;
	double minAngle;
	double maxAngle;
	RemeshOptions() : maxEdgeLength(1.0), minEdgeLength(0.3), minAngle(5.0), maxAngle(170.0) {}
};

// only for internal use, please do not use as client
struct TopoVert;
struct TopoEdge;
struct TopoTri;

typedef TopoVert* Vptr;
typedef TopoEdge* Eptr;
typedef TopoTri*  Tptr;
// end internal items

template<class VertData, class TriData>
class Mesh
{
public:
	Mesh() {};
	Mesh(Mesh &src);
	Mesh(const RawMesh<VertData, TriData> &raw);
	virtual ~Mesh() {};

	void operator=(Mesh &src);

	// validity check:
	//  - all numbers are well-defined and finite
	//  - all triangle vertex indices are in the right range
	bool valid() const;

	RawMesh<VertData, TriData> raw() const;

	inline int numVerts() const { return verts.size(); }
	inline int numTris() const { return tris.size(); }

	inline void mesh_for_tris()
	{
		for (int i = 0; i < tris.size(); i++)
			tris[i].data.bool_alg_data = 0;
	}

	inline void rhs_for_tris()
	{
		for (int i = 0; i < tris.size(); i++)
			tris[i].data.bool_alg_data = 1;
	}

	// form the disjoint union of two meshes
	void disjointUnion(const Mesh &cp);

	struct Isct {
		Ray3d   ray;
		bool    exists;

		uint    tri_id;
		Vec3d   isct;
		Vec3d   bary;
	};
	Isct pick(Ray3d ray);

	bool isClosed(); // checks if the mesh is closed

public: // REMESHING module
		// REQUIRES:
		//  - MinimalData
		//  - RemeshData
	void remesh();
	RemeshOptions remesh_options;

public: // ISCT (intersections) module
	void resolveIntersections(); // makes all intersections explicit
	bool isSelfIntersecting(); // is the mesh self-intersecting?
							   // TESTING
	void testingComputeStaticIsctPoints(std::vector<Vec3d> *points);
	void testingComputeStaticIsct(std::vector<Vec3d> *points, std::vector< std::pair<Vec3d, Vec3d> > *edges);

public: // BOOLean operation module
		// all of the form
		//      this = this OP rhs
	void boolUnion(Mesh &rhs);
	void boolDiff(Mesh &rhs);
	void boolIsct(Mesh &rhs);
	void boolXor(Mesh &rhs);

private:    // Internal Formats
	struct Tri {
		TriData data;
		union {
			struct { uint a, b, c; };// vertex ids  
			uint v[3];
		};
		inline Tri() {}
	};

	inline void merge_tris(uint tid_result, uint tid0, uint tid1) {
		tris[tid_result].data.merge(tris[tid0].data, tris[tid1].data);
	}

	inline void split_tris(uint t0ref, uint t1ref, uint t_orig_ref) {
		TriData::split(tris[t0ref].data, tris[t1ref].data, tris[t_orig_ref].data);
	}

	inline void move_tri(Tri &t_new, Tri &t_old)
	{
		t_new.data.move(t_old.data);
	}

	inline void subdivide_tri(uint t_piece_ref, uint t_parent_ref)
	{
		SubdivideTriInput<VertData, TriData>     input;
		input.pt = &(tris[t_parent_ref].data);
		for (uint k = 0; k < 3; k++) {
			input.pv[k] = &(verts[tris[t_parent_ref].v[k]]);
			input.v[k] = &(verts[tris[t_piece_ref].v[k]]);
		}
		tris[t_piece_ref].data.subdivide(input);
	}

private:    // DATA
	std::vector<Tri>        tris;
	std::vector<VertData>   verts;

private:    // caches
	struct NeighborEntry {
		uint vid;
		ShortVec<uint, 2> tids;
		inline NeighborEntry() {}
		inline NeighborEntry(uint vid_) : vid(vid_) {}
	};
	struct NeighborCache {
		std::vector< ShortVec<NeighborEntry, 8> > skeleton;
		inline NeighborEntry& operator()(uint i, uint j) {
			uint N = skeleton[i].size();
			for (uint k = 0; k < N; k++) {
				if (skeleton[i][k].vid == j)
					return skeleton[i][k];
			}
			skeleton[i].push_back(NeighborEntry(j));
			return skeleton[i][N];
		}
	};
	NeighborCache createNeighborCache();

	// parallel to vertex array
	std::vector<uint> getComponentIds();

	// like the neighbor cache, but more customizable
	template<class Edata>
	struct EGraphEntry {
		uint                vid;
		ShortVec<uint, 2>   tids;
		Edata               data;
		inline EGraphEntry() {}
		inline EGraphEntry(uint vid_) : vid(vid_) {}
	};
	template<class Edata>
	struct EGraphCache {
		std::vector< ShortVec<EGraphEntry<Edata>, 8> > skeleton;
		inline EGraphEntry<Edata> & operator()(uint i, uint j) {
			uint N = skeleton[i].size();
			for (uint k = 0; k < N; k++) {
				if (skeleton[i][k].vid == j)
					return skeleton[i][k];
			}
			skeleton[i].push_back(EGraphEntry<Edata>(j));
			return skeleton[i][N];
		}
	};
	template<class Edata>
	EGraphCache<Edata> createEGraphCache();

private:    // TopoCache Support
	struct TopoCache;
private:    // Isct Support
	class  IsctProblem; // implements intersection functionality
	class TriangleProblem; // support type for IsctProblem
	typedef TriangleProblem* Tprob;
private:    // Bool Support
	class BoolProblem;

private:    // Remeshing Support
	struct RemeshScratchpad;

	Eptr allocateRemeshEdge(RemeshScratchpad &);
	void deallocateRemeshEdge(RemeshScratchpad &, Eptr);

	void edgeSplit(RemeshScratchpad &, Eptr e_split);
	void edgeCollapse(RemeshScratchpad &, Eptr e_collapse, bool collapsing_tetrahedra_disappear);

	// Need edge scoring routines...
	void scoreAndEnqueue(std::set< std::pair<double, Eptr> > &queue, Eptr edge);
	void dequeue(std::set< std::pair<double, Eptr> > &queue, Eptr edge);
	double computeEdgeScore(Eptr edge);

	// support functions
	void populateTriFromTopoTri(Tptr t);

	// calls the first function once, then the second once for each triangle
	inline void edgeNeighborhood_GY(Eptr edge, double& edge_length, double& min_angle, double& max_angle)
	{
		Vptr       v0 = edge->verts[0];
		Vptr       v1 = edge->verts[1];
		VertData   &data0 = verts[v0->ref];
		VertData   &data1 = verts[v1->ref];

		// 源代码中的 once(data0, data1);
		edge_length = len(data1.pos - data0.pos);

		// for (Tptr tri : edge->tris)
		for (uint ii = 0; ii < edge->tris.size(); ii++)
		{
			Tptr tri = edge->tris[ii];
			TriData &tdata = tris[tri->ref].data;
			for (uint i = 0; i < 3; i++)
			{
				if (tri->edges[i] == edge)
				{
					Vptr vopp = tri->verts[i];
					VertData &dataopp = verts[vopp->ref];

					// 源代码中的 each_tri(data0, data1, dataopp, tdata);
					Vec3d e0 = data0.pos - dataopp.pos;
					Vec3d e1 = data1.pos - dataopp.pos;
					double cos_angle = dot(e0, e1) / (len(e0)*len(e1));
					double angle = rad2deg(acos(cos_angle));
					if (angle < min_angle) min_angle = angle;
					if (angle > max_angle) max_angle = angle;
				}
			}
		}
	}
};
