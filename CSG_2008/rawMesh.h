#pragma once


#include "vec.h"
#include <vector>

// vertex and triangle data must satisfy the following minimal specifications!
// These exact field names must be used!
struct MinimalVertexData
{
	Vec3d pos; // required for both RawMesh and Mesh
};

struct MinimalTriangleData
{
	int a, b, c; // Vertex Ids: (only used in raw mesh form)
};

// Raw mesh presents an exposed interface to a mesh.
// This allows for easy input/output of data to and from the more powerful Mesh data structure, as well as supporting more basic mesh applications
template<class VertData, class TriData>
struct RawMesh
{
	std::vector<VertData> vertices;
	std::vector<TriData> triangles;
};