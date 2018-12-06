
#pragma once

#ifndef uint
typedef unsigned int uint;
#endif

// if a mesh is taken as input, the client must manage the memory
// if a mesh is given as output, please use the provided function to free the allocated memory.
struct CorkTriMesh
{
	uint    n_triangles;
	uint    n_vertices;
	uint    *triangles;
	float   *vertices;
};

void freeCorkTriMesh(CorkTriMesh *mesh);

// the inputs to Boolean operations must be "solid":
//  -   closed (aka. watertight; see comment at bottom)
//  -   non-self-intersecting
// additionally, inputs should use a counter-clockwise convention
// for triangle facing.  If the triangles are presented in clockwise orientation, the object is interpreted as its unbounded complement

// This function will test whether or not a mesh is solid
bool isSolid(CorkTriMesh mesh);

// Boolean operations follow
void computeUnion(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);			// result = A U B

void computeDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);		// result = A - B

void computeIntersection(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);	// result = A ^ B

void computeSymmetricDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);	// result = A XOR B

																						// Not a Boolean operation, but related:
																						//  No portion of either surface is deleted.  However, the curve of intersection between the two surfaces is made explicit,
																						//  such that the two surfaces are now connected.
void resolveIntersections(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

