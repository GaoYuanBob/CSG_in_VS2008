
// GY
// This file contains a command line program that can be used to exercise Cork's functionality without having to write any code.

#include "files.h"
#include "cork.h"

#include <iostream>
#include <sstream>
#include <assert.h>

using namespace std;
using namespace Files;

#ifndef uint
typedef unsigned int uint;
#endif

void file2corktrimesh(const Files::FileMesh &in, CorkTriMesh *out)
{
	out->n_vertices = in.vertices.size();
	out->n_triangles = in.triangles.size();

	out->triangles = new uint[(out->n_triangles) * 3];
	out->vertices = new float[(out->n_vertices) * 3];

	for (uint i = 0; i < out->n_triangles; i++) {
		(out->triangles)[3 * i + 0] = in.triangles[i].a;
		(out->triangles)[3 * i + 1] = in.triangles[i].b;
		(out->triangles)[3 * i + 2] = in.triangles[i].c;
	}

	for (uint i = 0; i < out->n_vertices; i++) {
		(out->vertices)[3 * i + 0] = in.vertices[i].pos.x;
		(out->vertices)[3 * i + 1] = in.vertices[i].pos.y;
		(out->vertices)[3 * i + 2] = in.vertices[i].pos.z;
	}
}

void corktrimesh2file(CorkTriMesh in, Files::FileMesh &out)
{
	out.vertices.resize(in.n_vertices);
	out.triangles.resize(in.n_triangles);

	for (uint i = 0; i < in.n_triangles; i++) {
		out.triangles[i].a = in.triangles[3 * i + 0];
		out.triangles[i].b = in.triangles[3 * i + 1];
		out.triangles[i].c = in.triangles[3 * i + 2];
	}

	for (uint i = 0; i < in.n_vertices; i++) {
		out.vertices[i].pos.x = in.vertices[3 * i + 0];
		out.vertices[i].pos.y = in.vertices[3 * i + 1];
		out.vertices[i].pos.z = in.vertices[3 * i + 2];
	}
}

int main()
{
	FileMesh fm_res;
	FileMesh *fm1 = new FileMesh(), *fm2 = new FileMesh();

	CorkTriMesh *ctm1 = new CorkTriMesh(), *ctm2 = new CorkTriMesh(), *ctm_res = new CorkTriMesh();	// 不写这个new，生成的是临时变量？

	// int readOFF(string filename, FileMesh *data)
	readOFF("C:/Users/Administrator/Desktop/U_CSG_A.off", fm1);
	readOFF("C:/Users/Administrator/Desktop/U_CSG_B.off", fm2);

	// void file2corktrimesh(const Files::FileMesh &in, CorkTriMesh *out)
	file2corktrimesh(*fm1, ctm1);
	file2corktrimesh(*fm2, ctm2);

	assert(fm1->triangles.size() != 0 && fm2->triangles.size() != 0);

	// isSolid(CorkTriMesh mesh);
	//if (!isSolid(*ctm1))
	//{
	//	printf("第一个模型不是solid模型!\n");
	//	system("pause");
	//}
	//if (!isSolid(*ctm2))
	//{
	//	printf("第二个模型不是solid模型!\n");
	//	system("pause");
	//}

	// void computeUnion(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out)
	computeUnion(*ctm1, *ctm2, ctm_res);

	// void corktrimesh2file(CorkTriMesh in, Files::FileMesh &out)
	corktrimesh2file(*ctm_res, fm_res);

	// int writeOFF(string filename, FileMesh *data)
	writeOFF("C:/Users/Administrator/Desktop/U.off", &fm_res);

	printf("Complete!!!\n");
	system("pause");

	return 0;
}