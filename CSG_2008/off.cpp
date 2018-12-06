#include "files.h"

#include <iostream>
#include <fstream>
#include <vector>
using std::ifstream;
using std::ofstream;
using std::endl;

namespace Files {

	using std::string;
	using std::vector;

	int readOFF(string filename, FileMesh *data)
	{
		if (!data) return 1;

		ifstream in; // will close on exit from this read function
		in.open(filename.c_str());
		if (!in) return 1;

		// "OFF"
		string filetype;
		in >> filetype;
		if (filetype != "OFF") return 1;

		// counts of things
		int numvertices, numfaces, numedges;
		in >> numvertices >> numfaces >> numedges;
		if (!in) return 1;
		data->vertices.resize(numvertices);
		data->triangles.resize(numfaces);

		// vertex data
		for (uint i = 0; i < data->vertices.size(); i++)
		{
			Vec3d &p = data->vertices[i].pos;
			in >> p.x >> p.y >> p.z;
		}
		if (!in) return 1;

		// face data
		for (uint i = 0; i < data->triangles.size(); i++)
		{
			int polysize;
			in >> polysize;
			if (polysize != 3)   return 1;

			in >> data->triangles[i].a >> data->triangles[i].b >> data->triangles[i].c;
		}
		if (!in) return 1;

		return 0;
	}

	int writeOFF(string filename, FileMesh *data)
	{
		if (!data) return 1;

		ofstream out;
		out.open(filename.c_str());
		if (!out) return 1;

		// "OFF"
		out << "OFF" << endl;

		// numvertices, numfaces, numedges=0
		int numvertices = data->vertices.size();
		int numfaces = data->triangles.size();
		out << numvertices << ' ' << numfaces << ' ' << 0 << endl;

		// vertex data
		for (uint i = 0; i < data->vertices.size(); i++)
		{
			const Vec3d &p = data->vertices[i].pos;
			out << p.x << ' ' << p.y << ' ' << p.z << endl;
		}

		// face data
		for (uint i = 0; i < data->triangles.size(); i++)
		{
			out << "3 " << data->triangles[i].a << ' ' << data->triangles[i].b << ' ' << data->triangles[i].c << endl;
		}
		if (!out) return 1;

		return 0;
	}

} // end namespace Files
