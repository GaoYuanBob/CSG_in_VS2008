#pragma once

#include <string>

//#include "triMesh.h"

#include "rawMesh.h"

/*
*  Files provides a wrapper for different file types and a common
*  data view for the rest of the program.  This wrapper was introduced
*  to make it easier to support multiple file types using other people's file importer/exporter code
*/

namespace Files {

	// all functions with integer return values here are intended to return an error count as if they were a main f_unction

	struct FileVertex : public MinimalVertexData {};
	struct FileTriangle : public MinimalTriangleData {};
	typedef RawMesh<FileVertex, FileTriangle> FileMesh;
	//using FileMesh = RawMesh<FileVertex, FileTriangle>;

	// generic filetype functions
	// these detect which filetype to use by inspecting the filename

	int readOFF(std::string filename, FileMesh *mesh);
	int writeOFF(std::string filename, FileMesh *mesh);

} // end namespace Files
