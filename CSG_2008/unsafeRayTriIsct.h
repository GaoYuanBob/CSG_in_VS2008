
#pragma once

#include "ray.h"
#include "vec.h"

struct isct_ray_triangle_result
{
	double t;
	Vec3d bary;
};

inline
bool isct_ray_triangle(
	Ray3d               ray,
	Vec3d               va,
	Vec3d               vb,
	Vec3d               vc,
	double              *t,
	Vec3d               *bary
	) {
	ENSURE(t);
	ENSURE(bary);

	// re-center the problem at the base point of the ray
	va -= ray.p;
	vb -= ray.p;
	vc -= ray.p;

	// Then compute volumes of tetrahedra spanning
	//  * the base point / ray direction line segment
	//  * an edge of the triangle
	// Keeping orientations in mind...
	double volAB = det(va, vb, ray.r);
	double volBC = det(vb, vc, ray.r);
	double volCA = -det(va, vc, ray.r);
	// then also compute the volume of tet with the entire triangle as a face...
	double volABC = det(va, vb, vc);
	// if any of the signs of the edge tests
	// disagree with the sign of the whole triangle, then
	// the ray does not pass through the triangle
	if (volAB * volABC < 0 ||
		volBC * volABC < 0 ||
		volCA * volABC < 0)      return false;

	// otherwise, compute the t - value for the ray to intersect
	// if this is negative, then the client can detect that the
	// ray would have to travel backwards to hit the triangle in question.
	double edgeSum = volAB + volBC + volCA;
	if (edgeSum == 0)            return false;
	*t = volABC / (volAB + volBC + volCA);
	if (*t <= 0)                 return false;

	*bary = Vec3d(volBC / edgeSum, volCA / edgeSum, volAB / edgeSum);
	return true;
}
