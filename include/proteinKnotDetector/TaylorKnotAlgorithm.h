/*
 * Name        : Protein Knot Detector
 * Author      : Brad Lee
 * Version     : 1.00
 * License     : GNU LGPL v3
 * Description : Detects knots
 *
 * Works Cited:
 * Taylor, W. A deeply knotted protein structure and how it might fold.
 * Nature 406, 916–919 (2000) doi:10.1038/35022623
 */

#ifndef PKD_TAYLOR_KNOT_ALGORITHM_H
#define PKD_TAYLOR_KNOT_ALGORITHM_H

// c++17
#include <string>
#include <functional>
#include <iostream>
#include <optional>
#include <memory>
#include <limits>

// c
#include <stdio.h>
#include <string.h>

#include "amalgamated.h"

/*
 * PKD = Protein Knot Detector
 */
namespace PKD {

// Copyright (C) 2016 by Doug Baldwin.
// This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International
// License (http://creativecommons.org/licenses/by-sa/4.0/).
// High and low bounds on t values that are considered to represent a ray intersecting
// a triangle's plane.
#define tFar   10000.0f
#define tNear  0.0000001f
constexpr double lowest_double = std::numeric_limits<double>::lowest();
/* CROSS, DOT, and SUB3 Macros for 3-component vectors
 * used in original Moeller and Trumbore algorithm.
 *
 * Citation for macros from Moeller and Trumbore's original code
 * Created on: 1997
 * 		Author: Tomas Moeller -- Chalmers University of Technology
 * 		& Ben Trumbore -- Cornell University
 * 		Title:  Fast, Minimum Storage Ray /Triangle Intersection
 * 		Source: http://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
 */
#define CROSS(res, v1, v2)\
		res[0]=v1[1]*v2[2]-v1[2]*v2[1];\
		res[1]=v1[2]*v2[0]-v1[0]*v2[2];\
		res[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB3(res,v1,v2)\
		res[0]=v1[0]-v2[0];\
		res[1]=v1[1]-v2[1];\
		res[2]=v1[2]-v2[2];
inline bool MollerTrumboreRayTriangleIntersection(double *v0, double *v1,
		double *v2, double *rayOrigin, double *rayDirection) {
#ifdef TAYLOR_SMOOTH_DEBUG_INTERSECT
				if (i < TAYLOR_SMOOTH_DEBUG_DEPTH
						&& k < TAYLOR_SMOOTH_DEBUG_DEPTH) {
					printf(
							"TRI{(%.2f,%.2f,%.2f);(%.2f,%.2f,%.2f);(%.2f,%.2f,%.2f)} line{(%.2f,%.2f,%.2f);(%.2f,%.2f,%.2f)} TEST\n",
							i, k, v0[0], v0[1], v0[2], v1[0], v1[1], v1[2],
							v2[0], v2[1], v2[2], rayOrigin[0], rayOrigin[1],
							rayOrigin[2], rayDirection[0], rayDirection[1],
							rayDirection[2]);
				}
#endif
	double v0v1[3]; // Find vectors for two edges sharing vertex 0
	SUB3(v0v1, v1, v0);
	double v0v2[3];
	SUB3(v0v2, v2, v0);
	double pvec[3]; // Begin calculating determinant;
	CROSS(pvec, rayDirection, v0v2); // also used to calculate U parameter
	const double det = DOT(v0v1, pvec); // If determinant is near zero, ray lies in plane of triangle
	if (det > -1 * lowest_double && det < lowest_double) // Ray is parallel to this triangle.
		return false;
	const double inv_det = 1.0f / det;
	double tvec[3]; // Calculate vector from vertex to ray origin
	SUB3(tvec, rayOrigin, v0);
	const double u = DOT(tvec, pvec) * inv_det; // Calculate U parameter and test bounds
	if (u < 0.0f || u > 1.0f)
		return false;
	double qvec[3]; // Prepare to test V parameter
	CROSS(qvec, tvec, v0v1);
	const double v = DOT(rayDirection, qvec) * inv_det; // Calculate V parameter and test bounds
	if (v < 0.0f || u + v >= 1.0f)
		return false;
	//const double t = DOT(edge2, qvec) * inv_det; // Calculate t, final check to see if ray intersects triangle. Test to
	//if (t <= lowest_double || t >= 1-lowest_double) // see if t > tFar added for consistency with other algorithms in experiment.
	//	return false;
	// intersection found, don't move vertex
#ifdef TAYLOR_SMOOTH_DEBUG_INTERSECT
				if (i < TAYLOR_SMOOTH_DEBUG_DEPTH
						&& k < TAYLOR_SMOOTH_DEBUG_DEPTH) {
					printf(
							"TRI{(%.2f,%.2f,%.2f);(%.2f,%.2f,%.2f);(%.2f,%.2f,%.2f)} line{(%.2f,%.2f,%.2f);(%.2f,%.2f,%.2f)} INTERSECTION\n",
							i, k, v0[0], v0[1], v0[2], v1[0], v1[1], v1[2],
							v2[0], v2[1], v2[2], rayOrigin[0], rayOrigin[1],
							rayOrigin[2], rayDirection[0], rayDirection[1],
							rayDirection[2]);
				}
#endif
	return true;
}

/*
 * William R. Taylor Knot Detection Algorithm
 */
class TaylorKnotAlgorithm {
private:
	std::unique_ptr<DoubleMatrix> m;
public:
	std::unique_ptr<DoubleMatrix> getMatrix();
	void setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr);
	void smooth(unsigned int nRepeat, bool isSmoothStraightenCollinear);
	void smoothAuto();
	/* After about 50 iterations of smoothing,
	 * the knot now may be detected.
	 */
};

std::unique_ptr<DoubleMatrix> TaylorKnotAlgorithm::getMatrix() {
	return std::move(m);
}
void TaylorKnotAlgorithm::setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr) {
	m = std::move(matrixPtr);
}

/*
 * CRITICAL ALGORITHM
 * Function calls would create too much overhead so we use #define
 * for maximum computational efficiency
 */
void TaylorKnotAlgorithm::smooth(unsigned int nRepeat = 1,
		bool isSmoothStraightenCollinear = true) {
	double *x = m->m; // x is an alias for the vertex matrix
	double *v0, *v1, *v2, *v0a, *v1a, *v2a, *rayOrigin, *rayDirection;
	double v1p[3];
	/* v# are the operated vertexes
	 * v#a are the committed vertexes
	 * v#p are the prime vertexes (vertex after move)
	 */
	int n = m->n - 3;
	for (unsigned int j = 0; j < nRepeat; j++) {
		for (int i = 3; i < n; i += 3) {
			v0a = x + i - 3;
			v1a = x + i;
			v2a = x + i + 3;
			v1p[0] = ((v0a[0] + v2a[0]) / 2 + v1a[0]) / 2;
			v1p[1] = ((v0a[1] + v2a[1]) / 2 + v1a[1]) / 2;
			v1p[2] = ((v0a[2] + v2a[2]) / 2 + v1a[2]) / 2;
#ifdef TAYLOR_SMOOTH_DEBUG
			if (i < TAYLOR_SMOOTH_DEBUG_DEPTH) {
				printf(
						"i#%d i-1:(%.2f,%.2f,%.2f) i:(%.2f,%.2f,%.2f) i+1:(%.2f,%.2f,%.2f) i':(%.2f,%.2f,%.2f)\n",
						i, v0a[0], v0a[1], v0a[2], v1a[0], v1a[1], v1a[2],
						v2a[0], v2a[1], v2a[2], v1p[0], v1p[1], v1p[2]);
			}
#endif
			/* check that the triangles {i'-1,i,i'} and {i;i';i+1}
			 * did not intersect any line segment {j'-1;j'}(j<i) before the move point
			 * or any line {j;j+1}(j>i) following.
			 * implemented with the Möller–Trumbore intersection algorithm
			 */
			// triangle {i'-1,i,i'} and line {j'-1;j'}(j<i)
			v0 = v0a;
			v1 = v1a;
			v2 = v1p;
			for (int k = 3; k < i; k += 3) {
				rayOrigin = x + k - 3; // setup ray
				rayDirection = x + k;
#ifdef TAYLOR_SMOOTH_DEBUG
				if (i < TAYLOR_SMOOTH_DEBUG_DEPTH && k < TAYLOR_SMOOTH_DEBUG_DEPTH) {
					printf(
							"i#%d k#%d triangle {i'-1,i,i'} and line {j'-1;j'}(j<i)\n",
							i, k);
				}
#endif
				if (MollerTrumboreRayTriangleIntersection(v0, v1, v2, rayOrigin,
						rayDirection)) {
					goto intersect;
				} else {
					continue;
				}
			}
			// triangle {i;i';i+1} and line {j'-1;j'}(j<i)
			v0 = v1a;
			v1 = v1p;
			v2 = v2a;
			for (int k = 3; k < i; k += 3) {
				rayOrigin = x + k - 3;
				rayDirection = x + k;
#ifdef TAYLOR_SMOOTH_DEBUG
				if (i < TAYLOR_SMOOTH_DEBUG_DEPTH && k < TAYLOR_SMOOTH_DEBUG_DEPTH) {
					printf(
							"i#%d k#%d triangle {i;i';i+1} and line {j'-1;j'}(j<i)\n",
							i, k);
				}
#endif
				if (MollerTrumboreRayTriangleIntersection(v0, v1, v2, rayOrigin,
						rayDirection)) {
					goto intersect;
				} else {
					continue;
				}
			}
			// triangle {i'-1,i,i'} and line {j;j+1}(j>i)
			v0 = v0a;
			v1 = v1a;
			v2 = v1p;
			for (int k = i + 3; k < n; k += 3) {
				rayOrigin = x + k; // setup ray
				rayDirection = x + k + 3;
#ifdef TAYLOR_SMOOTH_DEBUG
				if (i < TAYLOR_SMOOTH_DEBUG_DEPTH && k < TAYLOR_SMOOTH_DEBUG_DEPTH) {
					printf(
							"i#%d k#%d triangle {i'-1,i,i'} and line {j;j+1}(j>i)\n",
							i, k);
				}
#endif
				if (MollerTrumboreRayTriangleIntersection(v0, v1, v2, rayOrigin,
						rayDirection)) {
					goto intersect;
				} else {
					continue;
				}
			}
			// triangle {i;i';i+1} and line {j;j+1}(j>i)
			v0 = v1a;
			v1 = v1p;
			v2 = v2a;
			for (int k = i + 3; k < n; k += 3) {
				rayOrigin = x + k;
				rayDirection = x + k + 3;
#ifdef TAYLOR_SMOOTH_DEBUG
				if (i < TAYLOR_SMOOTH_DEBUG_DEPTH && k < TAYLOR_SMOOTH_DEBUG_DEPTH) {
					printf(
							"i#%d k#%d triangle {i;i';i+1} and line {j;j+1}(j>i)\n",
							i, k);
				}
#endif
				if (MollerTrumboreRayTriangleIntersection(v0, v1, v2, rayOrigin,
						rayDirection)) {
					goto intersect;
				} else {
					continue;
				}
			}
			// both triangles don't intersect, commit vertex move
			v1a[0] = v1p[0];
			v1a[1] = v1p[1];
			v1a[2] = v1p[2];
			/* for algorithm efficiency we use goto instead of if else statements
			 * if and else statements and unnecessary comparisons
			 */
#ifdef TAYLOR_SMOOTH_DEBUG_INTERSECT
			goto nointersect;
#endif
			intersect: ;
#ifdef TAYLOR_SMOOTH_DEBUG_INTERSECT
			printf("i#%d INTERSECTION\n", i);
#endif
#ifdef TAYLOR_SMOOTH_DEBUG
			if (i < TAYLOR_SMOOTH_DEBUG_DEPTH) {
				printf("i#%d\n ", i);
			}
#endif
#ifdef TAYLOR_SMOOTH_DEBUG_INTERSECT
			nointersect: ;
#endif
		}
	}
}

void TaylorKnotAlgorithm::smoothAuto() {

}

} // namespace PKD

#endif
