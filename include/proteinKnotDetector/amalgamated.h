/*
 * @name Protein Knot Detector
 * @author Brad Lee
 * @version 1.00
 * @license GNU LGPL v3
 * @brief Core command line interface tools for detecting protein knots.
 * @details A library providing the core tools for detecting protein knots meant
 * for the command line interface.
 *
 * Works Cited:
 * Taylor, W. A deeply knotted protein structure and how it might fold.
 * Nature 406, 916–919 (2000) doi:10.1038/35022623
 */

#ifndef PKD_AMALGAMATED_H
#define PKD_AMALGAMATED_H

/*
 * Configuration
 */
#include <proteinKnotDetector/config.h>

// c++17
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <unordered_map>

// c
#include <stdio.h>
#include <string.h>

/*
 * PKD = Protein Knot Detector
 */
namespace PKD {

/*
 * our s x 3 matrix
 * s = amino acid chain length
 */
class DoubleMatrix {
public:
	/*
	 * 1D array is used instead of a 2D array to guarantee contiguous
	 * memory is used, so the algorithm will run faster.
	 */
	double *m;
	std::size_t n;
	std::size_t s;
	DoubleMatrix(std::size_t size) {
		s = size;
		n = s * 3;
		m = new double[n];
	}
	DoubleMatrix(const DoubleMatrix&) = delete;
	DoubleMatrix& operator=(const DoubleMatrix&) = delete;
	~DoubleMatrix() {
		delete[] m;
	}
	void printMatrix() {
		for (size_t i = 0; i < n; i += 3) {
			printf("%f %f %f\n", m[i], m[i + 1], m[i + 2]);
		}
	}
};

/*
 * Protein container
 */
class Protein {
public:
	std::unordered_map<std::string, std::unique_ptr<DoubleMatrix>> carbonAlphaMatrixMap;
	/* 0: Unknown
	 * 1: MMDB
	 * 2: PKA (This software's minimal management)
	 */
	int managerId;
};

} // namespace PKD

#endif
