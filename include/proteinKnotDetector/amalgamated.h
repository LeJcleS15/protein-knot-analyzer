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

#ifndef PKD_AMALGAMATED_H
#define PKD_AMALGAMATED_H

// c++17
#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <unordered_map>

// c
#include <stdio.h>
#include <string.h>

#ifdef IS_ANALYZER
/* mmdb 1.25.6.1
 * MMDB is a macromolecular coordinate library,
 * written by Eugene Krissinel primarily for use by CCP4 group.
 * The Coordinate Library is designed to assist CCP4 developers in working with coordinate files.
 * The Library features work with the primary file formats of the Protein Data Bank (PDB),
 * the PDB file format and the mmCIF file format
 * License: GNU LGPL v3
 * Documentation: https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_object.html
 */
#include <mmdb_manager.h>
#endif

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
#ifdef IS_ANALYZER
	std::unique_ptr<CMMDBManager> MMDB;
#endif
	std::unordered_map<std::string, std::unique_ptr<DoubleMatrix>> carbonAlphaMatrixMap;
	/* 0: Unknown
	 * 1: MMDB
	 * 2: PKA (This software's minimal management)
	 */
	int managerId;
};

} // namespace PKD

#endif
