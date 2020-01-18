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

} // namespace PKD

#endif
