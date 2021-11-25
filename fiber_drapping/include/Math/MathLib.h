#pragma once

#include <iostream>

#include <Render/RenderStructures.h>
#include <Math/Bspline.h>




///A - q*q square matrix
void LUDecomposition(double** A, size_t q, double*** L, double*** U);

/// b - right-handed sizde
double* LUForwardBackward(double** L, double** U, double* b, size_t q);