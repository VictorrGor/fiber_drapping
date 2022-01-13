#pragma once

#include <iostream>

#include <Render/RenderStructures.h>
#include <Math/Bspline.h>
#include <Math/MathStructures.h>
#include <math.h>




///A - q*q square matrix
void LUDecomposition(double** A, size_t q, double*** L, double*** U);

void LUDecomposition(double* A, size_t q, double* L, double* U);
/// b - right-handed sizde
double* LUForwardBackward(double** L, double** U, double* b, size_t q);
void LUForwardBackward(double* L, double* U, double* b, size_t q, double* res, double* y);


///A - q*q square matrix
void LUDecomposition(double* A, size_t q, double*** L, double*** U);



//b - right boundary. Left boundary is Pim1j.u (For different sign of funciton of left and right b)
double bisectionU(const surfInfo* sfI, const bSplinePt& Pij, const bSplinePt& Pim1j, double b, double eps, double A);
double bisectionV(const surfInfo* sfI, const bSplinePt& Pij, const bSplinePt& Pijm1, double b, double eps, double B);

//splineSurf, parametriacal initial coordinates, is U fixed
double getBsplineLineLength(surfInfo* sfi, double u, double v, bool isU);