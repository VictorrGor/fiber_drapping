#pragma once

#include "DataStructures.h"
#include "MathLib.h"
#include <DirectXMath.h>
#include <iomanip>

#define LOG_ON

extern size_t splineDegree;

//this func uses first derivation
splineInfo addInterpolationSpline(vertex* vtx, size_t pointCount);
void getInterpolationKnotVector(size_t pointCount, double* forwardU, double* backwardU, size_t knotVectorSize, vertex* vtx);
void BasisFuncs(size_t _i, double _u, size_t _p, double * _knots, double* pRes);
size_t FindSpan(size_t _n, size_t _p, double _u, double * _knot, size_t _knotSize);
double* makeKnotVector(size_t _ptCount, size_t _q);
vertex CurvePoint(splineInfo _spI, size_t _p, double _u);
vertex* makeBSpline(size_t _vxCount, size_t _p, splineInfo _spI);
vertex* getDerivatePoints(size_t _ptCount, vertex* vtx);
vertex* CurveDerivateAlg1(splineInfo spi, size_t p, double u, size_t d);

// (n+1)*(m+1) - data points sizes, Q - points 
//p, q - degree spline surface
//uk, vl - knot vectors
void SurfMeshParams(size_t n, size_t m, vertex** Q, double** uk, double** vl);

//A - q*q square matrix
void LUDecomposition(double** A, size_t q, double*** L, double*** U);

// b - right-handed sizde
double* LUForwardBackward(double** L, double** U, double* b, size_t q);


float* getVxCoordByPosNum(vertex& vx, size_t num);

//p - degree
//r - axises count
vertex* interpolateCurve(vertex* vtx, size_t vtx_ct, size_t p = 3, double* uk = nullptr, double* U = nullptr, size_t r = 3);


//size m = n + p
double* computeKnotVector(double* _u, size_t p, size_t n);

vertex* testSpline();

vertex** transposeMatrix(vertex** mx, size_t n, size_t m);

//mx - transposed m*n matrix, vec - m size vector;  idx -colomn index for saveing row
void saveAsTransponsed(vertex** mx, size_t n, size_t m, size_t idx, vertex* vec);

//Q - n*m vertex array; p, q - spline degree; 
surfInfo GenInterpBSplineSurface(size_t n, size_t m, vertex** Q, size_t p, size_t q);

vertex SurfacePoint(surfInfo* sfI, double u, double v);

double getSplineLen(double _left, double _right, vertex*(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p = 3, size_t n = 100);