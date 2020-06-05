#pragma once

#include "DataStructures.h"
#include "MathLib.h"

extern size_t splineDegree;

splineInfo addInterpolationSpline(vertex* vtx, size_t pointCount);
void getInterpolationKnotVector(size_t pointCount, double* forwardU, double* backwardU, size_t knotVectorSize, vertex* vtx);
void BasisFuncs(size_t _i, double _u, size_t _p, double * _knots, double* pRes);
size_t FindSpan(size_t _n, size_t _p, double _u, double * _knot, size_t _knotSize);
double* makeKnotVector(size_t _ptCount, size_t _q);
vertex CurvePoint(splineInfo _spI, size_t _p, double _u);
vertex * makeBSpline(size_t _vxCount, size_t _p, splineInfo _spI);
vertex* getDerivatePoints(size_t _ptCount, vertex* vtx);
vertex* CurveDerivateAlg1(splineInfo spi, size_t p, double u, size_t d);

double getSplineLen(double _left, double _right, vertex*(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p = 3, size_t n = 100);