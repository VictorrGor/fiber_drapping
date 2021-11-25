#pragma once

#include <DirectXMath.h>
#include <iomanip>
#include <Render/RenderStructures.h>
#include <Math/MathStructures.h>
#include "MathLib.h"
#include "fstream"

//#define LOG_ON


//this func uses first derivation
splineInfo addInterpolationSpline(d_vertex* vtx, size_t pointCount, size_t splineDegree);
void getInterpolationKnotVector(size_t pointCount, double* forwardU, double* backwardU, size_t knotVectorSize, d_vertex* vtx, size_t splineDegree);
void BasisFuncs(size_t _i, double _u, size_t _p, double * _knots, double* pRes);
//_n - control points count, +p - spline degree, _u - param, _knot - knot, _knotSize - knot size
size_t FindSpan(size_t _n, size_t _p, double _u, double * _knot, size_t _knotSize);
double* makeKnotVector(size_t _ptCount, size_t _q);
d_vertex CurvePoint(splineInfo _spI, size_t _p, double _u);
///Return vertex[_vxCount] array generated as BSpline gird
d_vertex* makeBSpline(size_t _vxCount, size_t _p, splineInfo _spI);
///@todo Deprocated. Set 0-derivates at point 
d_vertex* getDerivatePoints(size_t _ptCount, d_vertex* vtx);
double** DersBasisFuns(size_t i, double u, int p, int n, double* U);
d_vertex* CurveDerivateAlg1(splineInfo spi, size_t p, double u, size_t d);

// (n+1)*(m+1) - data points sizes, Q - points 
//p, q - degree spline surface
//uk, vl - knot vectors
void SurfMeshParams(size_t n, size_t m, d_vertex** Q, double** uk, double** vl);



double* getVxCoordByPosNum(d_vertex& vx, size_t num);

//p - degree
//r - axises count
d_vertex* interpolateCurve(d_vertex* vtx, size_t vtx_ct, size_t p = 3, double* uk = nullptr, double* U = nullptr, size_t r = 3);

//size m = n + p
double* computeKnotVector(double* _u, size_t p, size_t n);

//mx - transposed m*n matrix, vec - m size vector;  idx -colomn index for saveing row
void saveAsTransponsed(d_vertex** mx, size_t n, size_t m, size_t idx, d_vertex* vec);

//Q - n*m vertex array; p, q - spline degree; 
surfInfo GenInterpBSplineSurface(size_t n, size_t m, d_vertex** Q, size_t p, size_t q);

d_vertex SurfacePoint(surfInfo* sfI, double u, double v);

//d - derivaion degree
d_vertex** SurfaceDerivsAlg1(surfInfo* sfI, double u, double v, size_t d);
