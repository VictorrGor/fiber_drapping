#pragma once
#include <d3d11.h>
#include <DirectXMath.h>
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <thread>
#include <future>

#include <Render/Render.h>
#include <Render/Vertex.h>
#include <Render/RenderStructures.h>
#include <Math/Bspline.h>
#include <Math/MathStructures.h>
#include <Math/MathLib.h>
#include "../TimeBench.h"

#define NO_BENCHMARK_DRAPPING

struct drappingInit {
	surfInfo* sfI; // Spline surface info
	double u1;  //Initial u coordinate of first line
	double v1; // Initial v coordinate of first line
	bool isU1; // Is U coordinate fixed
	double u2; //Initial u coordinate of second line
	double v2; //Initial u coordinate of second line
	bool isU2; //Is U coordiante of second line fixed
	double A; //Edge size by u coordinate
	double B; //Edge size by v coordinate
	size_t gird_size;
};

struct drappingCell {
	bSplinePt* ptIJ; 
	bSplinePt* ptIm1J; 
	bSplinePt* ptIJm1;
	const drappingInit* si;
};

struct drapPointInit
{
	double* f;
	double* L;
	double* U;
	double* dx;
	double* LU_staff_y;
	double* W;
	double* invW;
	DerivationInit* derInit;
};


class RenderSys;

///All lines must start at once point
///drappingInit - all initial properties

//Return angle between lines on surface with nodes (ij, pq) (ij,sh)
double getAngle(const d_vertex* const* gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h);

//bool getBSplineDrapPoint(double** W, double** invW, drappingCell& cell);
void getJakobain(double** W, const bSplinePt* ptIJ, const bSplinePt* ptIm1J, const bSplinePt* ptIJm1, const d_vertex* const* IJder);

void makeDrappedGird(RenderSys* _rs, const drappingInit& _is);
drapPointInit* initDrapPointInitStruct(const surfInfo* sfI, size_t dim, size_t der_degree);
void releaseDrapPointInitStruct(const surfInfo* sfI, drapPointInit* obj, size_t der_degree);
void makeDrappedGird_paralled(RenderSys* _rs, const drappingInit& _is);