#pragma once
#include <d3d11.h>
#include <DirectXMath.h>
#include <Windows.h>
#include <iostream>
#include <fstream>

#include <Render/Render.h>
#include <Render/Vertex.h>
#include <Render/RenderStructures.h>
#include <Math/Bspline.h>
#include <Math/MathStructures.h>

///@todo Make variable length of edges A|B
#define CLOTH_WIDTH 1.
#define CLOTH_HEIGHT 1.
#define GIRD_SIZE 6

constexpr double A = 1. / ((GIRD_SIZE - 1));
constexpr double B = 1./ ((GIRD_SIZE - 1));


///@todo A B now used from global paramters. GIRD SIZE not specified in drapping initial state
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
};

class RenderSys;

///All lines must start at once point
///drappingInit - all initial properties
void drapping_part(RenderSys* _rs, const drappingInit& _is);

//Return angle between lines on surface with nodes (ij, pq) (ij,sh)
double getAngle(const d_vertex* const* gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h);

bool getBSplineDrapPoint(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI);
void getJakobain(double** W, const bSplinePt* ptIJ, const bSplinePt* ptIm1J, const bSplinePt* ptIJm1, const d_vertex* const* IJder);
