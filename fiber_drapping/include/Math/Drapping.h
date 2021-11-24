#pragma once
#include <d3d11.h>
#include <DirectXMath.h>
#include <Windows.h>
#include <iostream>
#include <fstream>

#include "../Render/Render.h"
#include "../Render/Vertex.h"
#include "../Render/DataStructures.h"
#include "Bspline.h"

//@todo Make variable length of edges A|B
#define CLOTH_WIDTH 1.
#define CLOTH_HEIGHT 1.
#define GIRD_SIZE 6

constexpr float R_sphere = 1;
//constexpr float A = DirectX::XM_PI/2 * R_sphere / ((GIRD_SIZE - 1));
//constexpr float B = DirectX::XM_PI/2 * R_sphere / ((GIRD_SIZE - 1));
constexpr float A = 1. / ((GIRD_SIZE - 1));
constexpr float B = 1./ ((GIRD_SIZE - 1));


struct drappingInit {
	surfInfo* sfI; // Spline surface info
	double u1;  //Initial u coordinate of first line
	double v1; // Initial v coordinate of first line
	bool isU1; // Is U coordinate fixed
	double u2; //Initial u coordinate of second line
	double v2; //Initial u coordinate of second line
	bool isU2; //Is U coordiante of second line fixed
	float A; //Edge size by u coordinate
	float B; //Edge size by v coordinate
};

class RenderSys;

//All lines must start at once point
void drapping_part(RenderSys* _rs, surfInfo* sfI, double u1, double v1, bool isU1,
	double u2, double v2, bool isU2);


//Retuns woven's gird
//vertex** makeGird();

//Return angle between lines on surface with nodes (ij, pq) (ij,sh)
float getAngle(vertex** gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h);


//bool getPt(float** W, float** invW, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1);

bool getBSplineDrapPoint(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI);
void getJakobain(double** W, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, vertex** IJder);

//Save points on all iterration, to make possible restore how algorithm work
bool getBSplineDrapPoint_with_trace(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI, vertex** traceMx, int& traceCnt);

