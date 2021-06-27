#pragma once
#include <d3d11.h>
#include "DataStructures.h"
#include <iostream>
#include "Bspline.h"

#define CLOTH_WIDTH 1.
#define CLOTH_HEIGHT 1.
#define GIRD_SIZE 100

constexpr float R = 1;
constexpr float A = DirectX::XM_PI/2 * R / ((GIRD_SIZE - 1));
constexpr float B = DirectX::XM_PI/2 * R / ((GIRD_SIZE - 1));


//Retuns woven's gird
//vertex** makeGird();

//Return angle between lines on surface with nodes (ij, pq) (ij,sh)
float getAngle(vertex** gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h);


//bool getPt(float** W, float** invW, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1);

bool getBSplineDrapPoint(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI);
void getJakobain(double** W, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, vertex** IJder);

//Save points on all iterration, to make possible restore how algorithm work
bool getBSplineDrapPoint_with_trace(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI, vertex** traceMx, int& traceCnt);
