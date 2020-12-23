#pragma once
#include <d3d11.h>
#include "DataStructures.h"
#include <iostream>

#define CLOTH_WIDTH 1.
#define CLOTH_HEIGHT 1.
#define GIRD_SIZE 91
//Убрать из макропараметров
constexpr float R = 1;
constexpr float A = DirectX::XM_PI * R / GIRD_SIZE;
constexpr float B = DirectX::XM_PI * R / GIRD_SIZE;


//Retuns woven's gird
vertex** makeGird();

//Return angle between lines on surface with nodes (ij, pq) (ij,sh)
float getAngle(vertex** gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h);