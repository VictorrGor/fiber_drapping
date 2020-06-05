#pragma once

#include "DataStructures.h"
#include "Bspline.h"

#include <iostream>


//Это нельзя использовать!!! Утечка памяти, ffunc возвращает масив точек, котороый после не алоцируется
//https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D1%82%D1%80%D0%B0%D0%BF%D0%B5%D1%86%D0%B8%D0%B9
vec3 integrate(double _left, double _right, vertex*(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p = 3, size_t n = 100);


double getArea(splineInfo* c1, splineInfo* c2, splineInfo* d1, splineInfo* d2);
//Coons patch area. Расчёт точки на элементарном кусочке площади 
vec3 coonsDerArea(splineInfo* c1, splineInfo* c2, splineInfo* d1, splineInfo* d2, double s, double t);

vertex derivateLcByT(splineInfo* c1, splineInfo* c2, double s);
vertex derivateLcByS(splineInfo* c1, splineInfo* c2, double s, double t);

vertex derivateLdByT(splineInfo* c1, splineInfo* c2, double s, double t);
vertex derivateLdByS(splineInfo* d1, splineInfo* d2, double t);
vertex derivateBbyT(splineInfo* c1, splineInfo* c2, double s, double t); 
vertex derivateBbyS(splineInfo* c1, splineInfo* c2, double s, double t);
