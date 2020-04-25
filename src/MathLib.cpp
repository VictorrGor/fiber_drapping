#include "MathLib.h"

//Это нельзя использовать!!! Утечка памяти, ffunc возвращает масив точек, котороый после не алоцируется
vec3 integrate(double _left, double _right, vertex*(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p, size_t n)
{
	double h = (_right - _left) / (n - 1);
	size_t derDegree = 2;

	vertex res = ffunc(_spi, _p, _left, derDegree)[1].getAbs() / 2. * (_left + h) + ffunc(_spi, _p, _right, derDegree)[1].getAbs() / 2. * (_right - (_left + (n - 1) * h));

	double xi = 0;
	for (int i = 1; i < n - 2; ++i)
	{
		xi = _left + i * h;
		res += ffunc(_spi, _p, xi, derDegree)[1].getAbs() * h;
	}
	return res.pos;
}

vec3 littleFunc(splineInfo * c1, splineInfo * c2, splineInfo * d1, splineInfo * d2, double s)
{
	size_t step_count = 100;
	double step = 1.0 / (step_count - 1);

	vec3 res = vec3(0, 0, 0);
	for (int j = 1; j < step_count - 1; j = j + 2)
	{
		double yj = j * step;

		res = res + coonsDerArea(c1, c2, d1, d2, s, yj - step);
		res = res + coonsDerArea(c1, c2, d1, d2, s, yj + step) * 4;
		res = res + coonsDerArea(c1, c2, d1, d2, s, yj + step);
	}
	res = res * (3. * step);
	return res;
}

double getArea(splineInfo * c1, splineInfo * c2, splineInfo * d1, splineInfo * d2)
{ 
	size_t step_count = 100;

	double step = 1.0 / (step_count - 1);

	vec3 result = vec3(0, 0, 0);
	for (int i = 1; i < step_count - 1; ++i)
	{
		double xi = i * step;
		result = result + littleFunc(c1, c2, d1, d2, xi);
	}
	result = result * (3. * step);
	return getLength(result);
}
//https://en.wikipedia.org/wiki/Coons_patch
//
//Учесть, что используюся разные с1 с2 д1 д2 и переопределить операци
vec3 coonsDerArea(splineInfo * c1, splineInfo * c2, splineInfo * d1, splineInfo * d2,  double s, double t)
{
	//Der by S
	vertex buf = derivateLcByS(c1, c2, s, t) + derivateLdByS(d1, d2, t) - derivateBbyS(c1, c2, s, t);
	vec3 E = buf.pos * buf.pos;

	//std::cout << "E is: " << E.x << " " << E.y << " " << E.z << "\n";
	
	//Der by T
	buf = derivateLcByT(c1, c2, s) + derivateLdByT(d1, d2, s, t) - derivateBbyT(c1, c2, s, t);
	vec3 G = buf.pos * buf.pos;

	//std::cout << "G is: " << G.x << " " << G.y << " " << G.z << "\n";
	
	vec3 F = derivateLcByS(c1, c2, s, t).pos * derivateLcByT(c1, c2, s).pos 
		+ derivateLdByS(d1, d2, t).pos * derivateBbyT(c1, c2, s, t).pos 
		- derivateBbyS(c1, c2, s, t).pos * derivateBbyT(c1, c2, s, t).pos;

	//std::cout << "F is: " << F.x << " " << F.y << " " << F.z << "\n";

	vec3 rr = E*G - F*F;
	//std::cout << "rr is: " << rr.x << " " << rr.y << " " << rr.z << "\n";

	double accuracy = 0.0001;
	if ((rr.x < 0) && (rr.x > -accuracy))
	{
		rr.x = 0;
	}
	if ((rr.y < 0) && (rr.y > -accuracy))
	{
		rr.y = 0;
	}

	if ((rr.z < 0) && (rr.z > -accuracy))
	{
		rr.z = 0;
	}

	if ((rr.x < 0) || (rr.y < 0) || (rr.y < 0))
	{
		std::cout << "rr is: " << rr.x << " " << rr.y << " " << rr.z << "\n";
		std::cout << "Error!";
	}
	return sqrt(rr);
}

vertex derivateLcByT(splineInfo * c1, splineInfo * c2, double s)
{	
	return CurvePoint(*c2, 3, s) - CurvePoint(*c1, 3, s);
}

vertex derivateLcByS(splineInfo* c1, splineInfo* c2, double s, double t)
{
	size_t derDegeree = 2;
	vertex* curveAr1 = CurveDerivateAlg1(*c1, 3, s, derDegeree);
	vertex* curveAr2 = CurveDerivateAlg1(*c2, 3, s, derDegeree);
	
	vertex result = (1 - t) * curveAr1[1] + t * curveAr2[1];
	delete[] curveAr1;
	delete[] curveAr2;
	return result;
}

vertex derivateLdByT(splineInfo* d1, splineInfo* d2, double s, double t)
{
	size_t derDegeree = 2;
	vertex* curveAr1 = CurveDerivateAlg1(*d1, 3, s, derDegeree);
	vertex* curveAr2 = CurveDerivateAlg1(*d2, 3, s, derDegeree);

	vertex result = (1 - s) * curveAr1[1] + s * curveAr2[1];
	delete[] curveAr1;
	delete[] curveAr2;
	return result;
}

vertex derivateLdByS(splineInfo* d1, splineInfo* d2, double t)
{
	return CurvePoint(*d2, 3, t) - CurvePoint(*d1, 3, t);
}

vertex derivateBbyT(splineInfo* c1, splineInfo* c2, double s, double t)
{
	return CurvePoint(*c1, 3, 0) * (s-1) - s * CurvePoint(*c1, 3, 1) + CurvePoint(*c2,3, 0) * (1-s) + CurvePoint(*c2, 3, 1) * s;
}

vertex derivateBbyS(splineInfo * c1, splineInfo * c2, double s, double t)
{
	return CurvePoint(*c1, 3, 0) * (t+1) + (1-t) * CurvePoint(*c1, 3, 1) - CurvePoint(*c2, 3, 0) * t + CurvePoint(*c2, 3, 1) * t;
}


//n -segmentation count
//double getArea(double _leftTeta, double _rightTeta, double _leftFi, double _rightFi, size_t n)
//{
//	double E=0, F=0, G=0;
//
//	double h = (_right - _left) / (n - 1);
//	size_t derDegree = 2;
//
//	vertex res = ffunc(_spi, _p, _left, derDegree)[1].getAbs() / 2. * (_left + h) + ffunc(_spi, _p, _right, derDegree)[1].getAbs() / 2. * (_right - (_left + (n - 1) * h));
//
//	double xi = 0;
//	for (int i = 1; i < n - 2; ++i)
//	{
//		xi = _left + i * h;
//		res += ffunc(_spi, _p, xi, derDegree)[1].getAbs() * h;
//	}
//	return res.pos;
//}
//
//
//double getArea(double _leftTeta, double _rightTeta, double _leftFi, double _rightFi, vertex*(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p, size_t n)
//{
//	vec3 integral = integrate(_leftFi, _rightFi, [](splineInfo c, size_t _p, double t, size_t _dummy) {vertex res = CurvePoint(c, 3, t); return &res;}, _spi, _p, n);
//	return sqrt(integral.x * integral.x + integral.y * integral.y + integral.z * integral.z) * (_rightFi - _leftFi);
//}