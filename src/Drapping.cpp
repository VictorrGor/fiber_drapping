#include "Drapping.h"

vec3** makeGird()
{

	float R = 1;
	vec3 startPoint = vec3(0, 0, R);

	vec3** res = new vec3 * [2];//Моделируются сейчас только оси
	for (int i = 0; i < 2; ++i)
	{
		res[i] = new vec3[GIRD_SIZE];
	}

	float fi = DirectX::XM_PIDIV2;
	float teta = 0;
	float step_teta = DirectX::XM_PI / GIRD_SIZE;
	float step_fi = DirectX::XM_PIDIV2;

	for (size_t i = 0; i < 2; ++i)
	{
		float FI = i * step_fi;
		for (size_t j = 0; j < GIRD_SIZE; ++j)
		{
			float TETA = step_teta * j;
			res[i][j] = vec3(R * sin(FI) * cos(TETA), R * cos(FI) * cos(TETA), R * sin(TETA));
		}
	}
	return res;
}

