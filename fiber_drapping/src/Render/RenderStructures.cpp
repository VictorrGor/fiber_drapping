// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Render/RenderStructures.h"




vertex* convert2DimArrayTo1(vertex** mx, UINT n, UINT m)
{
	vertex* res = new vertex[n * m];
	for (UINT i = 0; i < n; ++i)
		for (UINT j = 0; j < m; ++j)
			res[i * m + j] = mx[i][j];

	return res;
}

