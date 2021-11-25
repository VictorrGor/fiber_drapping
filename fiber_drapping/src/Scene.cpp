// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Scene.h"


void lighting_test(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	//Generate plane
	int plane_size = 100;
	vertex* plane = new vertex[plane_size * plane_size];

	for (int i = 0; i < plane_size; ++i)
	{
		for (int j = 0; j < plane_size; ++j)
		{
			plane[i * plane_size + j].pos = vec3(i * 1. / (plane_size - 1) - 0.5, 0, j * 1. / (plane_size - 1) - 0.5);
			plane[i * plane_size + j].Color = vec4(0, 1, 0.5, 1);
		}
	}

	UINT* indicies = new UINT[(plane_size - 1) * (plane_size - 1) * 6];
	UINT counter = 0;
	vertex* vx1, * vx2, * vx3;
	for (int i = 0; i < plane_size - 1; ++i)
		for (int j = 0; j < plane_size - 1; ++j)
		{
			indicies[counter] = i * plane_size + j;
			indicies[counter + 1] = i * plane_size + j + 1;
			indicies[counter + 2] = (i + 1) * plane_size + j;
			counter += 3;
			vx1 = &plane[i * plane_size + j];
			vx2 = &plane[i * plane_size + j + 1];
			vx3 = &plane[(i + 1) * plane_size + j];
			vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
			vx2->normal = vx1->normal;
			vx3->normal = vx1->normal;

			indicies[counter] = i * plane_size + j + 1;
			indicies[counter + 1] = (i + 1) * plane_size + j + 1;
			indicies[counter + 2] = (i + 1) * plane_size + j;

			vx1 = &plane[i * plane_size + j + 1];
			vx2 = &plane[(i + 1) * plane_size + j + 1];
			vx3 = &plane[(i + 1) * plane_size + j];
			vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
			vx2->normal = vx1->normal;
			vx3->normal = vx1->normal;

			counter += 3;
		}

	Object* plane_obj = new Object(_pDevice, _pVxSh, _pPxSh, plane_size * plane_size, plane, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST, counter, indicies);
	plane_obj->setMaterial(true, vec4(0., 1., 0., 1), vec4(0, 1, 0, 1), vec4(1, 0, 0, 1));
	_rs->pushObject(plane_obj);


	vertex* cube = new vertex[8];
	cube[0].pos = vec3(-0.1, 0, -0.1);
	cube[1].pos = vec3(0.1, 0, -0.1);
	cube[2].pos = vec3(0.1, 0, 0.1);
	cube[3].pos = vec3(-0.1, 0, 0.1);

	cube[4].pos = vec3(-0.1, 0.2, -0.1);
	cube[5].pos = vec3(0.1, 0.2, -0.1);
	cube[6].pos = vec3(0.1, 0.2, 0.1);
	cube[7].pos = vec3(-0.1, 0.2, 0.1);


	UINT* cube_indices = new UINT[36];
	for (int i = 0; i < 3; ++i)
	{
		cube_indices[0 + i * 6] = 0 + i;
		cube_indices[1 + i * 6] = 4 + i;
		cube_indices[2 + i * 6] = 1 + i;
		vx1 = &cube[0 + i];
		vx2 = &cube[4 + i];
		vx3 = &cube[1 + i];
		vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
		vx2->normal = vx1->normal;
		vx3->normal = vx1->normal;

		cube_indices[3 + i * 6] = 1 + i;
		cube_indices[4 + i * 6] = 4 + i;
		cube_indices[5 + i * 6] = 5 + i;
		vx1 = &cube[1 + i];
		vx2 = &cube[4 + i];
		vx3 = &cube[5 + i];
		vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
		vx2->normal = vx1->normal;
		vx3->normal = vx1->normal;
	}
	cube_indices[18] = 3;
	cube_indices[19] = 7;
	cube_indices[20] = 0;
	vx1 = &cube[3];
	vx2 = &cube[7];
	vx3 = &cube[0];
	vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
	vx2->normal = vx1->normal;
	vx3->normal = vx1->normal;

	cube_indices[21] = 7;
	cube_indices[22] = 4;
	cube_indices[23] = 0;
	vx1 = &cube[7];
	vx2 = &cube[4];
	vx3 = &cube[0];
	vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
	vx2->normal = vx1->normal;
	vx3->normal = vx1->normal;

	cube_indices[24] = 5;
	cube_indices[25] = 4;
	cube_indices[26] = 6;
	vx1 = &cube[5];
	vx2 = &cube[4];
	vx3 = &cube[6];
	vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
	vx2->normal = vx1->normal;
	vx3->normal = vx1->normal;
	cube_indices[27] = 4;
	cube_indices[28] = 7;
	cube_indices[29] = 6;
	vx1 = &cube[4];
	vx2 = &cube[7];
	vx3 = &cube[6];
	vx1->normal = calculateTriangleNormal(*vx1, *vx2, *vx3);
	vx2->normal = vx1->normal;
	vx3->normal = vx1->normal;

	Object* cube_obj = new Object(_pDevice, _pVxSh, _pPxSh, 8, cube, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST, 30, cube_indices);
	cube_obj->setMaterial(true, vec4(0., 0.7, 0., 1), vec4(0, 0.7, 0, 1), vec4(1, 0, 0, 1));
	_rs->pushObject(cube_obj);
}

//void test_surface(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
//{
//	UINT n = 10;
//	UINT m = 10;
//	vertex** Q = new vertex * [n];
//	for (UINT i = 0; i < n; ++i) Q[i] = new vertex[m];
//
//	double step_fi = XM_PI * 2 / (n - 1);
//	double step_teta = XM_PIDIV2 / (m - 1);
//	double R = 1;
//
//	for (UINT i = 0; i < n; ++i)
//		for (UINT j = 0; j < m; ++j)
//			Q[i][j].pos = vec3(R * cos(i * step_fi) * cos(step_teta * j), R * sin(step_teta * j), R * cos(j * step_teta) * sin(step_fi * i));
//
//
//	surfInfo sfI = GenInterpBSplineSurface(n, m, Q, 3, 3);
//
//
//	UINT size = 100;
//	vertex* res = new vertex[size * size];
//
//
//	for (UINT i = 0; i < size; ++i)
//		for (UINT j = 0; j < size; ++j)
//		{
//			res[i * size + j] = SurfacePoint(&sfI, 1. / (size - 1) * i, 1. / (size - 1) * j);
//			res[i * size + j].Color = vec4(0, 0, 0, 1);
//		}
//
//	Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, size * size, res, D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
//	//Object* obj = new Object(this, pVxSh, pPxSh, n*m, convert2DimArrayTo1(Q, n, n), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
//	//objects.push_back(obj);
//
//	//Drapping part
//	_rs->drawLineOnBSplineSurface(&sfI, 0, 0, false);
//	_rs->drawLineOnBSplineSurface(&sfI, 0.25, 0, true);
//	//drawLineOnBSplineSurface(&sfI, 0.5, 0, true);
//	//drawLineOnBSplineSurface(&sfI, 0.75, 0, true);
//
//	drapping_part(_rs, &sfI, 0, 0, true, 0.25, 0, true);
//	drapping_part(_rs, &sfI, 0.25, 0, true, 0.5, 0, true);
//	drapping_part(_rs, &sfI, 0.5, 0, true, 0.75, 0, true);
//	drapping_part(_rs, &sfI, 0.75, 0, true, 1, 0, true);
//}
//
//void generateSphere(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
//{
//	int step_cnt_fi = 100;
//	int step_cnt_teta = 100;
//	float R = 1.;
//	vertex** pts = new vertex * [step_cnt_fi];
//	for (int i = 0; i < step_cnt_fi; ++i)
//	{
//		pts[i] = new vertex[step_cnt_teta];
//	}
//
//	vertex bufVx;
//	float step_fi = XM_2PI / (step_cnt_fi - 1);
//	float step_teta = XM_PI / (step_cnt_teta - 1) / 2;
//	for (int i = 0; i < step_cnt_fi; ++i)
//	{
//		for (int j = 0; j < step_cnt_teta; ++j)
//		{
//			bufVx.pos = vec3(R * cos(i * step_fi) * sin(j * step_teta), R * sin(i * step_fi) * sin(j * step_teta), R * cos(j * step_teta));
//			pts[i][j] = bufVx;
//		}
//	}
//
//	Object* pts_raw = new Object(_pDevice, _pVxSh, _pPxSh, step_cnt_fi * step_cnt_teta, convert2DimArrayTo1(pts, step_cnt_fi, step_cnt_teta), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
//	//_rs->pushObject(pts_raw);
//
//	surfInfo sfi = GenInterpBSplineSurface(step_cnt_fi, step_cnt_teta, pts, 3, 3);
//
//	_rs->drawLineOnBSplineSurface(&sfi, 0, 0, true);
//	_rs->drawLineOnBSplineSurface(&sfi, 0.25, 0, true);
//	_rs->drawLineOnBSplineSurface(&sfi, 0.00457597, 0, true);
//
//	drapping_part(_rs, &sfi, 0, 0, true, 0.25, 0, true);
//	
//	for (int i = 0; i < step_cnt_fi; ++i)
//	{
//		delete[] pts[i];
//	}
//	delete[] pts;
//}
//
//void drawSinSurf(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
//{
//	const float R_sphere = 0.0315f; //meters
//	const float alf_cone = 80.5f; //degree
//	const float LA_height = 2.52f; //meters
//	const float jamma = atan(R_sphere / (LA_height - R_sphere));
//
//	float scale_coef = 10;
//
//	//Ellipsoid constants
//	const float a = 1.;
//	const float b = 1.;
//	const float c = 1.;
//
//	//float step_fi = DirectX::XM_PI / ((GIRD_SIZE - 1)); 
//	float fi = 0;
//	float teta = 0;
//	//float step_teta = DirectX::XM_PI / 2 / (GIRD_SIZE / 2 - 1);
//	float x = 0, y = 0, z = 0;
//
//	UINT n = 10;
//	UINT m = 10;
//	vertex** Q = new vertex * [n];
//	for (UINT i = 0; i < n; ++i) Q[i] = new vertex[m];
//
//	float step_fi = XM_PI / (n - 1);
//	float step_teta = XM_PI / (m - 1);
//	float R = 1;
//
//	srand((unsigned)time(NULL));
//	for (UINT i = 0; i < n; ++i)
//		for (UINT j = 0; j < m; ++j)
//		{
//			Q[i][j].pos = vec3(i * step_fi, sin(j) * sin(i) * 0.3, (j * step_teta));
//			//Q[i][j].pos = vec3(i * step_fi, 0, (j * step_teta));
//			Q[i][j].Color = vec4((float)(rand() % 101) / 100, (float)(rand() % 101) / 100, (float)(rand() % 101) / 100, 1);
//		}
//	Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, n * m, convert2DimArrayTo1(Q, n, m), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
//	_rs->pushObject(obj);
//
//
//	surfInfo sfI = GenInterpBSplineSurface(n, m, Q, 3, 3);
//	Q = new vertex * [n];
//	for (UINT i = 0; i < n; ++i) Q[i] = new vertex[m];
//
//	for (UINT i = 0; i < n; ++i)
//		for (UINT j = 0; j < m; ++j)
//		{
//			Q[i][j] = SurfacePoint(&sfI, (float)i / (n - 1), (float)j / (m - 1));
//			Q[i][j].Color = vec4((float)(rand() % 101) / 100, (float)(rand() % 101) / 100, (float)(rand() % 101) / 100, 1);
//		}
//	UINT* indices = new UINT[(n - 1) * (m - 1) * 6];
//	UINT counter = 0;
//	for (UINT i = 0; i < n - 1 - 3; ++i)
//		for (UINT j = 0; j < m - 1 - 3; ++j)
//		{
//			indices[counter] = i * m + j;
//			indices[counter + 1] = i * m + (j + 1);
//			indices[counter + 2] = (i + 1) * m + j;
//
//			indices[counter + 3] = i * m + (j + 1);
//			indices[counter + 4] = (i + 1) * m + (j + 1);
//			indices[counter + 5] = (i + 1) * m + j;
//			counter += 6;
//		}
//
//	obj = new Object(_pDevice, _pVxSh, _pPxSh, n * m, convert2DimArrayTo1(Q, n, m), D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST, counter, indices);
//	//obj->setMaterial(true, vec4(1, 0, 0, 1), vec4(1, 0.5, 0, 1), vec4(1, 0, 0, 1), vec4(1, 0, 0, 1));
//	//objects.push_back(obj);
//
//	//drawLineOnBSplineSurface(&sfI, 0, 0, true);
//	//drawLineOnBSplineSurface(&sfI, 0., 0., false);
//	//drawLineOnBSplineSurface(&sfI, 0., 0., true);
//	//drawLineOnBSplineSurface(&sfI, 0.75, 0., true);
//	drapping_part(_rs, &sfI, 0., 0, true, 0., 0, false);
//	/*drapping_part(&sfI, 0.25, 0, true, 0.5, 0, true);
//	drapping_part(&sfI, 0.5, 0, true, 0.75, 0, true);
//	drapping_part(&sfI, 0.75, 0, true, 1, 0, true);*/
//
//	//drapping_part(&sfI, 0.25, 0, true, 0.5, 0, true);
//	//drawLineOnBSplineSurface(&sfI, 1, 0, true);
//}
//
//
//void generateSurfaceByBSpline(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
//{
//	int size_n = 10;
//	int size_m = 10;
//	vertex** Q = new vertex * [size_n];
//	for (int i = 0; i < size_n; ++i) Q[i] = new vertex[size_m];
//	
//	float step = 1. / size_n;
//	float step_y = XM_PI / (size_n - 1);
//	for (int i = 0; i < size_n; ++i)
//	{
//		for (int j = 0; j < size_m; ++j)
//		{
//			Q[i][j].pos = vec3(step * i,  sin(step_y * i) * 0.1, step * j);
//			Q[i][j].Color = vec4(0, 0, 0, 1);
//		}
//	}
//
//	double* U, *V;
//	vertex** P;
//	int size_u = 100;
//	int size_v = 100;
//	surfInfo surf_inf = BSplineSurface(size_n, size_m, Q, 3, 3, &U, &V, &P, size_u, size_v);
//	//GlobalAproximationBSpline(size_n, size_m, Q, 3, 3, 20, 20, &U, &V, &P);
//	//Surfa
//
//	//Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, size_n * size_m, convert2DimArrayTo1(P, size_u, size_v), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
//	Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, size_u * size_v, convert2DimArrayTo1(P, size_u, size_v), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
//	_rs->pushObject(obj);
//	_rs->drawLineOnBSplineSurface(&surf_inf, 0, 0, false);
//	_rs->drawLineOnBSplineSurface(&surf_inf, 0, 0.5, true);
//	drapping_part(_rs, &surf_inf, 0, 0, false, 0, 0.5, true);
//	/*for (int i = 0; i < size_n; ++i)
//	{
//		delete[] Q[i];
//	}
//	delete[] Q;*/
//}


void testPlane(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	d_vertex** Q = new d_vertex *[2];
	for (int i = 0; i < 2; ++i)
	{
		Q[i] = new d_vertex[2];
	}
	Q[0][0] = { -0.5, 0, -0.5 };
	Q[0][1] = { -0.5, 0, 0.5 };
	Q[1][0] = { 0.5, 0, -0.5 };
	Q[1][1] = { 0.5, 0, 0.5 };
	double* uk = new double[4];
	double* vl = new double[4];
	uk[0] = 0;
	uk[1] = 0;
	uk[2] = 1;
	uk[3] = 1;
	vl[0] = 0;
	vl[1] = 0;
	vl[2] = 1;
	vl[3] = 1;

	surfInfo si;
	si.controlPoints = Q;
	si.m = 2;
	si.n = 2;
	si.p = 1;
	si.q = 1;
	si.Uk = uk;
	si.Vl = vl;

	_rs->drawLineOnBSplineSurface(&si, 0, 0, false);
	_rs->drawLineOnBSplineSurface(&si, 0, 0.5, true);

	size_t gird_size = 6;
	double A = 1. / (gird_size-1);
	double B = 1. / (gird_size-1);
	drappingInit is = { &si, 0, 0.5, true, 0, 0, false, A, B, gird_size};
	drapping_part(_rs, is);

}