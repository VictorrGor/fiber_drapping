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

void generateSphere(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	int step_cnt_fi = 100;
	int step_cnt_teta = 100;
	float R = 1.;
	d_vertex** pts = new d_vertex * [step_cnt_fi];
	for (int i = 0; i < step_cnt_fi; ++i)
	{
		pts[i] = new d_vertex[step_cnt_teta];
	}

	float step_fi = XM_2PI / (step_cnt_fi - 1);
	float step_teta = XM_PI / (step_cnt_teta - 1) / 2;
	for (int i = 0; i < step_cnt_fi; ++i)
	{
		for (int j = 0; j < step_cnt_teta; ++j)
		{
			pts[i][j] = { R * cos(i * step_fi) * sin(j * step_teta), R * sin(i * step_fi) * sin(j * step_teta), R * cos(j * step_teta) };
		}
	}

	//Object* pts_raw = new Object(_pDevice, _pVxSh, _pPxSh, step_cnt_fi * step_cnt_teta, convert2DimArrayTo1(pts, step_cnt_fi, step_cnt_teta), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	//_rs->pushObject(pts_raw);

	surfInfo sfi = GenInterpBSplineSurface(step_cnt_fi, step_cnt_teta, pts, 3, 3, CentripetalMeshParams);

	_rs->drawLineOnBSplineSurface(&sfi, 0, 0, true);
	_rs->drawLineOnBSplineSurface(&sfi, 0.25, 0, true);
	//_rs->drawLineOnBSplineSurface(&sfi, 0.00457597, 0, true);

	size_t gird_size = 100;
	double A = XM_PI * R * XM_PIDIV2 / (XM_PI  * gird_size);
	double B = A;
	drappingInit is = { &sfi, 0, 0, true, 0.25, 0, true, A, B, gird_size };

	std::cout << "Optimized version: \n";
	makeDrappedGird(_rs, is);
	std::cout << "Paralled optimized version: \n";
	makeDrappedGird(_rs, is);
	//makeDrappedGird_paralled(_rs, is);
	for (int i = 0; i < step_cnt_fi; ++i)
	{
		delete[] pts[i];
	}
	delete[] pts;
}


void generateSurfaceByBSpline(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	int size_n = 10;
	int size_m = 10;
	d_vertex** Q = new d_vertex * [size_n];
	
	float step = 1. / (size_n- 1);
	float step_y = XM_PI / (size_n - 1);
	for (int i = 0; i < size_n; ++i)
	{
		Q[i] = new d_vertex[size_m];
		for (int j = 0; j < size_m; ++j)
		{
			Q[i][j] = { step * j,  sin(step_y * j)  * 0.1, step * i };
		}
	}

	double* U, *V;
	vertex** P;
	int size_u = 100;
	int size_v = 100;
	surfInfo surf_inf = GenInterpBSplineSurface(size_n, size_m, Q, 2, 2, CentripetalMeshParams);//BSplineSurface(size_n, size_m, Q, 3, 3, &U, &V, &P, size_u, size_v);
	//GlobalAproximationBSpline(size_n, size_m, Q, 3, 3, 20, 20, &U, &V, &P);
	//Surfa

	//Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, size_n * size_m, convert2DimArrayTo1(P, size_u, size_v), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	//Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, size_u * size_v, convert2DimArrayTo1(P, size_u, size_v), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	//_rs->pushObject(obj);
	_rs->drawLineOnBSplineSurface(&surf_inf, 0, 0, false);
	_rs->drawLineOnBSplineSurface(&surf_inf, 0, 0, true);

	size_t gird_size = 100;
	double A = 1. / gird_size;
	double B = 1.00727 / gird_size;
	drappingInit is = { &surf_inf, 0, 0, false, 0, 0, true, A, B, gird_size };
	makeDrappedGird(_rs, is);
	/*for (int i = 0; i < size_n; ++i)
	{
		delete[] Q[i];
	}
	delete[] Q;*/
}


void testPlane(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	size_t pointCnt = 10;
	d_vertex** Q = new d_vertex *[pointCnt];
	for (int i = 0; i < pointCnt; ++i)
	{
		Q[i] = new d_vertex[pointCnt];
		for (int j = 0; j < pointCnt; ++j)
		{
			Q[i][j] = { -0.5 + i * 1. / (pointCnt - 1) , 0, -0.5 + j * 1. / (pointCnt - 1) };
		}
	}

	size_t p = 3;
	size_t q = 3;
	double* uk = new double[pointCnt + p];
	double* vl = new double[pointCnt + q];


	double* _uk = new double[pointCnt + p];
	double* _vl = new double[pointCnt + q];
	ChordalMeshParams(pointCnt, pointCnt, Q, &_uk, &_vl);
	double* Uk = computeKnotVector(_uk, p, pointCnt);
	double* Vl = computeKnotVector(_vl, q, pointCnt);
	size_t m = pointCnt + p + 1;
	std::cout << "Knot U: \n";
	for (int i = 0; i < m; ++i) std::cout << Uk[i] << " ";
	std::cout << "\n";
	std::cout << "Knot V: \n";
	for (int i = 0; i < m; ++i) std::cout << Vl[i] << " ";
	std::cout << "\n";

	//for (size_t i = 0; i <= p; ++i) uk[i] = 0;
	//for (size_t j = 1; j < pointCnt - p; ++j)
	//{
	//	uk[p + j] += (double)(j) / (pointCnt  - p + 1);
	//}
	//for (size_t i = m - p - 1; i < m; ++i) uk[i] = 1;
	////V knot vector init
	//size_t mv = pointCnt + q + 1;
	//for (size_t i = 0; i <= q; ++i) vl[i] = 0;
	//for (size_t j = 1; j < pointCnt - q; ++j)
	//{
	//	vl[q + j] += (double)(j) / (pointCnt - q + 1);
	//}
	//for (size_t i = m - q - 1; i < m; ++i) vl[i] = 1;


	surfInfo si;
	si.controlPoints = Q;
	si.m = pointCnt;
	si.n = pointCnt;
	si.p = p;
	si.q = q;
	si.Uk = Uk;
	si.Vl = Vl;

	_rs->drawLineOnBSplineSurface(&si, 0, 0, false);
	_rs->drawLineOnBSplineSurface(&si, 0, 0.5, true);

	size_t gird_size = 90;
	double A = 1. / (gird_size-1);
	double B = 1. / (gird_size-1);
	drappingInit is = { &si, 0, 0.5, true, 0, 0, false, A, B, gird_size};
	//makeDrappedGird(_rs, is);

}

void cornerDrapping(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	size_t gird_size = 120;
	size_t half_gird_size = gird_size / 2;
	size_t quater_gs = gird_size / 4;


	vertex** r_gird = new vertex* [gird_size];
	d_vertex** gird = new d_vertex* [gird_size];

	double step_i = 1. / (gird_size - 1);
	double step_j = 1. / (gird_size );

	double edge_size = 0.5;

	int N_TST = 4;
	if (N_TST == 1)
	{
		for (size_t i = 0; i < gird_size; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				r_gird[i][j].pos = vec3(j * edge_size / (gird_size - 1), 0, i * edge_size / (gird_size - 1));
				gird[i][j] = { j * edge_size / (gird_size - 1), 0, i * edge_size / (gird_size - 1) };
			}
		}
	}
	if (N_TST == 2)
	{
		if (gird_size % 2 != 0) throw "Incorrect size of gird_size! Must be a multiple of two \n";
		double x, y, z;
		for (size_t i = 0; i < gird_size / 2; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1) ;
				y = 0;
				z = i * edge_size / (gird_size / 2);
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = {x, y, z };
			}
		}
		for (size_t i = gird_size / 2; i < gird_size; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = (i - gird_size / 2 ) * edge_size / (gird_size / 2 - 1);
				z = edge_size;
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
	}
	if (N_TST == 3)
	{
		if (gird_size % 3 != 0) throw "Incorrect size of gird_size! Must be a multiple of three \n";
		double x, y, z;
		for (size_t i = 0; i < gird_size / 3; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = 0;
				z = i * edge_size / (gird_size / 3);
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
		for (size_t i = gird_size / 3; i < gird_size / 3 * 2; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = (i - gird_size / 3) * edge_size / (gird_size / 3);
				z = edge_size;
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
		for (size_t i = 2 * gird_size / 3; i < gird_size; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = edge_size;
				z = edge_size - (i - 2 * gird_size / 3) * edge_size / (gird_size / 3 - 1);
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
	}

	if (N_TST == 4)
	{
		if (gird_size % 4 != 0) throw "Incorrect size of gird_size! Must be a multiple of four \n";
		double x, y, z;
		for (size_t i = 0; i < quater_gs; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = 0;
				z = i * edge_size / (quater_gs);
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
		for (size_t i = quater_gs; i < 2 * quater_gs; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = (i - quater_gs) * edge_size / (quater_gs);
				z = edge_size;
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
		for (size_t i = 2 * quater_gs; i < 3 * quater_gs; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = edge_size;
				z = edge_size - (i - 2 * quater_gs) * edge_size / (quater_gs);
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
		for (size_t i = 3 * quater_gs; i < gird_size; ++i)
		{
			r_gird[i] = new vertex[gird_size];
			gird[i] = new d_vertex[gird_size];
			for (size_t j = 0; j < gird_size; ++j)
			{
				x = j * edge_size / (gird_size - 1);
				y = edge_size - (i - 3 * quater_gs) * edge_size / (quater_gs - 1);
				z = 0;
				r_gird[i][j].pos = vec3(x, y, z);
				gird[i][j] = { x, y, z };
			}
		}
	}

	



	for (int i = 0; i < gird_size; ++i)
	{
		for (int j = 0; j < gird_size; ++j)
		{
			//std::cout << i << " " << j << ": \n\t" << gird[i][j].x << " " << gird[i][j].y << " " << gird[i][j].z << "\n";
		}
	}

	//for (size_t i = 0; i < gird_size; ++i)
	//{
	//	//r_gird[i] = new vertex[gird_size];
	//	//gird[i] = new d_vertex[gird_size];
	//	
	//	for (int j = 0; j < half_gird_size; ++j)
	//	{
	//		r_gird[i][j].pos = vec3(-0.5 + i * step_i, 0.5 - j * 0.5 / (half_gird_size - 1), 0);
	//		gird[i][j] =           {-0.5 + i * step_i, 0.5 - j * 0.5 / (half_gird_size - 1), 0};
	//		std::cout << i << ", " << j << ": x: " << gird[i][j].x << " y: " << gird[i][j].y << " z: " << gird[i][j].z << ";\n";
	//	}
	//	for (int j = half_gird_size; j < gird_size; ++j)
	//	{
	//		r_gird[i][j].pos = vec3(-0.5 + i * step_i, 0, 0 + (j-half_gird_size+1) * step_j);
	//		gird[i][j]       =     {-0.5 + i * step_i, 0, 0 + (j-half_gird_size+1) * step_j};
	//		std::cout << i << ", " << j << ": x: " << gird[i][j].x << " y: " << gird[i][j].y << " z: " << gird[i][j].z << ";\n";
	//	}
	//}

	Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, gird_size * gird_size, convert2DimArrayTo1(r_gird, gird_size, gird_size), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	_rs->pushObject(obj);
	
	surfInfo si = GenInterpBSplineSurface(gird_size, gird_size, gird, 3, 3, CentripetalMeshParams);
	_rs->drawLineOnBSplineSurface(&si, 0, 0, false);
	double bigLineSize = getBsplineLineLength(&si, 0, 0, false);
	_rs->drawLineOnBSplineSurface(&si, 0, 0, true);
	double smallLineSize = getBsplineLineLength(&si, 0, 0, true);
	std::cout << "Big Line Length: " << bigLineSize << "\n";
	//
	size_t spline_gird_size = 1000;
	double A = bigLineSize / (spline_gird_size - 1);
	double B = smallLineSize / (spline_gird_size - 1);
	
	drappingInit is = { &si, 0, 0, true, 0, 0, false, A, B, spline_gird_size };
	makeDrappedGird(_rs, is);
}

void oldCornerDrapping(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	size_t gird_size = 90;
	size_t half_gird_size = gird_size / 2;
	vertex** r_gird = new vertex * [gird_size];
	d_vertex** gird = new d_vertex * [gird_size];

	double step_i = 1. / (gird_size - 1);
	double step_j = 1. / (gird_size);
	for (size_t i = 0; i < gird_size; ++i)
	{
		r_gird[i] = new vertex[gird_size];
		gird[i] = new d_vertex[gird_size];

		for (int j = 0; j < half_gird_size; ++j)
		{
			r_gird[i][j].pos = vec3(-0.5 + i * step_i, 0.5 - j * 0.5 / (half_gird_size - 1), 0);
			gird[i][j] = { -0.5 + i * step_i, 0.5 - j * 0.5 / (half_gird_size - 1), 0 };
			//std::cout << i << ", " << j << ": x: " << gird[i][j].x << " y: " << gird[i][j].y << " z: " << gird[i][j].z << ";\n";
		}
		for (int j = half_gird_size; j < gird_size; ++j)
		{
			r_gird[i][j].pos = vec3(-0.5 + i * step_i, 0, 0 + (j - half_gird_size + 1) * step_j);
			gird[i][j] = { -0.5 + i * step_i, 0, 0 + (j - half_gird_size + 1) * step_j };
			//std::cout << i << ", " << j << ": x: " << gird[i][j].x << " y: " << gird[i][j].y << " z: " << gird[i][j].z << ";\n";
		}
	}

	Object* obj = new Object(_pDevice, _pVxSh, _pPxSh, gird_size * gird_size, convert2DimArrayTo1(r_gird, gird_size, gird_size), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	_rs->pushObject(obj);

	surfInfo si = GenInterpBSplineSurface(gird_size, gird_size, gird, 3, 3, CentripetalMeshParams);
	_rs->drawLineOnBSplineSurface(&si, 0, 0, false);
	_rs->drawLineOnBSplineSurface(&si, 0, 0, true);
	//
	size_t spline_gird_size = 300;
	double A = 1. / (spline_gird_size - 1);
	double B = 1. / (spline_gird_size - 1);
	drappingInit is = { &si, 0, 0, true, 0, 0, false, A, B, spline_gird_size };
	makeDrappedGird_paralled(_rs, is);
}

void cylinderDrapping(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh)
{
	size_t gird_size = 15;
	vertex** r_gird = new vertex * [gird_size];
	d_vertex** gird = new d_vertex * [gird_size];
	double R = 0.7;
	double H = 2.;
	double shpangout_H = 0.2;
	double shpangout_W = 0.2;

	for (size_t i = 0; i < gird_size; ++i)
	{
		r_gird[i] = new vertex[gird_size];
		gird[i] = new d_vertex[gird_size];
	}

	if (gird_size % 5) throw "cylinderDrapping: Incorrect GirdSize!";
	double x = 0, y = 0, z = 0;
	for (size_t i = 0; i < gird_size; ++i)
	{
		for (size_t j = 0; j < gird_size/5; ++j)
		{
			x = R * cos(XM_PI / (gird_size - 1) * i);
			y = j * (H / 2 - shpangout_W) / (gird_size / 5);
			z = R * sin(XM_PI / (gird_size - 1) * i);

			r_gird[i][j].pos = vec3(x, y, z);
			gird[i][j] = { x, y, z };
		}
		double shpangout_coord = (H / 2 - shpangout_W);
		double actualR;
		for (size_t j = gird_size / 5; j < gird_size * 2 / 5; ++j)
		{
			actualR = R - shpangout_H * (j - gird_size / 5) / (gird_size / 5);
			x = actualR * cos(XM_PI / (gird_size - 1) * i);
			y = shpangout_coord;
			z = actualR * sin(XM_PI / (gird_size - 1) * i);

			r_gird[i][j].pos = vec3(x, y, z);
			gird[i][j] = { x, y, z };
		}
		actualR = R - shpangout_H;
		for (size_t j = 2 * gird_size / 5 ; j < gird_size * 3 / 5; ++j)
		{
			x = actualR * cos(XM_PI / (gird_size - 1) * i);
			y = shpangout_coord + shpangout_W * (j - 2 * gird_size / 5) / (gird_size / 5);
			z = actualR * sin(XM_PI / (gird_size - 1) * i);

			r_gird[i][j].pos = vec3(x, y, z);
			gird[i][j] = { x, y, z };
		}
		shpangout_coord += shpangout_W;
		for (size_t j = 3 * gird_size / 5; j < gird_size * 4 / 5; ++j)
		{
			actualR = R - shpangout_H * (j - 3 * gird_size / 5) / (gird_size / 5);
			x = actualR * cos(XM_PI / (gird_size - 1) * i);
			y = shpangout_coord;
			z = actualR * sin(XM_PI / (gird_size - 1) * i);

			r_gird[i][j].pos = vec3(x, y, z);
			gird[i][j] = { x, y, z };
		}
		for (size_t j = 4 * gird_size / 5; j < gird_size; ++j)
		{
			x = R * cos(XM_PI / (gird_size - 1) * i);
			y = shpangout_coord + (j - 4 * gird_size / 5) * (H / 2 - shpangout_W) / (gird_size / 5 - 1);
			z = R * sin(XM_PI / (gird_size - 1) * i);

			r_gird[i][j].pos = vec3(x, y, z);
			gird[i][j] = { x, y, z };
		}
	}

	Object* object = new Object(_pDevice, _pVxSh, _pPxSh, gird_size * gird_size, convert2DimArrayTo1(r_gird, gird_size, gird_size), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	_rs->pushObject(object);

	for (size_t i = 0; i < gird_size; ++i)
	{
		for (int j = 0; j < gird_size; ++j)
		{
			std::cout << i << ", " << j << ": x: " << gird[i][j].x << " y: " << gird[i][j].y << " z: " << gird[i][j].z << ";\n";
		}
	}


	surfInfo si = GenInterpBSplineSurface(gird_size, gird_size, gird, 3, 3, CentripetalMeshParams);
	double  line1_u = 0, 
			line1_v = 0.5,
			line2_u = 0,
			line2_v = 0;
	bool line1 = false, line2 = true;
	_rs->drawLineOnBSplineSurface(&si, line1_u, line1_v, line1);
	double bigLineSize = getBsplineLineLength(&si, line2_u, line2_v, line2);
	_rs->drawLineOnBSplineSurface(&si, line2_u, line2_v, line2);
	double smallLineSize = getBsplineLineLength(&si, line2_u, line2_v, line2);
	std::cout << "Big Line Length: " << bigLineSize << "\n";

	size_t spline_gird_size = 500;
	double A = bigLineSize / (spline_gird_size - 1);
	double B = smallLineSize / (spline_gird_size - 1);

	drappingInit is = { &si, line2_u, line2_v, line2, line1_u, line1_v, line1, A, B, spline_gird_size };
	//makeDrappedGird_paralled(_rs, is);

}
