#include "Drapping.h"
#include <fstream>

size_t cool_counter = 0;

//pt - Точка в которой считается Якобиан
void getJakobain(float** W, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1)
{
	W[0][0] = 2 * (ptIJ->x - ptIm1J->x);
	W[0][1] = 2 * (ptIJ->y - ptIm1J->y);

	W[1][0] = 2 * (ptIJ->x - ptIJm1->x);
	W[1][1] = 2 * (ptIJ->y - ptIJm1->y);
}

//
void getF(double* f, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1)
{
	f[0] = -(pow((ptIJ->x - ptIm1J->x), 2) + pow((ptIJ->y - ptIm1J->y), 2) + pow((ptIJ->z - ptIm1J->z), 2) - pow(A, 2));
	f[1] = -(pow((ptIJ->x - ptIJm1->x), 2) + pow((ptIJ->y - ptIJm1->y), 2) + pow((ptIJ->z - ptIJm1->z), 2) - pow(B, 2));
}


void calculateObr(float** W, float** invW)
{
	for (size_t i = 0; i < 3; i++)
	{
		memset(invW[i], 0, sizeof(float) * 3);
		invW[i][i] = 1;
	}

	float buf = 0;

	if (cool_counter == 45)
	{
		std::cout << std::endl;
		std::cout << std::endl;
		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				std::cout << W[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}

	//Прямой ход
	for (size_t j = 0; j < 3; ++j)
	{
		if (W[j][j] != 0)
		{
			buf = W[j][j];
			for (size_t i = 0; i < 3; i++)
			{
				invW[j][i] /= buf;
				W[j][i] /= buf;
			}

			for (size_t s = j + 1; s < 3; ++s)
			{
				buf = W[s][j];
				for (size_t q = 0; q < 3; ++q)
				{
					invW[s][q] -= buf * invW[j][q];
					W[s][q] -= buf * W[j][q];
				}
			}
		}
		else
		{
			if ( ( j+1 < 3) && (W[j + 1][j] != 0))
			{
				//Замена ряда на подходящий
				float* swapRow = W[j];
				W[j] = W[j + 1];
				W[j + 1] = swapRow;

				swapRow = invW[j];
				invW[j] = invW[j + 1];
				invW[j + 1] = swapRow;


				//Продолжаем
				buf = W[j][j];
				for (size_t i = 0; i < 3; i++)
				{
					invW[j][i] /= buf;
					W[j][i] /= buf;
				}

				for (size_t s = j + 1; s < 3; ++s)
				{
					buf = W[s][j];
					for (size_t q = 0; q < 3; ++q)
					{
						invW[s][q] -= buf * invW[j][q];
						W[s][q] -= buf * W[j][q];
					}
				}
			}
			else
			{
				if ((j + 2 < 3) && (W[j + 2][j] != 0))
				{
					//Замена ряда на подходящий
					float* swapRow = W[j];
					W[j] = W[j + 2];
					W[j + 2] = swapRow;

					swapRow = invW[j];
					invW[j] = invW[j + 2];
					invW[j + 2] = swapRow;

					//Продолжаем
					buf = W[j][j];
					for (size_t i = 0; i < 3; i++)
					{
						invW[j][i] /= buf;
						W[j][i] /= buf;
					}

					for (size_t s = j + 1; s < 3; ++s)
					{
						buf = W[s][j];
						for (size_t q = 0; q < 3; ++q)
						{
							invW[s][q] -= buf * invW[j][q];
							W[s][q] -= buf * W[j][q];
						}
					}
				}
				else
					throw "Cannot calculate inverse Matrix";
			}
		}
		if (cool_counter == 45)
		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				std::cout << W[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	//Обратный ход
	for (int j = 2; j > 0 ; --j)
	{		
		for (int s = j - 1; s >= 0; --s)
		{
			buf = W[s][j];
			for (size_t q = 0; q < 3; ++q)
			{
				invW[s][q] -= buf * invW[j][q];
				W[s][q] -= buf * W[j][q];
			}
		}
		
	}
}

void printMx(float** W)
{
	std::cout << "\n";
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			std::cout << W[i][j] << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}


void solveSLAU(float** W, float* f, float* dx)
{
	for (size_t i = 0; i < 3; i++)
	{
		dx[i] = f[i];
	}

	float buf = 0;
	float zeroEbs = 0.0001; //Меньше этого числа, число считается нулём

	//Прямой ход
	for (size_t j = 0; j < 3; ++j)
	{
		//printMx(W);
		if ( (W[j][j] != 0) && (fabs(W[j][j]) > zeroEbs))
		{
			buf = W[j][j];

			dx[j] /= buf;
			for (size_t i = j; i < 3; i++)
			{
				W[j][i] /= buf;
			}

			for (size_t s = j + 1; s < 3; ++s)
			{
				buf = W[s][j];

				dx[s] -= buf * dx[j];
				for (size_t q = 0; q < 3; ++q)
				{
					W[s][q] -= buf * W[j][q];
				}
			}
		}
		else
		{
			if ((j + 1 < 3) && (W[j+1][j] != 0) && (fabs(W[j+1][j]) > zeroEbs))
			{
				//Замена ряда на подходящий
				float* swapRow = W[j];
				W[j] = W[j + 1];
				W[j + 1] = swapRow;

				buf = dx[j];
				dx[j] = dx[j + 1];
				dx[j + 1] = buf;


				//Продолжаем
				buf = W[j][j];

				dx[j] /= buf;
				for (size_t i = j; i < 3; i++)
				{
					W[j][i] /= buf;
				}

				for (size_t s = j + 1; s < 3; ++s)
				{
					buf = W[s][j];

					dx[s] -= buf * dx[j];
					for (size_t q = 0; q < 3; ++q)
					{
						W[s][q] -= buf * W[j][q];
					}
				}
			}
			else
			{
				if ((j + 2 < 3) && (W[j+2][j] != 0) && (fabs(W[j+2][j]) > zeroEbs))
				{
					//Замена ряда на подходящий
					float* swapRow = W[j];
					W[j] = W[j + 2];
					W[j + 2] = swapRow;

					buf = dx[j];
					dx[j] = dx[j + 2];
					dx[j + 2] = buf;

					//Продолжаем
					buf = W[j][j];

					dx[j] /= buf;
					for (size_t i = j; i < 3; i++)
					{
						W[j][i] /= buf;
					}

					for (size_t s = j + 1; s < 3; ++s)
					{
						buf = W[s][j];

						dx[s] -= buf * dx[j];
						for (size_t q = 0; q < 3; ++q)
						{
							W[s][q] -= buf * W[j][q];
						}
					}
				}
				else
				{
					std::cout << "Got 0 on main diagonal!\n";
				}
			}
		}
	}
	//Обратный ход
	for (int j = 2; j > 0; --j)
	{
		for (int s = j - 1; s >= 0; --s)
		{
			buf = W[s][j];
			dx[s] -= buf * dx[j];
			for (size_t q = 0; q < 3; ++q)
			{
				W[s][q] -= buf * W[j][q];
			}
		}
	}

	//Проверка на зависимость переменных друг от друга:
	if ((W[0][0] != 1) || (W[1][1] != 1) || (W[2][2] != 1))
	{
		//printMx(W);
		
		if (W[0][0] != 1)
			dx[0] = 0;
		if (W[1][1] != 1)
			dx[1] = 0;
		if (W[2][2] != 1)
			dx[2] = 0;
	}

}

float getDet(float** W)
{
	return W[0][0] * (W[1][1] * W[2][2] - W[1][2] * W[2][1]) - 
		W[0][1] * (W[1][0] * W[2][2] - W[1][2] * W[2][0]) + 
		W[0][2] * (W[1][0] * W[2][1] - W[1][1] * W[2][0]);
}

//Depricated
void minorObr(float** W, float** invW)
{
	float detW = getDet(W);
	std::cout << std::endl;
	std::cout << std::endl;
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			std::cout << W[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	
	if (detW == 0)
	{
		throw "det == 0!\n";
	}
	//Создаём матрицу миноров
	invW[0][0] = W[0][0] * (W[1][1] * W[2][2] - W[1][2] * W[2][1]);
	invW[0][1] = -W[0][1] * (W[1][0] * W[2][2] - W[1][2] * W[2][2]) ;
	invW[0][2] = W[0][2] * (W[1][0] * W[2][1] - W[1][1] * W[2][0]) ;

	invW[1][0] = -W[1][0] * (W[0][1] * W[2][2] - W[0][2] * W[2][1]);
	invW[1][1] = W[1][1] * (W[0][0] * W[2][2] - W[0][2] * W[2][0]) ;
	invW[1][2] = -W[1][2] * (W[0][0] * W[2][1] - W[0][1] * W[2][0]);

	invW[2][0] = W[2][0] * (W[0][1] * W[1][2] - W[0][2] * W[1][1]);
	invW[2][1] = -W[2][1] * (W[0][0] * W[1][2] - W[0][2] * W[1][0]);
	invW[2][2] = W[2][2] * (W[0][0] * W[1][1] - W[0][1] * W[1][0]);
	//Меняем знаки, чтобы сделать матрицу алгебраических дополнений
	invW[0][1] = -invW[0][1];
	invW[1][0] = -invW[1][0];
	invW[1][2] = -invW[1][2];
	invW[2][1] = -invW[2][1];
	//Транспонируем полученную матрицу
	float buf = 0;

	buf = invW[1][0];
	invW[1][0] = invW[0][1];
	invW[0][1] = buf;

	buf = invW[2][0];
	invW[2][0] = invW[0][2];
	invW[0][2] = buf;

	buf = invW[2][1];
	invW[2][1] = invW[1][2];
	invW[1][2] = buf;


	for (size_t i = 0; i < 3; ++i)
		for (size_t j = 0; j < 3; ++j)
			invW[i][j] /= detW;

	std::cout << std::endl;
	std::cout << std::endl;
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			std::cout << invW[i][j] << "\t";
		}
		std::cout << std::endl;
	}

}
//Умножение матрицы на число
void scaleMx(float** mx, float num)
{
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			mx[i][j] *= num;
		}
	}
}

//std::ofstream ppfile;

//Поиск точки на поверхности методом Нюьтона
/*bool getPt(float** W, float** invW, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1)
{
	cool_counter++;
	
	const size_t maxIterCt = 1000;
	bool flag = true;
	size_t i = 0;
	size_t dim = 2;
	float* f = new float[dim];
	float* dx = new float[dim];

	float epsilon = 0.0001;
	float len = 100000;
	float lenMin = len;


	while (flag && (maxIterCt > i))
	{
		getJakobain(W, ptIJ, ptIm1J, ptIJm1);
		///@todo При расчёте обратной матрицы передавать уже готовую, а не каждый раз выделять память
		getF(f, ptIJ, ptIm1J, ptIJm1);
		solveSLAU(W, f, dx);

		
		len = sqrt(pow(dx[0], 2) + pow(dx[1], 2) + pow(dx[2], 2));
		if (len < lenMin) lenMin = len;

		if (len <= epsilon)
		{
			flag = false;
			ptIJ->x += dx[0];
			ptIJ->y += dx[1];
			ptIJ->z += dx[2];
			break;
		}
		ptIJ->x += dx[0];
		ptIJ->y += dx[1];
		ptIJ->z += dx[2];

		++i;
	}

	delete[] f;
	delete[] dx;
	

	return !flag;
}
*/

bool getBSplineDrapPoint(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI)
{
	const size_t maxIterCt = 10000;
	bool flag = true;
	size_t i = 0;

	size_t dim = 2;
	double* f = new double[dim];
	double* dx;

	float epsilon = 0.00001;
	float len = 100000;
	
	bool corrupted = false;

	int counter = 0;

	double** L = new double* [2];
	double** U = new double* [2];
	while (flag && (maxIterCt > i))
	{
		counter++;

		corrupted = false;
		//(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);

		vertex** IJder = SurfaceDerivsAlg1(sfI, ptIJ->u, ptIJ->v, 1);
		getJakobain(W, ptIJ, ptIm1J, ptIJm1, IJder);
		getF(f, &ptIJ->pt->pos, &ptIm1J->pt->pos, &ptIJm1->pt->pos);//@todo Optimize memory using. At each iteration memory allocation happens. It too much heavy
		
		//std::cout << "\n f is: " << f[0] << "\n" << f[1] << "\n";

		LUDecomposition(W, dim, &L, &U);
		dx = LUForwardBackward(L, U, f, 2);
		dx[0] *= 0.001;
		dx[1] *= 0.001;
		//if ((ptIJ->u + dx[0] < 0) || (ptIJ->v + dx[1] < 0) || (ptIJ->u + dx[0] > 1) || (ptIJ->v + dx[1] > 1)) corrupted = true;

		
		if (ptIJ->v + dx[1] < 0)
		{
			while (ptIJ->v + dx[1] < 0)
			{
				//dx[1] = 1 + dx[1];
				dx[1] = 0;
			}
		}
		else
			if (ptIJ->v + dx[1] > 1)
			{
				while (ptIJ->v + dx[1] > 1)
				{
					dx[1] = -1 + dx[1];
					dx[0] += 0.5;
				}
			}
		if (ptIJ->u + dx[0] < 0)
		{
			while (ptIJ->u + dx[0] < 0)
				dx[0] = 1 + dx[0];
		}
		else
			if (ptIJ->u + dx[0] > 1)
			{
				while (ptIJ->u + dx[0] > 1)
					dx[0] = -1. + dx[0];
			}

		vertex bufV = SurfacePoint(sfI, ptIJ->u + dx[0], ptIJ->v + dx[1]);
		len = sqrt(pow((ptIJ->pt->pos.x - bufV.pos.x), 2) + pow((ptIJ->pt->pos.y - bufV.pos.y), 2) + pow((ptIJ->pt->pos.z - bufV.pos.z), 2));
		/*std::cout << "Len is: " << len << "\n";
		std::cout << "A is: " << A << "\nB is: " << B
			<< "\ncaluclate A is: " << pow((ptIJ->pt->pos.x - ptIm1J->pt->pos.x), 2) + pow((ptIJ->pt->pos.y - ptIm1J->pt->pos.y), 2) +
			pow((ptIJ->pt->pos.z - ptIm1J->pt->pos.z), 2)
			<< "\ncaluclate B is: " << pow((ptIJ->pt->pos.x - ptIJm1->pt->pos.x), 2) + pow((ptIJ->pt->pos.y - ptIJm1->pt->pos.y), 2) +
			pow((ptIJ->pt->pos.z - ptIJm1->pt->pos.z), 2) << "\n";*/
			
		ptIJ->u += dx[0];
		ptIJ->v += dx[1];
		(*ptIJ->pt) = bufV;

		for (size_t ct1 = 0; ct1 < dim; ++ct1)
		{
			delete[] L[ct1];
			delete[] U[ct1];
			delete[] IJder[ct1];
		}
		delete[] dx;
		delete[] IJder;

		if ((len <= epsilon) && !corrupted) 
		//if ((fabs(f[0]) < 0.01) && (fabs(f[1]) < 0.01))
		{
			flag = false;
			break;
		}
		++i;

	}
	delete[] L;
	delete[] U;

	delete[] f;

	return !flag;
}

void getJakobain(double** W, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, vertex** IJder)
{	W[0][0] = 2 * (ptIJ->pt->pos.x - ptIm1J->pt->pos.x) * IJder[1][0].pos.x + 2 * (ptIJ->pt->pos.y - ptIm1J->pt->pos.y) * IJder[1][0].pos.y +
		+2 * (ptIJ->pt->pos.z - ptIm1J->pt->pos.z) * IJder[1][0].pos.z;
	W[0][1] = 2 * (ptIJ->pt->pos.x - ptIm1J->pt->pos.x) * IJder[0][1].pos.x + 2 * (ptIJ->pt->pos.y - ptIm1J->pt->pos.y) * IJder[0][1].pos.y +
		+2 * (ptIJ->pt->pos.z - ptIm1J->pt->pos.z) * IJder[0][1].pos.z;

	W[1][0] = 2 * (ptIJ->pt->pos.x - ptIJm1->pt->pos.x) * IJder[1][0].pos.x + 2 * (ptIJ->pt->pos.y - ptIJm1->pt->pos.y) * IJder[1][0].pos.y +
		+2 * (ptIJ->pt->pos.z - ptIJm1->pt->pos.z) * IJder[1][0].pos.z;
	W[1][1] = 2 * (ptIJ->pt->pos.x - ptIJm1->pt->pos.x) * IJder[0][1].pos.x + 2 * (ptIJ->pt->pos.y - ptIJm1->pt->pos.y) * IJder[0][1].pos.y +
		+2 * (ptIJ->pt->pos.z - ptIJm1->pt->pos.z) * IJder[0][1].pos.z;
}

/*
vertex** makeGird()
{

	float R = 1;

	vertex** res = new vertex * [GIRD_SIZE];//Моделируются сейчас только оси
	for (int i = 0; i < GIRD_SIZE; ++i)
	{
		res[i] = new vertex[GIRD_SIZE];
		memset(res[i], 0, sizeof(vertex) * GIRD_SIZE);
	}

	float fi = DirectX::XM_PIDIV2;
	float teta = 0;
	float step_teta = DirectX::XM_PI / (GIRD_SIZE - 1);
	float step_fi = DirectX::XM_PIDIV2;


	//Наложение главных осей на сетку
	float FI = 0;
	size_t index = (GIRD_SIZE - 1) / 2;
	for (size_t j = 0; j < GIRD_SIZE; ++j)
	{
		float TETA = step_teta * j;
		res[index][j] = vertex();
		res[index][j].Color = vec4(1, 0, 0, 1);
		res[index][j].pos = vec3(R * sin(FI) * cos(TETA), R * sin(TETA), R * cos(FI) * cos(TETA));
	}
	FI = DirectX::XM_PIDIV2;

	for (size_t j = 0; j < GIRD_SIZE; ++j)
	{
		float TETA = step_teta * j;
		res[j][index] = vertex();
		res[j][index].Color = vec4(1, 0, 0, 1);
		res[j][index].pos = vec3(R * sin(FI) * cos(TETA), R * sin(TETA), R * cos(FI) * cos(TETA));
	}
	//Попытка вычислить новую точку из трёх других на основе wang 1999.
	float** W = new float*[3]; //Jakobian
	float** invW = new float*[3]; //inverse
	for (size_t i = 0; i < 3; ++i)
	{
		W[i] = new float[3];
		invW[i] = new float[3];
	}

	size_t errCtr = 0;

	//ppfile.open("logg");

	//1
	for (int a_index = (GIRD_SIZE - 1) / 2; a_index >= 1; --a_index)
	{
		for (int b_index = (GIRD_SIZE - 1) / 2; b_index >= 1; --b_index)
		{
			vec3* ptIJ = &res[b_index - 1][a_index - 1].pos;
			vec3* ptIm1J = &res[b_index][a_index - 1].pos;
			vec3* ptIJm1 = &res[b_index - 1][a_index].pos;
			vec3* ptPrevIJ = &res[b_index][a_index].pos;


			ptIJ->x = ptIJm1->x;
			ptIJ->y = -ptIJm1->y ;
			ptIJ->z = ptIm1J->z;

			if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
			{
				++errCtr;
			}
		}
	}
	//2
	for (size_t a_index = (GIRD_SIZE - 1) / 2; a_index < GIRD_SIZE - 1; ++a_index)
	{

		for (int b_index = (GIRD_SIZE - 1) / 2; b_index >= 1; --b_index)
		{
			vec3* ptIJ = &res[b_index - 1][a_index + 1].pos;
			vec3* ptIm1J = &res[b_index][a_index + 1].pos;
			vec3* ptIJm1 = &res[b_index - 1][a_index].pos;
			vec3* ptPrevIJ = &res[b_index][a_index].pos;


			ptIJ->x = ptIJm1->x;
			ptIJ->y = -ptIJm1->y;
			ptIJ->z = ptIm1J->z;

			if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
			{
				++errCtr;
			}
		}
	}
	//3
	for (int a_index = (GIRD_SIZE - 1) / 2; a_index >= 1; --a_index)
	{

		for (size_t b_index = (GIRD_SIZE - 1) / 2; b_index < GIRD_SIZE - 1; ++b_index)
		{
			vec3* ptIJ = &res[b_index + 1][a_index - 1].pos;
			vec3* ptIm1J = &res[b_index][a_index - 1].pos;
			vec3* ptIJm1 = &res[b_index + 1][a_index].pos;
			vec3* ptPrevIJ = &res[b_index][a_index].pos;


			ptIJ->x = ptIJm1->x;
			ptIJ->y = -ptIJm1->y;
			ptIJ->z = ptIm1J->z;

			if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
			{
				++errCtr;
			}
		}
	}
	//4
	for (size_t a_index = (GIRD_SIZE - 1) / 2; a_index < GIRD_SIZE - 1; ++a_index)
	{

		for (size_t b_index = (GIRD_SIZE - 1) / 2; b_index < GIRD_SIZE - 1; ++b_index)
		{
			vec3* ptIJ = &res[b_index + 1][a_index + 1].pos;
			vec3* ptIm1J = &res[b_index][a_index + 1].pos;
			vec3* ptIJm1 = &res[b_index + 1][a_index].pos;
			vec3* ptPrevIJ = &res[b_index][a_index].pos;


			ptIJ->x = ptIJm1->x;
			ptIJ->y = -ptIJm1->y;
			ptIJ->z = ptIm1J->z;

			if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
			{
				++errCtr;
			}
		}
	}
	//ppfile.close();
	std::cout << "Errors count: " << errCtr << "\n";
	system("pause");

	for (size_t i = 0; i < 3; ++i)
	{
		delete[] W[i];
		delete[] invW[i];
	}
	delete[] W;
	delete[] invW;

	return res;
}
*/

//Конвертирует из декартовых в параметрические. Возвращает угол фи. 
float getFiAngle(vertex** gird, size_t i, size_t j)
{
	return atan(gird[i][j].pos.z / gird[i][j].pos.x);
}

//Конвертирует из декартовых в параметрические. Возвращает угол тета.
float getTetaAngle(vertex** gird, size_t i, size_t j)
{

	return atan(sqrt(pow(gird[i][j].pos.x, 2) + pow(gird[i][j].pos.z,2)) / gird[i][j].pos.y);
}

float getAngle(vertex** gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h)
{
	struct { float fi; float teta;} p0, p1, p2;
	XMFLOAT3 v1, v2;//Направление p0p1 и p0p2
	v1.x = gird[p][q].pos.x - gird[i][j].pos.x;
	v1.y = gird[p][q].pos.y - gird[i][j].pos.y;
	v1.z = gird[p][q].pos.z - gird[i][j].pos.z;


	v2.x = gird[s][h].pos.x - gird[i][j].pos.x;
	v2.y = gird[s][h].pos.y - gird[i][j].pos.y;
	v2.z = gird[s][h].pos.z - gird[i][j].pos.z;
	
	XMVECTOR angle = XMVector3AngleBetweenVectors(XMLoadFloat3(&v1), XMLoadFloat3(&v2));
	return angle.m128_f32[0] * 57.2958;
}

