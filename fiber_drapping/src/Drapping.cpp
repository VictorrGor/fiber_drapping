#include "Drapping.h"
#include <fstream>

size_t cool_counter = 0;

//pt - Точка в которой считается Якобиан
void getJakobain(float** W, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1)
{
	W[0][0] = 2 * (ptIJ->x - ptIm1J->x);
	W[0][1] = 2 * (ptIJ->y - ptIm1J->y);
	W[0][2] = 2 * (ptIJ->z - ptIm1J->z);

	W[1][0] = 2 * (ptIJ->x - ptIJm1->x);
	W[1][1] = 2 * (ptIJ->y - ptIJm1->y);
	W[1][2] = 2 * (ptIJ->z - ptIJm1->z);

	W[2][0] = 2 * ptIJ->x;
	W[2][1] = 2 * ptIJ->y;
	W[2][2] = 2 * ptIJ->z;

}

//
void getF(float* f, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1)
{
	f[0] = -(pow((ptIJ->x - ptIm1J->x), 2) + pow((ptIJ->y - ptIm1J->y), 2) + pow((ptIJ->z - ptIm1J->z), 2) - pow(A, 2));
	f[1] = -(pow((ptIJ->x - ptIJm1->x), 2) + pow((ptIJ->y - ptIJm1->y), 2) + pow((ptIJ->z - ptIJm1->z), 2) - pow(B, 2));
	f[2] = -(pow((ptIJ->x), 2) + pow((ptIJ->y), 2) + pow((ptIJ->z), 2) - pow(R, 2));
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
bool getPt(float** W, float** invW, vec3* ptIJ, vec3* ptIm1J, vec3* ptIJm1)
{
	cool_counter++;
	//ptIJ->x = ptIm1J->x;
	//ptIJ->y = ptIJm1->y;
	//ptIJ->z = (ptIm1J->z + ptIJm1->z) / 2;
	
	const size_t maxIterCt = 1000;
	bool flag = true;
	size_t i = 0;
	float* f = new float[3];
	float* dx = new float[3];

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
	
	//ppfile << lenMin << "\n";
	/*if (flag)
	{
		std::cout << "Len is: " << len << "\n";
	}*/

	return !flag;
}

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
	//for (size_t a_index = (GIRD_SIZE - 1) / 2; a_index < GIRD_SIZE - 1; ++a_index)
	//{

	//	for (int b_index = (GIRD_SIZE - 1) / 2; b_index >= 1; --b_index)
	//	{
	//		vec3* ptIJ = &res[b_index - 1][a_index + 1].pos;
	//		vec3* ptIm1J = &res[b_index][a_index + 1].pos;
	//		vec3* ptIJm1 = &res[b_index - 1][a_index].pos;
	//		vec3* ptPrevIJ = &res[b_index][a_index].pos;


	//		ptIJ->x = ptIJm1->x;
	//		ptIJ->y = -ptIJm1->y;
	//		ptIJ->z = ptIm1J->z;

	//		if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
	//		{
	//			++errCtr;
	//		}
	//	}
	//}
	////3
	//for (int a_index = (GIRD_SIZE - 1) / 2; a_index >= 1; --a_index)
	//{

	//	for (size_t b_index = (GIRD_SIZE - 1) / 2; b_index < GIRD_SIZE - 1; ++b_index)
	//	{
	//		vec3* ptIJ = &res[b_index + 1][a_index - 1].pos;
	//		vec3* ptIm1J = &res[b_index][a_index - 1].pos;
	//		vec3* ptIJm1 = &res[b_index + 1][a_index].pos;
	//		vec3* ptPrevIJ = &res[b_index][a_index].pos;


	//		ptIJ->x = ptIJm1->x;
	//		ptIJ->y = -ptIJm1->y;
	//		ptIJ->z = ptIm1J->z;

	//		if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
	//		{
	//			++errCtr;
	//		}
	//	}
	//}
	////4
	//for (size_t a_index = (GIRD_SIZE - 1) / 2; a_index < GIRD_SIZE - 1; ++a_index)
	//{

	//	for (size_t b_index = (GIRD_SIZE - 1) / 2; b_index < GIRD_SIZE - 1; ++b_index)
	//	{
	//		vec3* ptIJ = &res[b_index + 1][a_index + 1].pos;
	//		vec3* ptIm1J = &res[b_index][a_index + 1].pos;
	//		vec3* ptIJm1 = &res[b_index + 1][a_index].pos;
	//		vec3* ptPrevIJ = &res[b_index][a_index].pos;


	//		ptIJ->x = ptIJm1->x;
	//		ptIJ->y = -ptIJm1->y;
	//		ptIJ->z = ptIm1J->z;

	//		if (!getPt(W, invW, ptIJ, ptIm1J, ptIJm1))
	//		{
	//			++errCtr;
	//		}
	//	}
	//}
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

