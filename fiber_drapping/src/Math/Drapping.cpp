// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Math/Drapping.h"

size_t cool_counter = 0;

//pt - Òî÷êà â êîòîðîé ñ÷èòàåòñÿ ßêîáèàí
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

	//Ïðÿìîé õîä
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
				//Çàìåíà ðÿäà íà ïîäõîäÿùèé
				float* swapRow = W[j];
				W[j] = W[j + 1];
				W[j + 1] = swapRow;

				swapRow = invW[j];
				invW[j] = invW[j + 1];
				invW[j + 1] = swapRow;


				//Ïðîäîëæàåì
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
					//Çàìåíà ðÿäà íà ïîäõîäÿùèé
					float* swapRow = W[j];
					W[j] = W[j + 2];
					W[j + 2] = swapRow;

					swapRow = invW[j];
					invW[j] = invW[j + 2];
					invW[j + 2] = swapRow;

					//Ïðîäîëæàåì
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
	//Îáðàòíûé õîä
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
	float zeroEbs = 0.0001; //Ìåíüøå ýòîãî ÷èñëà, ÷èñëî ñ÷èòàåòñÿ íóë¸ì

	//Ïðÿìîé õîä
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
				//Çàìåíà ðÿäà íà ïîäõîäÿùèé
				float* swapRow = W[j];
				W[j] = W[j + 1];
				W[j + 1] = swapRow;

				buf = dx[j];
				dx[j] = dx[j + 1];
				dx[j + 1] = buf;


				//Ïðîäîëæàåì
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
					//Çàìåíà ðÿäà íà ïîäõîäÿùèé
					float* swapRow = W[j];
					W[j] = W[j + 2];
					W[j + 2] = swapRow;

					buf = dx[j];
					dx[j] = dx[j + 2];
					dx[j + 2] = buf;

					//Ïðîäîëæàåì
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
	//Îáðàòíûé õîä
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

	//Ïðîâåðêà íà çàâèñèìîñòü ïåðåìåííûõ äðóã îò äðóãà:
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
	//Ñîçäà¸ì ìàòðèöó ìèíîðîâ
	invW[0][0] = W[0][0] * (W[1][1] * W[2][2] - W[1][2] * W[2][1]);
	invW[0][1] = -W[0][1] * (W[1][0] * W[2][2] - W[1][2] * W[2][2]) ;
	invW[0][2] = W[0][2] * (W[1][0] * W[2][1] - W[1][1] * W[2][0]) ;

	invW[1][0] = -W[1][0] * (W[0][1] * W[2][2] - W[0][2] * W[2][1]);
	invW[1][1] = W[1][1] * (W[0][0] * W[2][2] - W[0][2] * W[2][0]) ;
	invW[1][2] = -W[1][2] * (W[0][0] * W[2][1] - W[0][1] * W[2][0]);

	invW[2][0] = W[2][0] * (W[0][1] * W[1][2] - W[0][2] * W[1][1]);
	invW[2][1] = -W[2][1] * (W[0][0] * W[1][2] - W[0][2] * W[1][0]);
	invW[2][2] = W[2][2] * (W[0][0] * W[1][1] - W[0][1] * W[1][0]);
	//Ìåíÿåì çíàêè, ÷òîáû ñäåëàòü ìàòðèöó àëãåáðàè÷åñêèõ äîïîëíåíèé
	invW[0][1] = -invW[0][1];
	invW[1][0] = -invW[1][0];
	invW[1][2] = -invW[1][2];
	invW[2][1] = -invW[2][1];
	//Òðàíñïîíèðóåì ïîëó÷åííóþ ìàòðèöó
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
//Óìíîæåíèå ìàòðèöû íà ÷èñëî
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

//Ïîèñê òî÷êè íà ïîâåðõíîñòè ìåòîäîì Íþüòîíà
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
		///@todo Ïðè ðàñ÷¸òå îáðàòíîé ìàòðèöû ïåðåäàâàòü óæå ãîòîâóþ, à íå êàæäûé ðàç âûäåëÿòü ïàìÿòü
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

	float epsilon = 0.001; // Варьировать сравнение ошибки в зависимости от накопленной ошибки. Т.е. если накопилась большая ошибка, то уменьшать макисмально взоможную, или прижимать в другую сторону
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
		vertex oldPoint = *(ptIJ->pt);
		
		if (ptIJ->v + dx[1] < 0)
		{
			dx[1] = 0;
			ptIJ->v = 0;
		}
		else
			if (ptIJ->v + dx[1] > 1)
			{
				dx[1] = 0;
				ptIJ->v = 1;
			}
		if (ptIJ->u + dx[0] < 0)
		{
			dx[0] = 0;
			ptIJ->u = 0;
		}
		else
			if (ptIJ->u + dx[0] > 1)
			{
				dx[0] = 0;
				ptIJ->u = 1;
			}

			
		ptIJ->u += dx[0];
		ptIJ->v += dx[1];
		(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);
		

		for (size_t ct1 = 0; ct1 < dim; ++ct1)
		{
			delete[] L[ct1];
			delete[] U[ct1];
			delete[] IJder[ct1];
		}
		delete[] dx;
		delete[] IJder;

		if ((f[0] <= A * epsilon) && (f[1] <= B * epsilon) && !corrupted)
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

bool getBSplineDrapPoint_with_trace(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI, vertex** traceMx, int& traceCnt)
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

	(*traceMx) = new vertex[maxIterCt];
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
		(*traceMx)[counter] = bufV;

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

	traceCnt = counter;
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


float getFiAngle(vertex** gird, size_t i, size_t j)
{
	return atan(gird[i][j].pos.z / gird[i][j].pos.x);
}


float getTetaAngle(vertex** gird, size_t i, size_t j)
{

	return atan(sqrt(pow(gird[i][j].pos.x, 2) + pow(gird[i][j].pos.z,2)) / gird[i][j].pos.y);
}

float getAngle(vertex** gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h)
{
	struct { float fi; float teta;} p0, p1, p2;
	XMFLOAT3 v1, v2;//Íàïðàâëåíèå p0p1 è p0p2
	v1.x = gird[p][q].pos.x - gird[i][j].pos.x;
	v1.y = gird[p][q].pos.y - gird[i][j].pos.y;
	v1.z = gird[p][q].pos.z - gird[i][j].pos.z;


	v2.x = gird[s][h].pos.x - gird[i][j].pos.x;
	v2.y = gird[s][h].pos.y - gird[i][j].pos.y;
	v2.z = gird[s][h].pos.z - gird[i][j].pos.z;
	
	XMVECTOR angle = XMVector3AngleBetweenVectors(XMLoadFloat3(&v1), XMLoadFloat3(&v2));
	return angle.m128_f32[0] * 57.2958;
}


void drapping_part(RenderSys* _rs, surfInfo* sfI, double u1, double v1, bool isU1, double u2, double v2, bool isU2)
{
	UINT size = GIRD_SIZE;
	bSplinePt** P = new bSplinePt * [size]; //points warper

	std::cout << "A is: " << A << "; B is: " << B << ";\n";
	vertex** Q = new vertex * [size];
	for (UINT i = 0; i < size; ++i)
	{
		P[i] = new bSplinePt[size];

		Q[i] = new vertex[size];
		for (UINT j = 0; j < size; ++j)
		{
			Q[i][j].Color = vec4(0, 1, 0, 1);
			P[i][j].pt = &Q[i][j];
		}
	}

	UINT dim = 2;
	double** W = new double* [dim]; //Jakobian
	double** invW = new double* [dim]; //inverse
	for (UINT i = 0; i < dim; ++i)
	{
		W[i] = new double[dim];
		invW[i] = new double[dim];
	}
	UINT err_ct = 0;

	double cycle_step = 1. / (size - 1);

	//Generating initial lines
	for (UINT i = 0; i < size; ++i)
	{
		if (isU1)
		{
			P[0][i].u = u1;
			P[0][i].v = i * cycle_step;
		}
		else
		{
			P[0][i].u = i * cycle_step;
			P[0][i].v = v1;
		}
		if (isU2)
		{
			P[i][0].u = u2;
			P[i][0].v = i * cycle_step;
		}
		else
		{
			P[i][0].u = i * cycle_step;
			P[i][0].v = v2;
		}
		Q[0][i] = SurfacePoint(sfI, P[0][i].u, P[0][i].v);
		Q[i][0] = SurfacePoint(sfI, P[i][0].u, P[i][0].v);
	}

	double delta_u = 0.01;

	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1/*size - 1*/; ++j)
		{
			if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
				continue;

			bSplinePt* ptIJ = &P[i + 1][j + 1];
			bSplinePt* ptIm1J = &P[i][j + 1];
			bSplinePt* ptIJm1 = &P[i + 1][j];

			ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;//ptIJm1->u;
			ptIJ->v = max(ptIJm1->v, ptIm1J->v);// ptIJm1->v - 2* delta_u;
			if ((ptIJ->u == ptIJm1->u) && (ptIJ->v == ptIJm1->v) || (ptIJ->u == ptIm1J->u) && (ptIJ->v == ptIm1J->v))
			{
				ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;
			}
			if (ptIJ->u > 1) ptIJ->u -= 2 * delta_u;
			if (ptIJ->v > 1) ptIJ->v -= 2 * delta_u;
			if (ptIJ->u < 0) ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;//0;
			if (ptIJ->v < 0) ptIJ->v = min(ptIJm1->v, ptIm1J->v);// 0;// 1 + ptIJ->v;

			(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);

#ifdef _DEBUG
			std::cout << "\ti:" << i << "; j:" << j << "\n";
#endif

			
			if (!getBSplineDrapPoint(W, invW, ptIJ, ptIm1J, ptIJm1, sfI))
			{
#ifdef _DEBUG
				std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
				std::cout << "\tx:" << ptIJ->pt->pos.x << "; y:" << ptIJ->pt->pos.y << "; z:" << ptIJ->pt->pos.z << "\n";
				std::cout << "Distance XJ, X_1J: " << getDistance(ptIJ->pt->pos, ptIm1J->pt->pos) << "\n";
				std::cout << "Distance XJ, XJ_1: " << getDistance(ptIJ->pt->pos, ptIJm1->pt->pos) << "\n";
#endif
				ptIJ->u = -1;
				ptIJ->v = -1;
				++err_ct;
				std::cout << err_ct << "\n";
				std::cout << "\ti:" << i << "; j:" << j << "\n";
			}
			else
			{ 
#ifdef _DEBUG
				std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
				std::cout << "\tx:" << ptIJ->pt->pos.x << "; y:" << ptIJ->pt->pos.y << "; z:" << ptIJ->pt->pos.z << "\n";
				std::cout << "Distance XJ, X_1J: " << getDistance(ptIJ->pt->pos, ptIm1J->pt->pos) << "\n";
				std::cout << "Distance XJ, XJ_1: " << getDistance(ptIJ->pt->pos, ptIJm1->pt->pos) << "\n";
#endif
			}
			ptIJ->pt->Color = vec4(0, 1, 0, 1);
		}
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";


	vertex* triangle = new vertex[3];
	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1; ++j)
		{
			if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
				continue;
			float angle = getAngle(Q, i, j, i + 1, j, i, j + 1);
			float red, green, blue, coeff;
			coeff = 1;
			if (angle < 90)
				red = 1 - angle * coeff / 90;
			else
				red = 0;
			if (angle > 90)
			{
				green = 1 - (angle - 90) * coeff / 90;
				blue = (angle - 90) * coeff / 90;
			}
			else
			{
				green = angle * coeff / 90;
				blue = 0;
			}

			triangle[0] = Q[i][j];
			triangle[1] = Q[i][j + 1];
			triangle[2] = Q[i + 1][j];
			triangle[0].normal = calculateTriangleNormal(Q[i][j], Q[i][j + 1], Q[i + 1][j]);
			triangle[1].normal = triangle[0].normal;
			triangle[2].normal = triangle[0].normal;

			triangle[0].Color = vec4(red, green, blue, 1);
			triangle[1].Color = vec4(red, green, blue, 1);
			triangle[2].Color = vec4(red, green, blue, 1);
			_rs->drawTriangle(triangle);
		}
	}
	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1; ++j)
		{
			float angle = getAngle(Q, i + 1, j + 1, i + 1, j, i, j + 1);
			float red, green, blue, coeff;
			coeff = 1;
			if (angle < 90)
				red = 1 - angle * coeff / 90;
			else
				red = 0;
			if (angle > 90)
			{
				green = 1 - (angle - 90) * coeff / 90;
				blue = (angle - 90) * coeff / 90;
			}
			else
			{
				green = angle * coeff / 90;
				blue = 0;
			}


			vertex* triangle = new vertex[3];
			if ((P[i + 1][j + 1].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
				continue;
			triangle[0] = Q[i + 1][j];
			triangle[1] = Q[i][j + 1];
			triangle[2] = Q[i + 1][j + 1];
			triangle[0].normal = calculateTriangleNormal(Q[i + 1][j], Q[i][j + 1], Q[i + 1][j + 1]);
			triangle[1].normal = triangle[0].normal;
			triangle[2].normal = triangle[0].normal;

			triangle[0].Color = vec4(red, green, blue, 1);
			triangle[1].Color = vec4(red, green, blue, 1);
			triangle[2].Color = vec4(red, green, blue, 1);
			_rs->drawTriangle(triangle);
		}
	}


	for (UINT i = 0; i < size; ++i)
	{
		delete[] Q[i];
		delete[] P[i];
	}
	for (UINT i = 0; i < dim; ++i)
	{
		delete W[i];
		delete invW[i];
	}
	delete[] Q;
	delete[] P;
	delete[] W;
	delete[] invW;
}
