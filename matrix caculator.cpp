#include"matrix caculator.h"


double a[M * N], b[M * N], res[M * N];
struct RANGE range[R];
//240 0 120
//150 60 240
//60 30 150
//120 60 240
//260 30 120?
struct PERSPECTIVE_1 per1 = {30, 100};//u, v;
struct PLOT plot = { 800, 1000, 20, 4, 8, 3, 25, 25, WHITE, BLUE, LIGHTBLUE, BLACK , BLACK, LIGHTGREEN, MAGENTA, false};
struct PLOT3D plot3d = { 800, 1000, 20, 4, 5, 3, 5, 8, 8, 8, 20, 20, 20, perspective_1, {0, 0, 0}, false, WHITE, BLUE, LIGHTBLUE, BLACK, BLUE, BLUE, BLUE, false };
struct OTMP otm = { 4 };
struct PLOT3D_DYNAMIC plot3d_dym = { 5, 5, 5, 0.5, 20 };

bool initmatrix(double *a, int m, int n)
{
	int i, j;

	if (m > N || n > N)
		return false;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i * n + j] = 0;
		}
	}
	return true;
}

bool identity_matrix(double *a, int n)
{
	int i, j;

	if (n > M || n > N)
		return false;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i * n + j] = 0;
		}
	}

	for (i = 0; i < n; i++)
		a[i * n + i] = 1;
	return true;
}


void matrix_copy(double *a, int m, int n, double * res)
{
	int i, j;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			res[i * n + j] = a[i * n + j];
		}
	}
}

void column_copy(double *a, int m, int n, int column, double *res, int res_column, int target_column)
{
	int i;

	res += target_column;
	a += column;
	for (i = 0; i < m; i++)
	{
		*res = *a;
		res += res_column;
		a += n;
	}
}

bool submatrix(double *a, int m, int n, int row_start, int row_length, int column_start, int column_length, double *res)
{

	if (m < row_start + row_length || n < column_start + column_length)
		return false;
	int i, j;
	a += row_start * n + column_start;
	
	for (i = 0; i < row_length; i++)
	{
		for (j = 0; j < column_length; j++)
		{
			res[i * column_length + j] = a[i * n + j];
		}
	}
	return true;
}


void output_matrix(double *a, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("% lf\t", a[i * n + j]);
		}
		printf("\n");
	}
}

void input_matrix(double *a, int m, int n)
{
	int i, j;
	printf("请输入%d个元素, 构成%dx%d矩阵的元素\n", m * n, m, n);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			scanf("%lf", &a[i * n + j]);
		}
	}
	printf("数据输入成功!\n");
}

bool finput_matrix(double *a, int m, int n, char loc[])
{
	FILE * fp;
	if ((fp = fopen(loc, "r")) == NULL)
		return false;

	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fscanf(fp, "%lf", &a[i * n + j]);
		}
	}
	printf("数据输入成功!\n");
	fclose(fp);
	return true;
}

void matrix_sum(double *a, double *b, int m, int n, double *res)
{
	int i, j;
	initmatrix(res, m, n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			res[i * n + j] = a[i * n + j] + b[i * n + j];
		}
	}
}

void matrix_number_multiply(double *a, int m, int n, int times, double *res)
{
	int i, j;
	initmatrix(res, m, n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			res[i * n + j] = times * a[i * n + j];
		}
	}
}

void matrix_multiply(double *a, double *b, int m, int s, int n, double *res)
{
	int i, j, k;
	double *temp1, *temp2;
	temp1 = (double *)malloc(m * s * sizeof(double));
	temp2 = (double *)malloc(s * n * sizeof(double));
	matrix_copy(a, m, s, temp1);
	matrix_copy(b, s, n, temp2);
	initmatrix(res, m, n);

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < s; k++)
			{
				res[i * n + j] += temp1[i * s + k] * temp2[k * n + j];
			}
		}
	}
	free(temp1);
	free(temp2);
}

void matrix_power(double *a, int n, int power, double *res)
{
	int i;
	double *temp = (double *)malloc(n * n * sizeof(double));
	matrix_copy(a, n, n, temp);
	for (i = 1; i < power; i++)
	{
		matrix_multiply(temp, a, n, n, n, temp);
	}
	matrix_copy(temp, n, n, res);
	free(temp);
}

static int factor;
int row_reduce(double *a, int m, int n, double *res, char simp, int model)
{
	matrix_copy(a, m, n, res);

	int i, j, k, s, t, flag, rank = 0, step = 1;
	double temp;
	factor = 1;
	for (i = 0, j = 0; i < m && j < n; j++)
	{
		//找j列非零元
		k = i;
		while (fabs(res[k * n + j]) < 1e-20 && ++k < m);
		if (k == m)
			continue;

		rank++;//存在一个零元, 矩阵的秩加1

		//找到最后一行的非零元后退出消元
		if (i == m - 1 && simp != 's')
			break;

		//将有非零元的一行交换到i行
		if (i != k)
		{
			for (s = j; s < n; s++)
			{
				temp = res[k * n + s];
				res[k * n + s] = res[i * n + s];
				res[i * n + s] = temp;
			}
			if (model == 1)
			{
				printf("step%d(S):\n", step++);
				output_matrix(res, m, n);
			}
			factor *= -1;//每次交换全局变量factor变号;
		}

		//利用主元将j列其他元素消为零
		flag = 0;
		if (simp == 's')
			s = 0;
		else s = i + 1;
		for (; s < m; s++)
		{
			if (s == i || fabs(temp = res[s * n + j] / res[i * n + j]) < 1e-20)
				continue;//保存主元乘的系数系数
			for (t = j; t < n; t++)
			{
				flag = 1;
				res[s * n + t] -= temp * res[i * n + t];
			}
		}

		if (model == 1 && flag == 1)
		{
			printf("step%d(T):\n", step++);
			output_matrix(res, m, n);
		}

		if (simp == 's')
		{
			temp = res[i * n + j];
			for (s = j; s < n; s++)
			res[i * n + s] /= temp;
		}

		i++;
	}
	return rank;
}

int column_reduce(double *a, int m, int n, double *res, char simp, int model)
{
	matrix_copy(a, m, n, res);

	int j, i, k, t, s, flag, rank = 0, step = 1;
	double temp;
	factor = 1;
	for (j = 0, i = 0; j < n && i < m; i++)
	{
		//找j列非零元
		k = j;
		while (fabs(res[k + i * n]) < 1e-20 && ++k < n);
		if (k == n)
			continue;

		rank++;//存在一个零元, 矩阵的秩加1

		//找到最后一行的非零元后退出消元
		if (j == n - 1 && simp != 's')
			break;

		//将有非零元的一行交换到i行
		if (j != k)
		{
			for (t = i; t < m; t++)
			{
				temp = res[k + t * n];
				res[k + t * n] = res[j + t * n];
				res[j + t * n] = temp;
			}
			if (model == 1)
			{
				printf("step%d(S):\n", step++);
				output_matrix(res, m, n);
			}
			factor *= -1;//每次交换全局变量factor变号;
		}

		//利用主元将j列其他元素消为零
		flag = 0;
		if (simp == 's')
			t = 0;
		else t = j + 1;
		for (; t < n; t++)
		{
			if (t == j || fabs(temp = res[t + i * n] / res[j + i * n]) < 1e-20)
				continue;//保存主元乘的系数系数
			for (s = i; s < m; s++)
			{
				flag = 1;
				res[t + s * n] -= temp * res[j + s * n];
			}
		}

		if (model == 1 && flag == 1)
		{
			printf("step%d(T):\n", step++);
			output_matrix(res, m, n);
		}

		if (simp == 's')
		{
			temp = res[j + i * n];
			for (t = i; t < m; t++)
				res[j + t * n] /= temp;
		}

		j++;
	}
	return rank;
}



int rank(double *a, int m, int n)
{
	double *temp;
	int r;
	temp = (double *)malloc(m * n * sizeof(double));
	r = row_reduce(a, m, n, temp, 0, 0);
	free(temp);
	return r;
}

int solve_liner_equation(double *a, int m, int n, double * res)
{
	int i, j, k, dim, rank, pos[N];
	double judge = 0;

	//创建temp矩阵, 将a化为简化阶梯型矩阵后存在temp中
	double *temp;
	temp = (double*)malloc(m * n * sizeof(double));
	if ((rank = row_reduce(a, m, n, temp, 's', 0)) == 0)
	{
		free(temp);
		if (identity_matrix(res, n) == false)
			return -1;
		else
			return n;
	}
	
	//无解
	for (i = 0; i < n - 1; i++)
		judge += fabs(temp[(rank - 1) * n + i]);
	if (judge < 1e-20)
	{
		free(temp);
		return 0;
	}

	//求特解, 解空间向量维数为n - 1
	//res大小为(n - 1) * (n - rank);
	if (initmatrix(res, n - 1, n - rank) == false)
	{
		free(temp);
		return -1;
	}
	//矩阵列满秩, 只有唯一解
	if (rank == n - 1)
	{
		for (i = 0; i < n - 1; i++)
			res[i] = temp[n * i + n - 1];
		free(temp);
		return 1;
	}

	//有多个解:
	//求特解
	dim = n - rank;
	k = 0;
	for (i = 0; i < N; i++)
		pos[i] = 0;
	for (i = 0, j = 0; i < m && j < n - 1; i++, j++)
	{
		while (j < n - 1 && fabs(temp[i * n + j]) < 1e-20)
		{
			res[j * dim] = 0;
			pos[k++] = j++;
		}
		res[j * dim] = temp[i * n + n - 1];
	}
	//将行满秩最后一行主元以后的自由变量记录位置
	while (j < n - 1)
		pos[k++] = j++;

	//求通解, 齐此方程解空间的维数为dim - 1;
	for (k = 1; k < dim; k++)
	{
		for (i = 0, j = 0; i < m && j < n - 1; i++, j++)
		{
			while (fabs(temp[i * n + j]) < 1e-20)
			{
				res[j * dim + k] = 0;
				j++;
			}
			res[j * dim + k] = -temp[i * n + pos[k - 1]];
		}
		res[pos[k - 1] * dim + k] = 1;

	}
	free(temp);
	return dim;
}

int least_square_method(double *a, int m, int n, double *res)
{
	double *at, *temp;
	at = (double*)malloc(m * (n - 1) * sizeof(double));
	temp = (double*)malloc(m * n * sizeof(double));
	int r;
	submatrix(a, m, n, 0, m, 0, n - 1, at);
	transform(at, m, n - 1, at);
	matrix_multiply(at, a, n - 1, m, n, temp);
	r = solve_liner_equation(temp, n - 1, n, res);
	free(at);
	free(temp);
	return r;
}

void output_solution(double *res, int m, int n)
{
	int i, j, t;
	for (j = 0; j < n; j++)
	{
		printf("(");
		for (i = 0; i < m; i++)
		{
			printf("% lf\t", res[i * n + j]);
		}
		printf(")\n");
		if (j < n)
		{
			for (t = 0; t < m / 2; t++)
				printf("  \t\t");
			if (j == 0)
				printf("+\n");
			else if (j == n - 1)
				printf("*T%d\n", j);
			else
				printf("*T%d+\n", j);
		}
	}
}

double laplace_expension(double *a, int n)
{
	int i, j;
	double t = 0;
	double *p = (double *)malloc(sizeof(double)*n*n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			p[i * n + j] = 1;
		}
	}
	for (i = 0; i < n; i++)
		t += ((i % 2 == 0) ? 1 : -1)*a[i] * laplace(a, p, n, 1, i);
	free(p);
	return t;
}

double laplace(double *a, double *p1, int n, int b, int lastreduce)
{

	int i, j = 0;
	double t = 0;
	double *p = (double *)malloc(sizeof(double)*n*n);
	matrix_copy(p1, n, n, p);
	for (i = 0; i < n; i++)
		p[b*n - n + i] = 0;
	for (i = 0; i < n; i++)
		p[lastreduce + i * n] = 0;
	if (b == n - 1)
	{
		for (i = 0; i < n; i++)
		{
			if (p[b*n + i] == 1)
			{
				free(p);
				return a[b*n + i];
			}
		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			if (p[b*n + i] == 1)
			{
				t += ((j % 2 == 0) ? 1 : -1)*a[b*n + i] * laplace(a, p, n, b + 1, i);
				j++;
			}
		}
	}
	free(p);
	return t;
}

bool grammer(double *a, int n, double *res)//Chage
{
	int i, j;
	double *p = (double *)malloc(sizeof(double) * n *n);
	submatrix(a, n, n + 1, 0, n, 0, n, p);
	double b = laplace_expension(p, n);
	if (b < 1e-20)
	{
		free(p);
		return false;
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			p[i + j * n] = a[n + j * (n + 1)];
		res[i] = laplace_expension(p, n) / b;
		for (j = 0; j < n; j++)
			p[i + j * n] = a[i + j * (n + 1)];
	}
	free(p);
	return true;
}

void transform(double *a,int m, int n, double *res)
{

	int i, j;
	double *temp;
	temp = (double *)malloc(m * n * sizeof(double));
	matrix_copy(a, m, n, temp);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			res[j * m + i] = temp[i * n + j];
		}
	}
	free(temp);
}


int offset_standard_by_guass(double *a, int m, int n, double *res)
{
	int i, r;
	double *temp;
	temp = (double *)malloc(m * n * sizeof(double));
	transform(a, m, n, temp);
	row_reduce(temp, n, m, temp, 0, 0);
	transform(temp, n, m, res);
	r = row_reduce(res, m, n, res, 0, 0);

	for (i = 0; i < r; i++)
	{
		res[i * n + i] = 1;
	}

	free(temp);
	return r;
}

int offset_standard(double *a, int m, int n, double *res)
{
	int r, i;
	r = rank(a, m, n);
	initmatrix(res, m, n);

	for (i = 0; i < r; i++)
	{
		res[i * n + i] = 1;
	}
	return r;
}

bool judge_offset(double *a, double *b, int m, int n)
{
	if (rank(a, m, n) == rank(b, m, n))
		return true;
	else
		return false;
}

bool inverse_by_guass(double *a, int n, double *res)
{

	int rank = 0;
	double *augm;
	augm = (double *)malloc(2 * n * n * sizeof(double));
	matrix_copy(a, n, n, augm);
	if (identity_matrix(augm + n * n, n) == false)
	{
		free(augm);
		return false;
	}

	if ((rank = column_reduce(augm, 2 * n, n, augm, 's', 0)) < n)
	{
		free(augm);
		return false;
	}

	matrix_copy(augm + n * n, n, n, res);
	free(augm);
	if (rank == n)
		return true;
}


double determination_by_guass(double *a, int n)
{
	int i;
	double times, *temp;
	temp = (double *)malloc(n * n * sizeof(double));
	factor = 1;
	times = 1;

	if (row_reduce(a, n, n, temp, 0, 0) < n)
	{
		free(temp);
		return 0;
	}

	for (i = 0; i < n; i++)
		times *= temp[i * n + i];
	free(temp);
	return times * factor;
}

int jacobian_method(double *a, int n, double accuracy, double *res)
{
	int i, j, x, y, conter;
	double g0, g, max, s, t, cos_a, sin_a, cos_2a, *temp, *temp1;
	temp = (double*)malloc(n * n * sizeof(double));
	temp1 = (double*)malloc(n * n * sizeof(double));
	matrix_copy(a, n, n, temp);

	s = 0;
	for (i = 0; i < n; i++)
	{
		for (j = i + 1; j < n; j++)
			s += fabs(temp[i * n + j]);
	}
	g = g0 = sqrt(2 * s);
	conter = 0;
	do
	{
		conter++;
		g /= n;
		//找到最大非对角元素
		x = 0, y = 1;
		max = fabs(temp[1]);
		for (i = 0; i < n; i++)
		{
			for (j = i + 1; j < n; j++)
			{
				t = fabs(temp[i * n + j]);
				if (max < t)
				{
					x = i, y = j;
					max = t;
				}
			}
		}
		while (max < g && g >= accuracy * g0 / n)
			g /= n;
		
		//对x行y列元素分别左乘和右乘旋转矩阵
		cos_2a = (temp[y * n + y] - temp[x * n + x]) / sqrt(4 * temp[x * n + y] * temp[x * n + y] + (temp[x * n + x] - temp[y * n + y]) * (temp[x * n + x] - temp[y * n + y]));
		cos_a = sqrt((1 + cos_2a) / 2);
		sin_a = sqrt((1 - cos_2a) / 2);

		matrix_copy(temp, n, n, temp1);
		for (j = 0; j < n; j++)
			temp1[x * n + j] = temp[x * n + j] * cos_a - temp[y * n + j] * sin_a;
		for (j = 0; j < n; j++)
			temp1[y * n + j] = temp[y * n + j] * cos_a + temp[x * n + j] * sin_a;
		matrix_copy(temp1, n, n, temp);
		for (i = 0; i < n; i++)
			temp[i * n + x] = temp1[i * n + x] * cos_a - temp1[i * n + y] * sin_a;
		for (i = 0; i < n; i++)
			temp[i * n + y] = temp1[i * n + y] * cos_a + temp1[i * n + x] * sin_a;

	} while (g >= accuracy * g0 / n );
	matrix_copy(temp, n, n, res);
	free(temp);
	free(temp1);
	return conter;

}

bool qr_decomposition(double *a, int n, double *q, double *r)
{
	if (rank(a, n, n) < n)
		return false;
	int j, t, i;
	double tem, tem1;
	//检查输出地址是否覆盖
	for (j = 0; j < n * n; j++)
	{
		if (q + j == r)
			return false;
	}
	for (j = 0; j < n * n; j++)
	{
		if (r + j == q)
			return false;
	}
	double *temp = (double*)malloc(n * n * sizeof(double));
	double *temp1 = (double*)malloc(n * n * sizeof(double));
	matrix_copy(a, n, n, temp);
	for (j = 0; j < n; j++)
	{
		for (t = 0; t < n; t++)
		{
			*(q + t * n + j) = *(temp + t * n + j);
		}
		for (i = 0; i < j; i++)
		{
			//计算ai * q
			for (tem = 0, t = 0; t < n; t++)
			{
				tem += *(temp + t * n + j) *  *(q + t * n + i);
			}
			//计算q^2;
			for (tem1 = 0, t = 0; t < n; t++)
			{
				tem1 += *(q + t * n + i) *  *(q + t * n + i);
			}
			tem /= tem1;
			//计算qj = |ai * q| / |q^2| * q (第i列)
			for (t = 0; t < n; t++)
			{
				*(q + t * n + j) -= tem * *(q + t * n + i);
			}
		}
		for (tem = 0, t = 0; t < n; t++)
		{
			tem += *(q + t * n + j) * *(q + t * n + j);
		}
		tem = sqrt(tem);
		for (t = 0; t < n; t++)
		{
			*(q + t * n + j) /= tem;
		}
	}
	
	inverse_by_guass(q, n, temp1);
	matrix_multiply(temp1, temp, n, n, n, r);
	free(temp);
	free(temp1);
	return true;
}



/*向量*/

double dot_product(double *vec1, double *vec2, int n)
{
	double ans = 0;
	int i = 0;
	while (i < n)
		ans += vec1[i] * vec2[i++];
	return ans;
}

double cross_product(double *vec1, double *vec2, double *res)
{
	double temp1[N], temp2[N];
	matrix_copy(vec1, 3, 1, temp1);
	matrix_copy(vec2, 3, 1, temp2);
	res[0] = temp1[1] * temp2[2] - temp1[2] * temp2[1];
	res[1] = temp1[2] * temp2[0] - temp1[0] * temp2[2];
	res[2] = temp1[0] * temp2[1] - temp1[1] * temp2[0];
	return vector_module(res, 3);
}

double vector_module(double *vec, int n)
{
	int i = 0;
	double ans = 0;
	while (i < n)
		ans += vec[i] * vec[i++];
	return sqrt(ans);
}

double vector_angle(double *vec1, double *vec2, int n)
{
	const double pi = 2 * acos(0);
	return acos(dot_product(vec1, vec2, n) / vector_module(vec1, n) / vector_module(vec2, n)) / pi * 180;
}

double vector_projection(double *vec1, double *vec2, int n)
{
	return dot_product(vec1, vec2, n) / vector_module(vec2, n);
}




/*应用*/

double fibonacci_sequence(int n)
{
	if (n < 1)
		return 0;
	double b[2] = { 1, 0 };
	double a[4] = { 1, 1, 1, 0 };
	matrix_power(a, 2, n - 1, a);
	matrix_multiply(a, b, 2, 2, 1, b);
	return b[0];
}

double rotate(double *a, int n, double deg, double *res)
{
	const double pi = 2 * acos(0);
	double rot[4] = { cos(deg * 180 / pi), sin(deg * 180 / pi), -sin(deg * 180 / pi), cos(deg * 180 / pi) };
	matrix_multiply(rot, a, 2, 2, n, res);
	return deg;
}

double total_mass(double(*pm)(double *r), struct RANGE range)
{
	register double r[3], m = 0, v = ACCURACY * ACCURACY * ACCURACY;
	register double x1 = range.x1, y1 = range.y1, z1 = range.z1;
	register double x_step = ACCURACY * (range.x1 - range.x0);
	register double y_step = ACCURACY * (range.y1 - range.y0);
	register double z_step = ACCURACY * (range.z1 - range.z0);
	printf("caculating total_mass...\n");
	for (r[0] = range.x0; r[0] <= x1; r[0] += x_step)
	{
		for (r[1] = range.y0; r[1] <= y1; r[1] += y_step)
		{
			for (r[2] = range.z0; r[2] <= z1; r[2] += z_step)
			{
				m += pm(r) * v;
			}
		}
	}
	printf("totla_mass caculation finished!\n");
	return m;
}

void centroid(double m, double(*pm)(double *r), struct RANGE range, double *res)
{
	register double r[3], dm, x = 0, y = 0, z = 0, v = ACCURACY * ACCURACY * ACCURACY;
	register double x1 = range.x1, y1 = range.y1, z1 = range.z1;
	register double x_step = ACCURACY * (range.x1 - range.x0);
	register double y_step = ACCURACY * (range.y1 - range.y0);
	register double z_step = ACCURACY * (range.z1 - range.z0);
	printf("caculating coordination...\n");
	for (r[0] = range.x0; r[0] <= x1; r[0] += x_step)
	{
		for (r[1] = range.y0; r[1] <= y1; r[1] += y_step)
		{
			for (r[2] = range.z0; r[2] <= z1; r[2] += z_step)
			{
				dm = pm(r);
				x += r[0] * dm * v;
				y += r[1] * dm * v;
				z += r[2] * dm * v;
			}
		}
	}
	printf("coordination caculation finished!\n");
	
	res[0] = x / m;
	res[1] = y / m;
	res[2] = z / m;
}

void inertia_tensor(double *cen, double(*pm)(double *r), struct RANGE range, double *res)
{
	register double m, v = ACCURACY * ACCURACY * ACCURACY * (range.x1 - range.x0) * (range.y1 - range.y0) * (range.z1 - range.z0);
	register double temp[9] = { 0 }, r[3];
	register double x1 = range.x1 - cen[0];
	register double y1 = range.y1 - cen[1];
	register double z1 = range.z1 - cen[2];
	register double x_step = ACCURACY * (range.x1 - range.x0);
	register double y_step = ACCURACY * (range.y1 - range.y0);
	register double z_step = ACCURACY * (range.z1 - range.z0);

	printf("caculating inertia_tensor...\n");
	for (r[0] = range.x0 - cen[0]; r[0] <= x1; r[0] += x_step)
	{
		for (r[1] = range.y0 - cen[1]; r[1] <= y1; r[1] += y_step)
		{
			for (r[2] = range.z0 - cen[2]; r[2] <= z1; r[2] += z_step)
			{
				m = pm(r) * v;
				temp[0] += (r[1] * r[1] + r[2] * r[2]) * m;
				temp[1] += r[0] * r[1] * m,					temp[2] += (r[0] * r[0] + r[2] * r[2]) * m;
				temp[3] += r[0] * r[2] * m,					temp[4] += r[1] * r[2] * m,					temp[5] += (r[0] * r[0] + r[1] * r[1]) * m;
			}
		}
	}
	res[0] = temp[0];
	res[1] = res[3] = -temp[1];
	res[2] = res[6] = -temp[3];
	res[4] = temp[2];
	res[5] = res[7] = -temp[4];
	res[8] = temp[5];
	printf("inertia_tensor caculation finished!\n");

}


double x_rotate(double *a, int n, double deg, double *res)
{
	const double pi = 2 * acos(0);

	double rot[9] = { 1,				0,					0,
					   0,		cos(deg * pi / 180), -sin(deg * pi / 180),
					   0,		sin(deg * pi / 180), cos(deg * pi / 180) };
	matrix_multiply(rot, a, 3, 3, n, res);

	return deg;
}

double y_rotate(double *a, int n, double deg, double *res)
{
	const double pi = 2 * acos(0);

	double rot[9] = { cos(deg * pi / 180),	0,	sin(deg * pi / 180),
							 0,				1,			0,
					  -sin(deg * pi / 180),	0,	cos(deg * pi / 180) };
	matrix_multiply(rot, a, 3, 3, n, res);

	return deg;

}

double z_rotate(double *a, int n, double deg, double *res)
{
	const double pi = 2 * acos(0);

	double rot[9] = { cos(deg * pi / 180),	-sin(deg * pi / 180),	0, 
					  sin(deg * pi / 180),	cos(deg * pi / 180),	0,
							  0,					0,				1 };
	matrix_multiply(rot, a, 3, 3, n, res);

	return deg;

}

double rotate3D(double *a, int n, char pos1, double deg1, char pos2, double deg2, char pos3, double deg3, double *res)
{
	int i;
	const double pi = 2 * acos(0);
	char pos[3] = { pos1, pos2, pos3 };
	double deg[3] = { deg1, deg2, deg3 };
	double temp[3 * N], temp1[2 * N];
	matrix_copy(a, 3, n, temp);
	for (i = 0; i < 3; i++)
	{
		switch (pos[i])
		{
		case'x':
			x_rotate(temp, n, deg[i], temp);
			break;
		case'y':
			y_rotate(temp, n, deg[i], temp);
			break;
		case'z':
			z_rotate(temp, n, deg[i], temp);
			break;
		default:
			break;
		}
	}
	column_copy(temp, 3, n, 0, temp1, 1, 0);
	column_copy(a, 3, n, 0, temp1 + 3, 1, 0);
	matrix_copy(temp, 3, n, res);
	return vector_angle(temp1, temp1 + 3, 3);
}

/*可视化*/
//2维
void plot_coordinate(void)
{
	initgraph(plot.L, plot.H);
	setbkcolor(plot.bk_color);
	cleardevice();
	setlinestyle(PS_SOLID, plot.line_thickness);
	setlinecolor(plot.xoy_color);
	//画坐标轴
	settextstyle(plot.text_size << 1, plot.text_size, _T("courier"));
	settextcolor(plot.indicator_color);
	line(plot.padding, plot.H >> 1, plot.L - plot.padding, plot.H >> 1);
	line(plot.L >> 1, plot.padding, plot.L >> 1, plot.H - plot.padding);
	outtextxy(plot.L - plot.text_size / 2 - plot.padding, plot.H / 2 - 2 * plot.text_size - plot.line_thickness, 'x');
	outtextxy(plot.L / 2 + plot.text_size / 2 + plot.line_thickness / 2 + 5, plot.padding, 'y');
	outtextxy(plot.L / 2 + plot.text_size, plot.H / 2 + plot.text_size / 2, '0');
	//标刻度
	settextcolor(plot.scale_color);
	int i;
	TCHAR s[10];
	for (i = 1; i < (plot.L / 2 - plot.padding) / plot.x_scale; i++)
	{
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(plot.L / 2 + i * plot.x_scale - plot.text_size), (int)(plot.H / 2 + plot.text_size / 2), s);
		line((int)(plot.L / 2 + i * plot.x_scale), (int)(plot.H / 2), (int)(plot.L / 2 + i * plot.x_scale), (int)(plot.H / 2 - plot.line_thickness));
	}
	for (i = -1; i > (-plot.L / 2 + plot.padding) / plot.x_scale; i--)
	{
		_stprintf(s, _T("% d"), i);
		outtextxy((int)(plot.L / 2 + i * plot.x_scale - 2 * plot.text_size), (int)(plot.H / 2 + plot.text_size / 2), s);
		line((int)(plot.L / 2 + i * plot.x_scale), (int)(plot.H / 2), (int)(plot.L / 2 + i * plot.x_scale), (int)(plot.H / 2 - plot.line_thickness));
	}
	for (i = 1; i < (plot.H / 2 - plot.padding) / plot.y_scale; i++)
	{
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(plot.L / 2 - (2 + (int)log10(i)) * plot.text_size - plot.line_thickness), (int)(plot.H / 2 - i * plot.y_scale - plot.line_thickness), s);
		line((int)(plot.L / 2), (int)(plot.H / 2 - i * plot.y_scale), (int)(plot.L / 2 + plot.line_thickness), (int)(plot.H / 2 - i * plot.y_scale));
	}
	for (i = -1; i > (-plot.H / 2 + plot.padding) / plot.y_scale; i--)
	{
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(plot.L / 2 - (3 + (int)log10(-i)) * plot.text_size - plot.line_thickness), (int)(plot.H / 2 - i * plot.y_scale - plot.line_thickness), s);
		line((int)(plot.L / 2), (int)(plot.H / 2 - i * plot.y_scale), (int)(plot.L / 2 + plot.line_thickness), (int)(plot.H / 2 - i * plot.y_scale));
	}
	setorigin(plot.L / 2, plot.H / 2);
	setaspectratio(1, -1);
	plot.coordinate = true;
}

void plot_vector(double *a, int n)
{
	if (false == plot.coordinate)
	{
		printf("请先创建坐标轴...\n");
		return;
	}
	setlinecolor(plot.vector_color);
	int i;
	const double pi = 2 * acos(0);
	double x, y, xar, yar, theta = pi / 6, d = 4 * plot.line_thickness;
	for (i = 0; i < n; i++)
	{
		x = *(a + i) * plot.x_scale;
		y = *(a + n + i) * plot.y_scale;
		line(0, 0, (int)x, (int)y);
		xar = x - d * (x * cos(theta) + y * sin(theta)) / sqrt(x * x + y * y);
		yar = y + d * (x * sin(theta) - y * cos(theta)) / sqrt(x * x + y * y);
		line((int)x, (int)y, (int)xar, (int)yar);
		xar = x - d * (x * cos(theta) - y * sin(theta)) / sqrt(x * x + y * y);
		yar = y + d * (-x * sin(theta) - y * cos(theta)) / sqrt(x * x + y * y);
		line((int)x, (int)y, (int)xar, (int)yar);
	}
}

void plot_point(double *a, int n)
{
	if (false == plot.coordinate)
	{
		printf("请先创建坐标轴...\n");
		return;
	}
	setfillcolor(plot.point_color);
	int i;
	double x, y;
	for (i = 0; i < n; i++)
	{
		x = *(a + i) * plot.x_scale;
		y = *(a + n + i) * plot.y_scale;
		solidcircle((int)x, (int)y, plot.point_size);
	}
}

void plot_function(double(*pf)(double x), double x0, double x1)
{
	register double x_step = ACCURACY * (x1 - x0);
	register double x;
	register int t;
	for (x = x0; x < x1 ; x += x_step)
	{
		for (t = -plot.line_thickness / 2; t <= plot.line_thickness / 2; t++)
			putpixel((int)(x * plot.x_scale), (int)(pf(x) * plot.y_scale + t), plot.function_color);
	}
}

void plot_parametric_function(double(*pf)(double t), double(*pg)(double t), double t0, double t1)
{
	register double step = ACCURACY * (t1 - t0);
	register double x;
	register double t;
	for (x = t0; x < t1; x += step)
	{
		for (t = -plot.line_thickness / 2; t <= plot.line_thickness / 2; t++)
			putpixel((int)(pf(x) * plot.x_scale), (int)(pg(x) * plot.y_scale + t), plot.function_color);
	}
}

void plot_implicit_function(double(*pf)(double x, double y), double x0, double x1, double y0, double y1)
{
	srand((unsigned)time(NULL));
	register double x, y;
	register double x_length = x1 - x0;
	register double y_length = y1 - y0;
	register int end = (int)MAX_POINT * (x1 - x0) * (y1 - y0);
	register int t;
	setfillcolor(plot.function_color);
	for (t = 0; t < end; t++)
	{
		x = x_length * rand() / RAND_MAX + x0;
		y = y_length * rand() / RAND_MAX + y0;
		if (fabs(pf(x, y)) < plot3d.line_thickness * 0.01)
		{;
			solidcircle((int)(x * plot.x_scale), (int)(y * plot.y_scale), plot3d.line_thickness / 4 + 1);
		}
	}
}


//3维
void initplot3D(void)
{
	initgraph(plot3d.L, plot3d.H);
	setbkcolor(plot3d.bk_color);
	cleardevice();
	setorigin(plot3d.L >> 1, plot3d.H >> 1);
	plot3d.initplot3d = true;
}

/*
*透视变换
*规定物接近像时透视变换的像的第三个分量为增大
*/
void otmp(double *a, int n, double *res)
{
	double *temp = (double *)malloc(3 * n * sizeof(double));
	rotate3D(a, n, 'z', plot3d.angle3d.phi, 'x', plot3d.angle3d.theta, 'z', plot3d.angle3d.psai, temp);
	double otmp[9] = {  otm.scale / sqrt(2), otm.scale, 0,
						otm.scale / sqrt(2), 0, -otm.scale,
							otm.scale,		 0,  0};
	matrix_multiply(otmp, temp, 3, 3, n, res);
	free(temp);
}

void perspective_1(double *a, int n, double *res)
{
	int i;
	double index;
	double *temp = (double *)malloc(3 * n * sizeof(double));
	rotate3D(a, n, 'z', plot3d.angle3d.phi, 'x', plot3d.angle3d.theta, 'z', plot3d.angle3d.psai, temp);
	for (i = 0; i < n; i++)
	{
		index = (per1.u + per1.v) / (per1.v - *(temp + i));
		*(res + i) = index * *(temp + n + i);
		*(res + n + i) = -index * *(temp + (n << 1)+ i);
		*(res + (n << 1) + i) = *(temp + i);
	}
	free(temp);
}


void plot3D_coordinate(void)
{
	if (false == plot3d.initplot3d)
	{
		printf("请先初始化三维画图...\n");
		return;
	}
	setlinestyle(PS_SOLID, plot3d.line_thickness);
	settextcolor(plot3d.scale_color);
	setlinecolor(plot3d.xoy_color);
	setfillcolor(plot3d.xoy_color);
	int i;
	TCHAR s[10];
	double a[3], r[3] = { 0 };
	//x
	r[1] = 0;
	r[2] = 0;
	int xl = (int)plot3d.x_length;
	int yl = (int)plot3d.y_length;
	int zl = (int)plot3d.z_length;
	if (plot3d.udf_color == false)
	{
		setfillcolor(RED);
		setlinecolor(RED);
	}
	for (i = 1; i < xl; i++)
	{
		r[0] = i * plot3d.x_scale;
		plot3d.per(r, 1, a);
		if (fabs(a[0]) > plot3d.L / 2 - plot3d.padding || fabs(a[1]) > plot3d.H / 2 - plot3d.padding)
			break;
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(a[0] - 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] + 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), s);
		fillcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2);
	}
	r[0] = --i * plot3d.x_scale;
	plot3d.per(r, 1, a);
	line(0, 0, (int)a[0], (int)a[1]);
	settextcolor(plot3d.indicator_color);
	outtextxy((int)(a[0] + 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] - 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), 'x');
	settextcolor(plot3d.scale_color);
	for (i = -1; i > -xl; i--)
	{
		r[0] = i * plot3d.x_scale;
		plot3d.per(r, 1, a);
		if (fabs(a[0]) > plot3d.L / 2 - plot3d.padding || fabs(a[1]) > plot3d.H / 2 - plot3d.padding)
			break;
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(a[0] + 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] - 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), s);
		fillcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2);
	}
	r[0] = ++i * plot3d.x_scale;
	plot3d.per(r, 1, a);
	line(0, 0, (int)a[0], (int)a[1]);
	//y
	r[0] = 0;
	r[2] = 0;
	if (plot3d.udf_color == false)
	{
		setfillcolor(GREEN);
		setlinecolor(GREEN);
	}
	for (i = 1; i < yl; i++)
	{
		r[1] = i * plot3d.y_scale;
		plot3d.per(r, 1, a);
		if (fabs(a[0]) > plot3d.L / 2 - plot3d.padding || fabs(a[1]) > plot3d.H / 2 - plot3d.padding)
			break;
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(a[0] - 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] + 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), s);
		fillcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2);
	}
	r[1] = --i * plot3d.y_scale;
	plot3d.per(r, 1, a);
	line(0, 0, (int)a[0], (int)a[1]);
	settextcolor(plot3d.indicator_color);
	outtextxy((int)(a[0] + 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] - 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), 'y');
	settextcolor(plot3d.scale_color); 
	for (i = -1; i > -yl; i--)
	{
		r[1] = i * plot3d.y_scale;
		plot3d.per(r, 1, a);
		if (fabs(a[0]) > plot3d.L / 2 - plot3d.padding || fabs(a[1]) > plot3d.H / 2 - plot3d.padding)
			break;
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(a[0] + 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] - 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), s);
		fillcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2);
	}
	r[1] = ++i * plot3d.y_scale;
	plot3d.per(r, 1, a);
	line(0, 0, (int)a[0], (int)a[1]);
	//z
	r[0] = 0;
	r[1] = 0;
	if (plot3d.udf_color == false)
	{
		setfillcolor(BLUE);
		setlinecolor(BLUE);
	}
	for (i = 1; i < zl; i++)
	{
		r[2] = i * plot3d.z_scale;
		plot3d.per(r, 1, a);
		if (fabs(a[0]) > plot3d.L / 2 - plot3d.padding || fabs(a[1]) > plot3d.H / 2 - plot3d.padding)
			break;
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(a[0] - 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] + 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), s);
		fillcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2);
	}
	r[2] = --i * plot3d.z_scale;
	plot3d.per(r, 1, a);
	line(0, 0, (int)a[0], (int)a[1]);
	settextcolor(plot3d.indicator_color);
	outtextxy((int)(a[0] + 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] - 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), 'z');
	settextcolor(plot3d.scale_color);
	for (i = -1; i > -zl; i--)
	{
		r[2] = i * plot3d.z_scale;
		plot3d.per(r, 1, a);
		if (fabs(a[0]) > plot3d.L / 2 - plot3d.padding || fabs(a[1]) > plot3d.H / 2 - plot3d.padding)
			break;
		_stprintf(s, _T("%d"), i);
		outtextxy((int)(a[0] + 3 * plot3d.text_size * a[1] / sqrt(a[0] * a[0] + a[1] * a[1])), (int)(a[1] - 3 * plot3d.text_size * a[0] / sqrt(a[0] * a[0] + a[1] * a[1])), s);
		fillcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2);
	}
	r[2] = ++i * plot3d.z_scale;
	plot3d.per(r, 1, a);
	line(0, 0, (int)a[0], (int)a[1]);

}

void plot3D_makerline(double *a, int n)
{
	if (false == plot3d.initplot3d)
	{
		printf("请先初始化三维画图...\n");
		return;
	}
	int i;
	double *temp = (double *)malloc(3 * n * sizeof(double));
	double *am = (double *)malloc(9 * n * sizeof(double));
	double scale[9] = { plot3d.x_scale, 0, 0,
						0, plot3d.y_scale, 0,
						0, 0, plot3d.z_scale };
	matrix_multiply(scale, a, 3, 3, n, temp);
	submatrix(temp, 3, n, 0, 2, 0, n, am);//拷贝(x, y, 0)
	submatrix(temp, 3, n, 0, 1, 0, n, am + 3 * n);//拷贝(x, 0, 0)
	submatrix(temp, 3, n, 1, 1, 0, n, am + 7 * n);//拷贝(0, y, 0)
	for (i = 0; i < n; i++)
	{
		*(am + (n << 1) + i) = *(am + (n << 2) + i) = *(am + 5 * n + i) = *(am + 6 * n + i) = *(am + (n << 3) + i) = 0;
	}
	plot3d.per(temp, n, temp);
	plot3d.per(am, n, am);
	plot3d.per(am + 3 * n, n, am + 3 * n);
	plot3d.per(am + 6 * n, n, am + 6 * n);
	setlinestyle(PS_DASH, plot3d.line_thickness / 2);
	setlinecolor(plot3d.makerline_color);

	for (i = 0; i < n; i++)
	{
		line((int)*(temp + i), (int)*(temp + n + i), (int)*(am + i), (int)*(am + n + i));
		line((int)*(am + i), (int)*(am + n + i), (int)*(am + 6 * n + i), (int)*(am + 7 * n + i));
		line((int)*(am + i), (int)*(am + n + i), (int)*(am + 3 * n + i), (int)*(am + 4 * n + i));
	}
	free(am);
	free(temp);
}

void plot3D_vector(double *a, int n)
{
	if (false == plot3d.initplot3d)
	{
		printf("请先初始化三维画图...\n");
		return;
	}
	double *temp = (double *)malloc(3 * n * sizeof(double));
	double scale[9] = { plot3d.x_scale, 0, 0,
						0, plot3d.y_scale, 0,
						0, 0, plot3d.z_scale };
	matrix_multiply(scale, a, 3, 3, n, temp);
	plot3d.per(temp, n, temp);
	const double pi = 2 * acos(0);
	int i;
	double x, y, xar, yar, theta = pi / 6, d = 4 * plot3d.line_thickness;
	setlinestyle(PS_SOLID, plot3d.line_thickness);
	setlinecolor(plot3d.vector_color);
	for (i = 0; i < n; i++)
	{
		x = *(temp + i);
		y = *(temp + n + i);
		line(0, 0, (int)x, (int)y);
		xar = x - d * (x * cos(theta) + y * sin(theta)) / sqrt(x * x + y * y);
		yar = y + d * (x * sin(theta) - y * cos(theta)) / sqrt(x * x + y * y);
		line((int)x, (int)y, (int)xar, (int)yar);
		xar = x - d * (x * cos(theta) - y * sin(theta)) / sqrt(x * x + y * y);
		yar = y + d * (-x * sin(theta) - y * cos(theta)) / sqrt(x * x + y * y);
		line((int)x, (int)y, (int)xar, (int)yar);
	}
	free(temp);
}

void plot3D_point(double *a, int n)
{
	if (false == plot3d.initplot3d)
	{
		printf("请先初始化三维画图...\n");
		return;
	}
	int i;
	double *temp = (double *)malloc(3 * n * sizeof(double));
	double scale[9] = { plot3d.x_scale, 0, 0,
						0, plot3d.y_scale, 0,
						0, 0, plot3d.z_scale };
	matrix_multiply(scale, a, 3, 3, n, temp);
	plot3d.per(temp, n, temp);
	setfillcolor(plot3d.point_color);
	for (i = 0; i < n; i++)
		solidcircle((int)*(temp + i), (int)*(temp + n + i), plot3d.point_size);
	free(temp);

}

void plot3D_render_binary_function(double(*pf)(double x, double y), double x0, double x1, double y0, double y1)
{
	if (false == plot3d.initplot3d)
	{
		printf("请先初始化三维画图...\n");
		return;
	}
	register double r[3];
	register double am[3];
	register double x_length = (x1 - x0) * plot3d.x_scale;
	register double y_length = (y1 - y0) * plot3d.y_scale;
	x0 *= plot3d.x_scale;
	y0 *= plot3d.y_scale;
	register int t;
	register int end = (int)MAX_POINT * (x1 - x0) * (y1 - y0);
	srand((unsigned)time(NULL));	
	for (t = 0; t < end; t++)
	{
		r[0] = x_length * rand() / RAND_MAX + x0;
		r[1] = x_length * rand() / RAND_MAX + y0;
		r[2] = pf(r[0], r[1]);
		plot3d.per(r, 1, am);
		setfillcolor((int)fabs(plot3d.function_color + am[2]) % 255);
		solidcircle((int)am[0], (int)am[1], plot3d.line_thickness);
	}

}

void plot3D_binary_function(double(*pf)(double x, double y), double x0, double x1, double y0, double y1)
{
	if (false == plot3d.initplot3d)
	{
		printf("请先初始化三维画图...\n");
		return;
	}
	register double r[3], a[3];
	register double x_step = (double)plot3d.grid_spacing / plot3d.x_scale;
	register double y_step = (double)PLOTSTEP / plot3d.y_scale;
	register double xstart, ystart, xend, yend;
	double judge[12] = { x0, x1, 0, 0, 0, 0, y0, y1, 0, 0 , 0, 0};
	plot3d.per(judge, 4, judge);
	if (judge[8] > judge[9])
	{
		xstart = x1;
		xend = x0;
	}
	else
	{
		xstart = x0;
		xend = x1;
	}
	if (judge[10] > judge[11])
	{
		ystart = y1;
		yend = y0;
	}
	else
	{
		ystart = y0;
		yend = y1;
	}
	if (xstart > xend)
		x_step *= -1;
	if (ystart > yend)
		y_step *= -1;

	//register float H, S, L, m;
	//RGBtoHSL(plot3d.function_color, &H, &S, &L);
	r[2] = 0;
	for (r[0] = xstart; xend > xstart && r[0] < xend || xend < xstart && r[0] > xend; r[0] += x_step)
	{
		for (r[1] = ystart; yend > ystart && r[1] < yend || yend < ystart && r[1] > yend; r[1] += y_step)
		{
			a[0] = r[0] * plot3d.x_scale;
			a[1] = r[1] * plot3d.y_scale;
			a[2] = pf(r[0], r[1]) * plot3d.z_scale;
			plot3d.per(a, 1, a);
			setfillcolor((int)fabs(plot3d.function_color + a[2]) % 255);
			solidcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2 + 1);
		}
	}
	for (r[0] = xend, r[1] = ystart; yend > ystart && r[1] < yend || yend < ystart && r[1] > yend; r[1] += y_step)
	{
		a[0] = r[0] * plot3d.x_scale;
		a[1] = r[1] * plot3d.y_scale;
		a[2] = pf(r[0], r[1]) * plot3d.z_scale;
		plot3d.per(a, 1, a);
		setfillcolor((int)fabs(plot3d.function_color + a[2]) % 255);
		solidcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2 + 1);
	}
	x_step = (double)PLOTSTEP / plot3d.x_scale;
	y_step = (double)plot3d.grid_spacing / plot3d.y_scale;
	if (xstart > xend)
		x_step *= -1;
	if (ystart > yend)
		y_step *= -1;
	for (r[1] = ystart; yend > ystart && r[1] < yend || yend < ystart && r[1] > yend; r[1] += y_step)
	{
		for (r[0] = xstart;  xend > xstart && r[0] < xend || xend < xstart && r[0] > xend; r[0] += x_step)
		{
			a[0] = r[0] * plot3d.x_scale;
			a[1] = r[1] * plot3d.y_scale;
			a[2] = pf(r[0], r[1]) * plot3d.z_scale;
			plot3d.per(a, 1, a);
			setfillcolor((int)fabs(plot3d.function_color + a[2]) % 255);
			solidcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2 + 1);
		}
	}
	for (r[0] = xstart, r[1] = yend;  xend > xstart && r[0] < xend || xend < xstart && r[0] > xend; r[0] += x_step)
	{
		a[0] = r[0] * plot3d.x_scale;
		a[1] = r[1] * plot3d.y_scale;
		a[2] = pf(r[0], r[1]) * plot3d.z_scale;
		plot3d.per(a, 1, a);
		setfillcolor((int)fabs(plot3d.function_color + a[2]) % 255);
		solidcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2 + 1);
	}

}

void plot3D_binary_function_onxyplane(double(*pf)(double x, double y), double x0, double x1, double y0, double y1)
{
	register double r[3], a[3], t;
	srand((unsigned)time(NULL));
	register double x_length = x1 - x0;
	register double y_length = y1 - y0;
	register int end = MAX_POINT * (x1 - x0) * (y1 - y0);
	for (t = 0; t < end; t++)
	{
		r[0] = x_length * rand() / RAND_MAX + x0;
		r[1] = y_length * rand() / RAND_MAX + y0;
		if (fabs(pf(r[0], r[1])) < plot3d.line_thickness * 0.01)
		{
			a[0] = r[0] * plot3d.x_scale;
			a[1] = r[1] * plot3d.y_scale;
			a[2] = pf(r[0], r[1]) * plot3d.z_scale;
			plot3d.per(a, 1, a);
			setfillcolor((int)fabs(plot3d.function_color + a[2]) % 255);
			solidcircle((int)a[0], (int)a[1], plot3d.line_thickness / 2 + 1);
		}
	}
}

int plot3d_dynamic(void)
{
	int ch;
	switch (ch = _getch())
	{
	case 72:
		per1.u += plot3d_dym.per1_v;
		otm.scale += plot3d_dym.otmp_scale;
		break;
	case 80:
		per1.u -= plot3d_dym.per1_v;
		otm.scale -= plot3d_dym.otmp_scale;
		break;
	case 75:
		plot3d.angle3d.phi -= plot3d_dym.phi;
		break;
	case 77:
		plot3d.angle3d.phi += plot3d_dym.phi;
		break;
	case 's':
		plot3d.angle3d.theta += plot3d_dym.theta;
		break;
	case'w':
		plot3d.angle3d.theta -= plot3d_dym.theta;
		break;
	case'd':
		plot3d.angle3d.psai += plot3d_dym.psai;
		break;
	case'a':
		plot3d.angle3d.psai -= plot3d_dym.psai;
		break;
	default:
		break;
	}
	return ch;
}