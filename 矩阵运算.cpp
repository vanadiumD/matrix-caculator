#include"matrix caculator.h"


extern double a[M * N], b[M * N], res[M * N];
extern struct RANGE range[R];
extern struct PERSPECTIVE_1 per1;
extern struct PLOT plot;
extern struct PLOT3D plot3d;


int main()
{
	char loc1[] = "E:\\manuscript\\test.txt";
	char loc2[] = "E:\\manuscript\\test1.txt";

	int m = 5, n = 5, s;

	finput_matrix(a, m, n, loc1);
	printf("%d\n",finput_matrix(b, m, n, loc2));
#if 0
	output_matrix(a, m, n);
	output_matrix(b, m, n);


	row_reduce(a, m, n, res,0, 1);
	printf("rank = %d\n",rank(a, m, n));
	identity_matrix(res, n);
	output_matrix(res, m, n);

#endif

	//高斯消元
#if 0

	finput_matrix(a, m, n, loc1);
	printf("a = \n");
	output_matrix(a, m, n);
	printf("解空间维度 + 1 = %d\nres = \n",s = solve_liner_equation(a, m, n, res));
	//output_matrix(res, n - 1, s);
	printf("solution = \n");
	output_solution(res, n - 1, s);
	//printf("dim = %d\n", s = least_square_method(a, m, n, res));
	//output_matrix(res, n - 1, s);
	//printf("least_square_method_solution = \n");
	//output_solution(res, n - 1, s);
	//printf("rank = %d, a = \n", row_reduce(a, m, n, a, 's', 0));
	//output_matrix(a, m, n);

#endif
#if 0

	finput_matrix(a, m, n, loc1);
	printf("a = \n");
	output_matrix(a, m, n);
	printf("bool = %d\n", inverse_by_guass(a, m, res));
	printf("a-1=\n");
	output_matrix(res, m, n);
	matrix_multiply(res, a, n, n, n, b);
	printf("a * a- =\n");
	output_matrix(b, n, n);
	matrix_multiply(a, res, n, n, n, res);
	printf("a- * a =\n");
	output_matrix(res, n, n);
	printf("%lf", determination_by_guass(a, n));

#endif

#if 0

	finput_matrix(a, m, n, loc1);
	output_matrix(a, 3, m);
	printf("degree = %lf\n", rotate3D(a, m, 'x', 45, 'y', 60, 'z', 45, a));
	output_matrix(a, 3, m);
	submatrix(a, 3, m, 1, 2, 1, 2, b);
	output_matrix(b, 2, 2);

#endif

#if 0
	double mass;
	range[0].x0 = -5;
	range[0].x1 = 5;
	range[0].y0 = -4;
	range[0].y1 = 4;
	range[0].z0 = -3;
	range[0].z1 = 3;
	mass = total_mass(fm, range[0]);
	printf("%lf\ncentriod=\n",mass);
	centroid(mass, fm, range[0], res);
	output_matrix(res, 3, 1);
	double r[3] = { 0 };
	inertia_tensor(r, fm, range[0], res);
	output_matrix(res, 3, 3);

#endif
#if 0
	input_matrix(a, 4, 4);
	output_matrix(a, 4, 4);
	if (false == inverse_by_guass(a, 3, b))
		printf("矩阵不可逆\n");
	output_matrix(b, 3, 3);
	matrix_multiply(b, a, 3, 3, 3, res);
	output_matrix(res, 3, 3);

#endif

#if 0
	output_matrix(a, m, m);
	printf("number of iterations is %d\n", jacobian_method(a, m, 1e-100, b));
	output_matrix(b, m, m);
#endif

#if 0
	int i;
	for (i = 0; i < 2000; i++)
		printf("%.0lf\t", fibonacci_sequence(i));
#endif
#if 0
	printf("a = \n");
	output_matrix(a, n, n);
	qr_decomposition(a, n, b, res);
	printf("q = \n");
	output_matrix(b, n, n);
	printf("r = \n");
	output_matrix(res, n, n);
	matrix_multiply(b, res, n, n, n, a);
	printf("qr = \n");
	output_matrix(a, n, n);
	transform(b, n, n, a);
	matrix_multiply(a, b, n, n, n, a);
	printf("qtq = \n");
	output_matrix(a, n, n);
#endif
	/*应用测试*/

#if 0

	finput_matrix(a, m, n, loc1);
	printf("a = \n");
	output_matrix(a, m, n);
	rotate3D(a, n, 'x', 30, 'y', 0, 'z', 0, a);
	printf("a'=\n");
	output_matrix(a, m, n);
	rotate3D(a, n, 'x', 90, 'y', 0, 'z', 0, a);
	printf("a'=\n");
	output_matrix(a, m, n);
	rotate3D(a, n, 'x', 0, 'y', 90, 'z', 0, a);
	printf("a'=\n");
	output_matrix(a, m, n);
	rotate3D(a, n, 'x', 0, 'y', 0, 'z', 90, a);
	printf("a'=\n");
	output_matrix(a, m, n);
	system("pause");
#endif
#if 0
	plot_coordinate();
	int i;
	/*for (i = 0; i < 360; i++)
	{
		rotate(a, 1, i, a);
		plot_vector(a, 1);
		Sleep(200);
	}*/
	//plot_function(fx, -20, 20);
	//plot_function(gx, -20, 20);
	//plot_implicit_function(fxy, -8, 8, -8, 8);
	//plot_vector(a, 6);
	//plot_point(a, 6);
	plot_parametric_function(fx, gx, 0, 2 * 3.141592653598979);
	system("pause");
#endif

#if 0
	int i, j, k;
	//TCHAR c[10];
	initplot3D();
	plot3d.per = otmp;
			//	_stprintf(c, _T("%d, %d, %d"), i, j, k);
			//	outtextxy(-4, -4, c);
			//	plot3D_point(a, 3);
	while (1)
	{
		while (1)
		{
			//plot3D_vector(a, 3);
			//plot3D_makerline(a, 3);
			//plot3D_binary_function(fxy, -5, 5, -5, 5);
			plot3D_coordinate();
			if (plot3d_dynamic() == 'p')
				break;;
			clearcliprgn();

			//per1.psai++;
			//per1.theta ++;
			//plot3D_vector(a, 3);
			//plot3D_makerline(a, 3);
			//plot3D_binary_function(fxy, -5, 5, -5, 5);
			plot3D_coordinate();
			Sleep(100);
		}

		//plot3D_binary_function(gxy, -3, 3, -5, 5);
		//plot3D_binary_function_onxyplane(gxy, -3, 3, -5, 5);
		plot3D_binary_function(fxy, -5, 5, -5, 5);
		plot3D_coordinate();
		system("pause");
	}



#endif

# if 1
	initplot3D();
	plot3d.per = otmp;
	while (1)
	{
		plot3D_coordinate();
		cuboid(-30, -40, -50, 30, 40, 50);
		plot3d_dynamic();
		clearcliprgn();
		plot3D_coordinate();
		cuboid(-30, -40, -50, 30, 40, 50);
		//per1.psai++;
		//per1.theta ++;
		Sleep(100);
	}

# endif

#if 0
	finput_matrix(a, m, n, loc1);
	printf("a = \n");
	output_matrix(a, m, n);
	perspective_1(a, n, res);
	printf("res = \n");
	output_matrix(res, 2, n);
	system("pause");
#endif
	//external test
#if 0

	double a[20], c;
	//double b[4] = { 1,2,3,4 };
	finput_matrix(a, 4, 5, loc1);
	printf("a = \n");
	output_matrix(a, 4, 5);
	printf("bool = %d\n", grammer(a, 4, b));
	printf("b = \n");
	output_matrix(b, 4, 1);
	solve_liner_equation(a, 4, 5, b);
	printf("liner_equation b = \n");
	output_matrix(b, 4, 1);
	submatrix(a, 4, 5, 0, 4, 0, 4, a);
	printf("a = \n");
	output_matrix(a, 4, 4);
	matrix_multiply(a, b, 4, 4, 1, b);
	printf("ab = \n");
	output_matrix(b, 1, 4);
	printf("laplace - guass = %lf\n", laplace_expension(a, 4) - determination_by_guass(a, 4));

#endif

#if 0
	// 设置随机种子
	srand((unsigned)time(NULL));

	// 初始化图形模式
	initgraph(640, 480);

	int  x, y;
	char c;

	settextstyle(16, 8, _T("Courier"));	// 设置字体

	// 设置颜色
	settextcolor(GREEN);
	setlinecolor(BLACK);

	for (int i = 0; i <= 479; i++)
	{
		// 在随机位置显示三个随机字母
		for (int j = 0; j < 3; j++)
		{
			x = (rand() % 80) * 8;
			y = (rand() % 20) * 24;
			c = (rand() % 26) + 65;
			outtextxy(x, y, c);
		}

		// 画线擦掉一个像素行
		line(0, i, 639, i);

		Sleep(10);					// 延时
		if (i >= 479)	i = -1;
		if (_kbhit())	break;		// 按任意键退出
	}

	// 关闭图形模式
	closegraph();

#endif
	return 0;
}

/*自定义函数*/
double fm(double * r)
{
	double x = r[0], y = r[1], z = r[2];
	if (x * x / 25 + y * y / 16 + z * z / 9 > 1)
		return 0;
	else return 1;
}

double fx(double x)
{
	return //x * sqrt(asin(0.00001 - cos(x * x * x * x)));
		16 * sin(x) * sin(x) * sin(x);
}

double gx(double x)
{
	return  //-x * sqrt(asin(0.00001 - cos(x * x * x * x)));
		13 * cos(x) - 5 * cos(2 * x) - 2 * cos(3 * x) - cos(4 * x);
}

double fxy(double x, double y)
{
	return //4 / sqrt(x * x + y * y) - 1 / sqrt(x * x + (y - 3) * (y - 3)) + 2 / sqrt((x - 2) * (x - 2) + y * y);
		//sqrt(1 - x * x / 9 - y * y / 25) * 4;
		//x * x - y * y;
		//pow(x * x + y * y - 1, 3) - x * x * y * y;
		//-10 / sqrt(x * x + y * y);
		sqrt(x * x + y * y + 1);
}

double gxy(double x, double y)
{
	return -sqrt(1 - x * x / 9 - y * y / 25) * 4;
}