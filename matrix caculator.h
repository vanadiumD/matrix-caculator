#ifndef MATRIX_CACULATOR_H
#define MATRIX_CACULATOR_H

#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include<stdlib.h>
#include<graphics.h>
#include<conio.h>
#include<time.h>
#define M 30
#define N 30
#define R 10
#define ACCURACY 0.0001
#define PLOTSTEP 0.1
#define MAX_POINT 1000000

typedef struct node
{
	double a[N];
	struct node *next;
}NODE;


struct RANGE
{
	double x0;
	double x1;
	double y0;
	double y1;
	double z0;
	double z1;
	double t0;
	double t1;
};
extern struct RANGE range[R];

struct PLOT
{
	int H;
	int L;
	int padding;
	int line_thickness;
	int text_size;
	int point_size;
	float x_scale;
	float y_scale;
	COLORREF bk_color;
	COLORREF xoy_color;
	COLORREF indicator_color;
	COLORREF scale_color;
	COLORREF point_color;
	COLORREF vector_color;
	COLORREF function_color;
	bool coordinate;
};
extern struct PLOT plot;

struct OTMP
{
	double scale;
};
extern struct OTMP otm;

struct PERSPECTIVE_1
{
	double u;//坐标原点到光源距离
	double v;//坐标原点到近裁剪平面距离;
};
extern PERSPECTIVE_1 per1;

struct PLOT3D_ANGLE
{
	double theta;
	double psai;
	double phi;
};

struct PLOT3D
{
	int H;
	int L;
	int padding;
	int line_thickness;
	int text_size;
	int point_size;
	int grid_spacing;
	double x_length;
	double y_length;
	double z_length;
	double x_scale;
	double y_scale;
	double z_scale;
	void (*per)(double *a, int n, double *res);
	struct PLOT3D_ANGLE angle3d;
	bool udf_color;
	COLORREF bk_color;
	COLORREF xoy_color;
	COLORREF indicator_color;
	COLORREF scale_color;
	COLORREF point_color;
	COLORREF vector_color;
	COLORREF function_color;
	COLORREF makerline_color;
	bool initplot3d;
};
extern struct PLOT3D plot3d;

struct PLOT3D_DYNAMIC
{
	double theta;
	double psai;
	double phi;
	double otmp_scale;
	double per1_v;
};
extern struct PLOT3D_DYNAMIC plot3d_dym;

/*自定义函数*/

double fm(double * r);

double fx(double x);

double gx(double x);

double fxy(double x, double y);

double gxy(double x, double y);

double udf_parametric_function(double t, double * res);

double udf_scalar_field(double * r, double t);

double udf_vector_field(double * r, double t, double *res);

/*基本函数*/
/*m行数, n列数, a, b操作矩阵, res:储存结果的矩阵*/

bool initmatrix(double *a, int m, int n);

bool identity_matrix(double *a, int n);

void input_matrix(double *a, int m, int n);

void matrix_copy(double *a, int m, int n, double * res);

void column_copy(double *a, int m, int n, int column, double *res, int s, int rescolumn);

bool submatrix(double *a, int m, int n, int rstart, int rlength, int cstart, int clength, double *res);

bool finput_matrix(double *a, int m, int n, char loc[]);

void output_matrix(double *a, int m, int n);

void matrix_sum(double *a, double *b, int m, int n, double *res);

void matrix_number_multiply(double *a, int m, int n, int times, double *res);

void matrix_multiply(double *a, double *b, int m, int s, int n, double *res);

void matrix_power(double *a, int n, int power, double *res);

int row_reduce(double *a, int m, int n, double *res, char simp, int model);

int column_reduce(double *a, int m, int n, double *res, char simp, int model);

int offset_standard_by_guass(double *a, int m, int n, double *res);

int offset_standard(double *a, int m, int n, double *res);

int solve_liner_equation(double *a, int m, int n, double * res);

int least_square_method(double *a, int m, int n, double *res);

void output_solution(double *res, int m, int n);

bool grammer(double *a, int n, double *res);

double laplace(double *a, double *p1, int n, int b, int lastreduce);

double laplace_expension(double *a, int n);

double determination_by_guass(double *a, int n);

bool inverse_by_guass(double *a, int n, double *res);

int rank(double *a, int m, int n);

void transform(double *a, int m, int n, double *res);

void gh_decomposition(double *a, int m, int n, double *g, double *h);//unfinished

bool qr_decomposition(double *a, int n, double *q, double *r);

bool judge_offset(double *a, double *b, int m, int n);

int jacobian_method(double *a, int n, double accuracy, double *res);

int characteristic_value(double *a, int n, double *res);

/*向量*/
double dot_product(double *vec1, double *vec2, int n);

double cross_product(double *vec1, double *vec2, double *res);

double vector_angle(double *vec1, double *vec2, int n);

double vector_module(double *vec, int n);

//向量投影//将vec1投影到vec2上
double vector_projection(double *vec1, double *vec2, int n);

void printlist(NODE *head);

double mochang(double *a);

void group_sort(NODE * head);

void inner_sort(NODE * head);

NODE *rear_creat(void);

NODE *head_creat(void);

NODE *positioninsert(NODE *head);

NODE *orderinsert(NODE * head);

/*应用*/
double fibonacci_sequence(int n);

double rotate(double *a, int n, double deg, double *res);

double total_mass(double(*pm)(double *r), struct RANGE range);

void centroid(double m, double(*pm)(double *r), struct RANGE range, double *res);

void inertia_tensor(double *cen,  double(*pm)(double *r), struct RANGE range, double *res);

double x_rotate(double *a, int n, double deg, double *res);

double y_rotate(double *a, int n, double deg, double *res);

double z_rotate(double *a, int n, double deg, double *res);

double rotate3D(double *a, int n, char pos1, double deg1, char pos2, double deg2, char pos3, double deg3, double *res);

/*可视化*/
/*2维*/
void plot_coordinate(void);

void plot_vector(double *a, int n);

void plot_point(double *a, int n);

void plot_function(double(*pf)(double x), double x0, double x1);

void plot_parametric_function(double(*pf)(double t), double(*pg)(double t), double t0, double t1);

void plot_implicit_function(double(*pf)(double x, double y), double x0, double x1, double y0, double y1);
/*三维*/
void initplot3D(void);

//一点透视
void perspective_1(double *a, int n, double *res);

//斜二测画法
void otmp(double *a, int n, double *res);

void plot3D_coordinate(void);

void plot3D_makerline(double *a, int n);

void plot3D_vector(double *a, int n);

void plot3D_point(double *a, int n);

void plot3D_binary_function(double(*pf)(double x, double y), double x0, double x1, double y0, double y1);

void plot3D_render_binary_function(double(*pf)(double x, double y), double x0, double x1, double y0, double y1);

void plot3D_binary_function_onxyplane(double(*pf)(double x, double y), double x0, double x1, double y0, double y1);

int plot3d_dynamic(void);

/*graphics3d*/

void line3d(double x0, double y0, double z0, double x1, double y1, double z1);

void cuboid(double x0, double y0, double z0, double x1, double y1, double z1);


#endif