#include"matrix caculator.h"

void line3d(double x0, double y0, double z0, double x1, double y1, double z1)
{
	double a[6] = { x0, x1, y0, y1, z0, z1 };
	plot3d.per(a, 2, a);
	line((int)a[0], (int)a[3], (int)a[1], (int)a[4]);
}

void cuboid(double x0, double y0, double z0, double x1, double y1, double z1)
{
	double a[24] = { x0, x1, x1, x0, x0, x1, x1, x0, y0, y0, y1, y1, y0, y0, y1, y1, z0, z0, z0, z0, z1, z1, z1, z1 };
	plot3d.per(a, 8, a);
	line((int)a[0], (int)a[8], (int)a[1], (int)a[9]);
	line((int)a[1], (int)a[9], (int)a[2], (int)a[10]);
	line((int)a[2], (int)a[10], (int)a[3], (int)a[11]);
	line((int)a[3], (int)a[11], (int)a[0], (int)a[8]);
	line((int)a[0], (int)a[8], (int)a[4], (int)a[12]);
	line((int)a[1], (int)a[9], (int)a[5], (int)a[13]);
	line((int)a[2], (int)a[10], (int)a[6], (int)a[14]);
	line((int)a[3], (int)a[11], (int)a[7], (int)a[15]);
	line((int)a[4], (int)a[12], (int)a[5], (int)a[13]);
	line((int)a[5], (int)a[13], (int)a[6], (int)a[14]);
	line((int)a[6], (int)a[14], (int)a[7], (int)a[15]);
	line((int)a[7], (int)a[15], (int)a[4], (int)a[12]);

}