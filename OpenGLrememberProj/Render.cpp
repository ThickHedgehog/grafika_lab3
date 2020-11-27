

#include "Render.h"
#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <math.h>
#include <cmath>
#include <array>
#include <vector>
#include <random>
#include <algorithm>
#define PI 3.14159265

double t_max = 0;
bool Loop = true;
int N = 3; // строки
int M = 3; // столбцы

// Контрольные точки поверхности
std::vector <std::vector <std::vector <double>>> surfacePoints = { { {0, 0, 0.1}, {0, 5, 7.5}, {0, 10, 9}, {0, 15, 2} },
																	{ {5, 0, 6.1}, {5, 5, 7.8}, {5, 10, 9.6}, {5, 15, 2} },
																	{ {10, 0, 2}, {10, 5, 0}, {10, 10, 9.6}, {10, 15, 2} },
																	{ {15, 0, 9}, {15, 5, 9.4}, {15, 10, 9.6}, {15, 15, 2} } };

// Формула Безье 2 порядка
double f(double p1, double p2, double p3, double t)
{
	return p1 * (1 - t) * (1 - t) + 2 * p2 * t * (1 - t) + p3 * t * t; //посчитанная формула
}

// Формула Безье 3 порядка
double f2(double p0, double p1, double p2, double p3, double t)
{
	return p0 * (1 - t) * (1 - t) * (1 - t) + 3 * t * p1 * (1 - t) * (1 - t) + 3 * t * t * p2 * (1 - t) + t * t * t * p3; //Безье 3го порядка
}

// Формула Эрмита
double fer(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (3 * t * t - 2 * t * t * t) + r1 * (t * t * t - 2 * t * t + t) + r4 * (t * t * t - t * t); //Эрмит
}

// Расчет координат вектора
void CalcVector(double x1[3], double x2[3], double x3[3], double x4[3], double *n1, double *n2) {
	n1[0] = x2[0] - x1[0];
	n1[1] = x2[1] - x1[1];
	n1[2] = x2[2] - x1[2];

	n2[0] = x4[0] - x3[0];
	n2[1] = x4[1] - x3[1];
	n2[2] = x4[2] - x3[2];
}

// Рисование фигуры
void DrawFigure(double* Point) {

	glPushMatrix();
	glTranslated(Point[0], Point[1], Point[2]);

	double A[] = { 1, -1, 0 };
	double B[] = { 1, 1, 0 };
	double C[] = { 1, 1, 1 };
	double D[] = { 1, -1, 1 };

	double A2[] = { 1, 1, 0 };
	double B2[] = { -1, 1, 0 };
	double C2[] = { -1, 1, 1 };
	double D2[] = { 1, 1, 1 };

	double A3[] = { -1, 1, 0 };
	double B3[] = { -1, 1, 1 };
	double C3[] = { -1, -1, 1 };
	double D3[] = { -1, -1, 0 };

	double A4[] = { 1, -1, 0 };
	double B4[] = { 1, -1, 1 };
	double C4[] = { -1, -1, 1 };
	double D4[] = { -1, -1, 0 };

	double A5[] = { 1, 1, 1 };
	double B5[] = { 0, 0, 3 };
	double C5[] = { -1, 1, 1 };

	double A6[] = { -1, 1, 1 };
	double B6[] = { 0, 0, 3 };
	double C6[] = { -1, -1, 1 };

	double A7[] = { -1, -1, 1 };
	double B7[] = { 0, 0, 3 };
	double C7[] = { 1, -1, 1 };

	double A8[] = { 1, -1, 1 };
	double B8[] = { 0, 0, 3 };
	double C8[] = { 1, 1, 1 };

	double A9[] = { 1, 1, 0 };
	double B9[] = { 1, -1, 0 };
	double C9[] = { -1, -1, 0 };
	double D9[] = { -1, 1, 0 };

	glBegin(GL_QUADS);
	glColor3d(1, 0, 0);
	glVertex3dv(A);
	glVertex3dv(B);
	glVertex3dv(C);
	glVertex3dv(D);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1, 0, 1);
	glVertex3dv(A2);
	glVertex3dv(B2);
	glVertex3dv(C2);
	glVertex3dv(D2);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1, 0, 0);
	glVertex3dv(A3);
	glVertex3dv(B3);
	glVertex3dv(C3);
	glVertex3dv(D3);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1, 0, 1);
	glVertex3dv(A4);
	glVertex3dv(B4);
	glVertex3dv(C4);
	glVertex3dv(D4);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0, 1, 0);
	glVertex3dv(A5);
	glVertex3dv(B5);
	glVertex3dv(C5);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0, 1, 1);
	glVertex3dv(A6);
	glVertex3dv(B6);
	glVertex3dv(C6);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0, 1, 0);
	glVertex3dv(A7);
	glVertex3dv(B7);
	glVertex3dv(C7);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0, 1, 1);
	glVertex3dv(A8);
	glVertex3dv(B8);
	glVertex3dv(C8);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(0, 0, 1);
	glVertex3dv(A9);
	glVertex3dv(B9);
	glVertex3dv(C9);
	glVertex3dv(D9);
	glEnd();

	glPopMatrix();
}

// Рисование фигуры с поворотом для линии Безье 3 порядка
void DrawFigureRotate(double* Point, double* P1, double* P2, double* P3, double* P4, double t) {

	double PDer[6], PDelta[3];

	glPushMatrix();
	glTranslated(Point[0], Point[1], Point[2]);
	PDer[0] = f2(P1[0], P2[0], P3[0], P4[0], t);
	PDer[1] = f2(P1[1], P2[1], P3[1], P4[1], t);
	PDer[2] = f2(P1[2], P2[2], P3[2], P4[2], t);

	PDer[3] = f2(P1[0], P2[0], P3[0], P4[0], t + 0.01);
	PDer[4] = f2(P1[1], P2[1], P3[1], P4[1], t + 0.01);
	PDer[5] = f2(P1[2], P2[2], P3[2], P4[2], t + 0.01);

	PDelta[0] = PDer[3] - PDer[0];
	PDelta[1] = PDer[4] - PDer[1];
	PDelta[2] = PDer[5] - PDer[2];

	double Orig[3] = { 1, 0, 0 };
	double RotateX[3] = { PDelta[0], PDelta[1], 0 };
	// normalize
	double length = sqrt(pow(RotateX[0], 2) + pow(RotateX[1], 2) + pow(RotateX[2], 2));
	RotateX[0] /= length;
	RotateX[1] /= length;
	RotateX[2] /= length;

	double cos = RotateX[0] * Orig[0] + RotateX[1] * Orig[1] + RotateX[2] * Orig[2];
	double VecProizv[3] = {
		(Orig[2] * RotateX[1] - Orig[1] * RotateX[2]),
		(Orig[0] * RotateX[2] - Orig[2] * RotateX[0]),
		(Orig[1] * RotateX[0] - Orig[0] * RotateX[1]) };

	double Sign = VecProizv[2] / abs(VecProizv[2]);
	double Angle1 = -acos(cos) * 180.0 / PI * Sign;
	double Angle2 = acos(PDelta[2] / sqrt(pow(PDelta[0], 2) + pow(PDelta[1], 2) + pow(PDelta[2], 2))) * 180.0 / PI - 90.0;

	glRotated(Angle1, 0, 0, 1);	// x
	glRotated(Angle2, 0, 1, 0);	// y

	double A[] = { 1, -1, -1 };
	double B[] = { 1, 1, -1 };
	double C[] = { 1, 1, 1 };
	double D[] = { 1, -1, 1 };

	double A2[] = { 1, 1, -1 };
	double B2[] = { -1, 1, -1 };
	double C2[] = { -1, 1, 1 };
	double D2[] = { 1, 1, 1 };

	double A3[] = { -1, 1, -1 };
	double B3[] = { -1, 1, 1 };
	double C3[] = { -1, -1, 1 };
	double D3[] = { -1, -1, -1 };

	double A4[] = { 1, -1, -1 };
	double B4[] = { 1, -1, 1 };
	double C4[] = { -1, -1, 1 };
	double D4[] = { -1, -1, -1 };

	double A5[] = { 1, 1, 1 };
	double B5[] = { 0, 0, 3 };
	double C5[] = { -1, 1, 1 };

	double A6[] = { -1, 1, 1 };
	double B6[] = { 0, 0, 3 };
	double C6[] = { -1, -1, 1 };

	double A7[] = { -1, -1, 1 };
	double B7[] = { 0, 0, 3 };
	double C7[] = { 1, -1, 1 };

	double A8[] = { 1, -1, 1 };
	double B8[] = { 0, 0, 3 };
	double C8[] = { 1, 1, 1 };

	double A9[] = { 1, 1, -1 };
	double B9[] = { 1, -1, -1 };
	double C9[] = { -1, -1, -1 };
	double D9[] = { -1, 1, -1 };

	double A10[] = { 1, 1, 1 };
	double B10[] = { 1, -1, 1 };
	double C10[] = { -1, -1, 1 };
	double D10[] = { -1, 1, 1 };

	glBegin(GL_QUADS);
	glColor3d(1, 0, 0);
	glVertex3dv(A);
	glVertex3dv(B);
	glVertex3dv(C);
	glVertex3dv(D);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1, 0, 1);
	glVertex3dv(A2);
	glVertex3dv(B2);
	glVertex3dv(C2);
	glVertex3dv(D2);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1, 0, 0);
	glVertex3dv(A3);
	glVertex3dv(B3);
	glVertex3dv(C3);
	glVertex3dv(D3);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1, 0, 1);
	glVertex3dv(A4);
	glVertex3dv(B4);
	glVertex3dv(C4);
	glVertex3dv(D4);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(0, 0, 1);
	glVertex3dv(A9);
	glVertex3dv(B9);
	glVertex3dv(C9);
	glVertex3dv(D9);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(0, 0, 1);
	glVertex3dv(A10);
	glVertex3dv(B10);
	glVertex3dv(C10);
	glVertex3dv(D10);
	glEnd();
	
	glPopMatrix();
}

// Рисование линии Безье 2 порядка
void DrawLineBez2(double* P1, double* P2, double* P3) {
	double P[3];
	glBegin(GL_LINES); //построим отрезки P1P2 и P2P3
	glVertex3dv(P1);
	glVertex3dv(P2);
	glEnd();
	glBegin(GL_LINES);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glEnd();
	glLineWidth(3); //ширина линии
	glColor3d(0, 1, 0);
	// рисуем линию
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		P[0] = f(P1[0], P2[0], P3[0], t);
		P[1] = f(P1[1], P2[1], P3[1], t);
		P[2] = f(P1[2], P2[2], P3[2], t);
		glVertex3dv(P); //Рисуем точку P
	}
	glEnd();
	glColor3d(1, 0, 1);
	glLineWidth(1); //возвращаем ширину линии = 1
	// рисуем пирамидку
	for (double t = 0; t <= t_max; t += 0.01)
	{
		P[0] = f(P1[0], P2[0], P3[0], t);
		P[1] = f(P1[1], P2[1], P3[1], t);
		P[2] = f(P1[2], P2[2], P3[2], t);
	}
	DrawFigure(P);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	glVertex3dv(P1);
	glVertex3dv(P3);
	glEnd();
}

// Рисование линии Безье 3 порядка
void DrawLineBez3(double* P1, double* P2, double* P3, double* P4) {
	double P_2[3];

	glBegin(GL_LINES); //построим отрезки P1P2 и P2P3 и Р3Р4
	glVertex3dv(P1);
	glVertex3dv(P2);
	glEnd();
	glBegin(GL_LINES);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glEnd();
	glBegin(GL_LINES);
	glVertex3dv(P3);
	glVertex3dv(P4);
	glEnd();
	glLineWidth(3); //ширина линии
	glColor3d(0, 1, 0);
	// рисуем линию
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		P_2[0] = f2(P1[0], P2[0], P3[0], P4[0], t);
		P_2[1] = f2(P1[1], P2[1], P3[1], P4[1], t);
		P_2[2] = f2(P1[2], P2[2], P3[2], P4[2], t);
		glVertex3dv(P_2);
		//Рисуем точку P
	}
	glEnd();
	glColor3d(1, 0, 1);
	glLineWidth(1); //возвращаем ширину линии = 1
	for (double t = 0; t <= t_max; t += 0.01)
	{
		P_2[0] = f2(P1[0], P2[0], P3[0], P4[0], t);
		P_2[1] = f2(P1[1], P2[1], P3[1], P4[1], t);
		P_2[2] = f2(P1[2], P2[2], P3[2], P4[2], t);
	}
	// рисуем пирамидку
	double t = t_max;
	DrawFigureRotate(P_2, P1, P2, P3, P4, t);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	glVertex3dv(P1);
	glVertex3dv(P4);
	glEnd();
}

// Рисование линии Эрмита
void DrawLineErmit(double* Per1, double* Per2, double* Per3, double* Per4) {
	double Rer1[3], Rer2[3];
	CalcVector(Per1, Per2, Per3, Per4, Rer1, Rer2);
	double Per[4];
	glColor3d(0, 0, 1);
	glBegin(GL_LINES); //построим отрезки P1P2 и P3P4
	glVertex3dv(Per1);
	glVertex3dv(Per2);
	glEnd();
	glBegin(GL_LINES);
	glVertex3dv(Per3);
	glVertex3dv(Per4);
	glEnd();
	glLineWidth(3); //ширина линии
	glColor3d(0, 1, 0);
	// рисуем линию
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		Per[0] = fer(Per1[0], Per3[0], Rer1[0], Rer2[0], t);
		Per[1] = fer(Per1[1], Per3[1], Rer1[1], Rer2[1], t);
		Per[2] = fer(Per1[2], Per3[2], Rer1[2], Rer2[2], t);
		glVertex3dv(Per); //Рисуем точку P
	}
	glEnd();
	glColor3d(1, 0, 1);
	glLineWidth(1); //возвращаем ширину линии = 1
	// рисуем пирамидку
	for (double t = 0; t <= t_max; t += 0.01)
	{
		Per[0] = fer(Per1[0], Per3[0], Rer1[0], Rer2[0], t);
		Per[1] = fer(Per1[1], Per3[1], Rer1[1], Rer2[1], t);
		Per[2] = fer(Per1[2], Per3[2], Rer1[2], Rer2[2], t);
	}
	DrawFigure(Per);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	glVertex3dv(Per1);
	glVertex3dv(Per2);
	glVertex3dv(Per3);
	glVertex3dv(Per4);
	glEnd();
}

// Расчет нормали
float* FindNormal(float d1[], float d2[], float d3[])
{
	float v1[] = { d1[0] - d2[0], d1[1] - d2[1], d1[2] - d2[2] };
	float v2[] = { d3[0] - d2[0], d3[1] - d2[1], d3[2] - d2[2] };
	float normal[3];

	normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
	normal[1] = -1 * (v1[0] * v2[2] - v1[2] * v2[0]);
	normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

	return normal;
}

// Расчет факториала
float Fact(int n)
{
	if (n == 0) return 1;
	else return n * Fact(n - 1);
}

// Расчет многочлена Берштейна
double calculateBernsteinPolynomials(int i, int n, double u) {
	n -= 1;
	return (Fact(n) / (Fact(i) * Fact(n - i))) * pow(u, i) * pow(1.0 - u, (double)n - i);
}

// Расчет поверхности Безье
void CalcSurface(double* point, double u, double v) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			double bernsteinPolynomials_N = calculateBernsteinPolynomials(i, N, u);
			double bernsteinPolynomials_M = calculateBernsteinPolynomials(j, M, v);

			point[0] += bernsteinPolynomials_N * bernsteinPolynomials_M * surfacePoints[i][j][0];
			point[1] += bernsteinPolynomials_N * bernsteinPolynomials_M * surfacePoints[i][j][1];
			point[2] += bernsteinPolynomials_N * bernsteinPolynomials_M * surfacePoints[i][j][2];
		}
}

// Рисование поверхности Безье из точек
void DrawSurfacePoints() {
	glPointSize(5);

	glBegin(GL_POINTS);
	for (double u = 0; u < 1.0; u += 0.1) {
		for (double v = 0; v < 1.0; v += 0.1) {
			double point[3] = { 0, 0, 0 };
			CalcSurface(point, u, v);
			glVertex3dv(point);
		}	
	}
	glEnd();
}

// Рисование поверхности Безье из линий
void DrawSurfaceLines() {
	glColor3d(0, 0, 0);

	for (double u = 0; u < 1.0; u += 0.1) {
		glBegin(GL_LINE_STRIP);
		for (double v = 0; v < 1.0; v += 0.1) {
			double point[3] = { 0, 0, 0 };
			CalcSurface(point, u, v);
			glVertex3dv(point);
		}
		glEnd();
	}
	for (double u = 0; u < 1.0; u += 0.1) {
		glBegin(GL_LINE_STRIP);
		for (double v = 0; v < 1.0; v += 0.1) {
			double point[3] = { 0, 0, 0 };
			CalcSurface(point, v, u);
			glVertex3dv(point);
		}
		glEnd();
	}
}

void Render(double delta_time)
{    
	 //t_max становится = 1 за 5 секунд
	if (t_max < 1 && Loop) 
	{
		t_max += delta_time / 5;
	}
	else
	{
		Loop = false;
		t_max -= delta_time / 5;
		if (t_max <= 0) {
			Loop = true;
		}
	}//после обнуляется

	double P1[] = { 0,0,0 }; //Наши точки
	double P2[] = { -4,6,7 };
	double P3[] = { 10,10,0 };
	/*DrawLineBez2(P1, P2, P3);*/
	//
	//2 линия
	//
	double P1_1[] = { 0,0,0 }; //Наши точки
	double P2_2[] = { 7,7,-20 };
	double P3_3[] = { 20,20,10 };
	double P4_4[] = { 0,0,0 };

	//double P1_1[] = { 0,0,0 }; //Наши точки
	//double P2_2[] = { 1,7,3 };
	//double P3_3[] = { 12,1,3 };
	//double P4_4[] = { 5,5,0 };
	/*DrawLineBez3(P1_1, P2_2, P3_3, P4_4);*/
	//
	//3 линия (Эрмита)
	//
	double Per1[] = { 0,0,0 }; //Наши точки
	double Per2[] = { 8,3,0 };
	double Per3[] = { 3,5,0 };
	double Per4[] = { 8,0,0 };
	/*DrawLineErmit(Per1, Per2, Per3, Per4);*/
	//
	//4 линия (Эрмита)
	//
	double Per1_1[] = { 0,0,0 }; //Наши точки
	double Per2_2[] = { 10,15,8 };
	double Per3_3[] = { 14,10,10 };
	double Per4_4[] = { 7,8,0 };
	/*DrawLineErmit(Per1_1, Per2_2, Per3_3, Per4_4);*/

	DrawSurfacePoints();

	DrawSurfaceLines();
}   