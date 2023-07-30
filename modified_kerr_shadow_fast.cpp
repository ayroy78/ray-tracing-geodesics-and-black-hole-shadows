#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>

// these functions calculate metric values
double g_tt(double r, double rh, double th, double f0, double f2, double w)
{
	return -(1.0 - rh / r) * exp(2.0 * f0) + exp(2.0 * f2) * pow(sin(th) * w, 2.0);
}

double g_tphi(double r, double th, double f2, double w)
{
	return -w * r * exp(2.0 * f2) * pow(sin(th), 2.0);
}

double g_rr(double r, double rh, double f1)
{
	return exp(2.0 * f1) / (1.0 - rh / r);
}

double g_thth(double r, double f1)
{
	return exp(2.0 * f1) * pow(r, 2.0);
}

double g_phiphi(double r, double th, double f2)
{
	return exp(2.0 * f2) * pow(r * sin(th), 2.0);
}

// functions for derivatives for radial equation of motion
double g_rr_th(double r, double rh, double f1, double f1_th)
{
	return (r / (r - rh)) * exp(2.0 * f1) * 2.0 * f1_th;
}

double g_rr_r(double r, double rh, double f1, double f1_r)
{
	return -exp(2.0 * f1) * rh / pow((r - rh), 2.0) + (r / (r - rh)) * exp(2.0 * f1) * 2.0 * f1_r;
}

double g_tt_r(double r, double rh, double th, double f0, double f2, double w, double f0_r, double f2_r, double w_r)
{
	return -exp(2.0 * f0) * rh / pow(r, 2.0) - (1.0 - rh / r) * exp(2.0 * f0) * 2.0 * f0_r + pow(sin(th), 2) * (2.0 * w * w_r * exp(2.0 * f2) + w * w * 2.0 * exp(2.0 * f2) * f2_r);
}

double g_tphi_r(double r, double th, double f2, double w, double f2_r, double w_r)
{
	return -w * pow(sin(th), 2.0) * exp(2.0 * f2) - w_r * r * pow(sin(th), 2) * exp(2.0 * f2) - w * r * 2.0 * pow(sin(th), 2) * exp(2.0 * f2) * f2_r;
}

double g_phiphi_r(double r, double th, double f2, double f2_r)
{
	return 2.0 * r * pow(sin(th), 2.0) * exp(2.0 * f2) + pow(r * sin(th), 2) * 2.0 * exp(2.0 * f2) * f2_r;
}

double g_thth_r(double r, double f1, double f1_r)
{
	return 2.0 * r * exp(2.0 * f1) + r * r * 2.0 * exp(2.0 * f1) * f1_r;
}

// derivatives for the theta equation of motion

double g_tt_th(double r, double rh, double th, double f0, double f2, double w, double f0_th, double f2_th, double w_th)
{
	return -(1.0 - rh / r) * exp(2.0 * f0) * f0_th + exp(2.0 * f2) * pow(w, 2.0) * sin(2.0 * th) + pow(sin(th), 2) * (2.0 * w * exp(2.0 * f2) * w_th + 2.0 * w * w * exp(2.0 * f2) * f2_th);
}

double g_tphi_th(double r, double th, double f2, double w, double f2_th, double w_th)
{
	return -w_th * r * exp(2.0 * f2) * pow(sin(th), 2) - w * r * exp(2.0 * f2) * sin(2.0 * th) - w * r * 2.0 * exp(2.0 * f2) * pow(sin(th), 2) * f2_th;
}

double g_thth_th(double r, double f1, double f1_th)
{
	return r * r * exp(2.0 * f1) * 2.0 * f1_th;
}

double g_phiphi_th(double r, double th, double f2, double f2_th)
{
	return exp(2.0 * f2) * pow(r, 2.0) * sin(2.0 * th) + pow(r * sin(th), 2) * 2.0 * exp(2.0 * f2) * f2_th;
}

// this is a function to create linearly separated points.
void linspace(double start, double end, double numpoints, double *points)
{
	// double * points  = new double [numpoints];
	double dx = (end - start) / (numpoints - 1);
	for (int i = 0; i < numpoints; i++)
	{
		points[i] = (start + i * dx);
	}
}

class Stationary_axisymmetric_data
{
public:
	double **data = new double *[6];
	double ***coeff = new double **[4];
	int xsize;
	int ysize;
	double rh;
	double thm;

	Stationary_axisymmetric_data(int x_size, int y_size, double r_h)
	{
		xsize = x_size;
		ysize = y_size;
		rh = r_h;
		thm =0.0;

		for (int i = 0; i < 6; i++)
		{
			data[i] = new double[xsize * ysize];
		}

		for (int i = 0; i < 4; i++)
		{
			coeff[i] = new double *[xsize * ysize];
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < xsize * ysize; j++)
			{
				coeff[i][j] = new double[16];
			}
		}
	}

	void make_data()
	{
		int counter = 0;
		std::ifstream myFile("kerr_rh=0.0359508947168857_omh=0.956.txt");
		double x;
		double y;
		double f0;
		double f1;
		double f2;
		double w;
		double phi;

		while (!myFile.eof())
		{
			myFile >> x >> y >> f0 >> f1 >> f2 >> w >> phi;
			//myFile >> x >> y >> f1 >> f2 >> f0 >> W;

			data[0][counter] = x;
			data[1][counter] = y;
			data[2][counter] = f0;
			data[3][counter] = f1;
			data[4][counter] = f2;
			data[5][counter] = -w;

			counter++;
		}

		myFile.close();
		thm = data[1][xsize*ysize-1];
		
	}

	int sorter(double x, double y)
	{	
		int j = 0;
		int i = 0;
		int icount = 0;
		double diffx = 0;
		double diffy = 0;
		double yd;
		double xd;

		if (M_PI_2 <= y && y <=M_PI){
			yd = M_PI - y;
		}

		else if (y<= 0 && y>= -M_PI_2){
			yd = fabs(y);
		}
		else if (y>=M_PI && y<= 1.5*M_PI){
			yd = y - M_PI;
		}

		else if (y<=-M_PI_2 && y>=-M_PI){
			yd = fabs(-M_PI -y);
		}

		else if (y>= 1.5*M_PI && y<= 2.0*M_PI){
			yd = 2.0*M_PI - y;
		}

		else if (y<=-M_PI && y>= -1.5*M_PI){
			yd = fabs(y +M_PI);
		}

		else if (y>=2.0*M_PI && y<= 2.5*M_PI){
			yd = y - 2.0*M_PI;
		}
		
		else {
			yd = y;
		}

		if ((y>=thm && y<=M_PI_2) || (yd>=thm && yd<=M_PI_2)){
			yd = thm-pow(10,-6);
		}

		xd = sqrt(x*x - rh*rh)/(sqrt(x*x - rh*rh) + 1);
		
		while (diffx>=0){
			diffx = xd - data[0][j];
			j = j + 1;
		}
		

		while (diffy>=0){
			diffy = yd - data[1][i];
			i = i + xsize;
			icount = icount + 1;
		}
		
		return (icount-2)*xsize + (j-2);

		
	}

	double diff_x(int index, int index2)
	{ // index is the point in data array, index2 is the interpolation term

		double f2;
		double f1;
		double f0;
		double h1;
		double h2;

		if (index % xsize == 0)
		{	
			f0 = data[index2][index];
			f1 = data[index2][index+1];
			f2 = data[index2][index+2];
			h1 = data[0][index+1] - data[0][index];
			h2 = data[0][index+2] - data[0][index];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h2 - h1));
		}
		else if ((index + 1) % xsize == 0)
		{
			f0 = data[index2][index];
			f1 = data[index2][index-1];
			f2 = data[index2][index-2];
			h1 = data[0][index] - data[0][index-1];
			h2 = data[0][index] - data[0][index-2];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h1 - h2));
		}
		else
		{	
			f0 = data[index2][index];
			f1 = data[index2][index+1];
			f2 = data[index2][index-1];
			h1 = data[0][index+1] - data[0][index];
			h2 = data[0][index] - data[0][index-1];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h2 + h1));
		}
	}

	double diff_y(int index, int index2)
	{ // index is the point in data array, index2 is the metric term

		double f2;
		double f1;
		double f0;
		double h1;
		double h2;

		if (0 <= index && index <= xsize - 1)
		{	
			f0 = data[index2][index];
			f1 = data[index2][index+xsize];
			f2 = data[index2][index+2*xsize];
			h1 = data[1][index+xsize] - data[1][index];
			h2 = data[1][index+2*xsize] - data[1][index];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h2 - h1));
		}
		else if (xsize * ysize - xsize <= index && index <= xsize * ysize - 1)
		{	
			f0 = data[index2][index];
			f1 = data[index2][index-xsize];
			f2 = data[index2][index-2*xsize];
			h1 = data[1][index] - data[1][index-xsize];
			h2 = data[1][index] - data[1][index-2*xsize];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h1 - h2));
		}
		else
		{
			f0 = data[index2][index];
			f1 = data[index2][index+xsize];
			f2 = data[index2][index-xsize];
			h1 = data[1][index+xsize] - data[1][index];
			h2 = data[1][index] - data[1][index-xsize];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h2 + h1));
		}
	}

	double diff_xy(int index, int index2)
	{
		double f2;
		double f1;
		double f0;
		double h1;
		double h2;

		if ( 0<= index && index <=xsize -1)
		{	
			f0 = diff_x(index, index2);
			f1 = diff_x(index + xsize, index2);
			f2 = diff_x(index + 2*xsize, index2);
			h1 = data[1][index+xsize] - data[1][index];
			h2 = data[1][index +2*xsize] - data[1][index];

			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h2 - h1));
		}
		
		else if (xsize *ysize - xsize <= index && index <= xsize*ysize -1){
			f0 = diff_x(index, index2);
			f1 = diff_x(index - xsize, index2);
			f2 = diff_x(index - 2*xsize, index2);
			h1 = data[1][index] - data[1][index - xsize];
			h2 = data[1][index] - data[1][index - 2*xsize];

			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h1 - h2));
		}

		else
		{	
			f0 = diff_x(index, index2);
			f1 = diff_x(index + xsize, index2);
			f2 = diff_x(index - xsize, index2);
			h1 = data[1][index+xsize] - data[1][index];
			h2 = data[1][index] - data[1][index-xsize];
			return (f1*pow(h2,2) - f2*pow(h1,2) - f0*(pow(h2,2) - pow(h1,2)))/(h1*h2*(h2 + h1));
		}
	}

	void bicubic_spline_coefficient()
	{
		int index;
		int lower_right_index;
		int upper_left_index;
		int upper_right_index;
		double dx;
		double dy;

		for (int i = 0; i < xsize * ysize; i++)
		{
			if ((i + 1) % xsize == 0 || (xsize * ysize - xsize <= i && i <= xsize * ysize - 1))
			{
				continue;
			}
			// this is different from the grid visualization chosen for the spline interpolation with analytical kerr
			index = i;							   //(0,0) point
			lower_right_index = index + 1;		   // step in x
			upper_left_index = index + xsize;	   // step in y
			upper_right_index = index + xsize + 1; // step in x and y
			dx = data[0][index+1] - data[0][index];
			dy = data[1][index+ xsize] - data[1][index];

			for (int j = 2; j < 6; j++)
			{

				double c00 = data[j][index];
				double c01 = data[j][upper_left_index];
				double c02 = dy * diff_y(index, j);
				double c03 = dy * diff_y(upper_left_index, j);
				double c10 = data[j][lower_right_index];
				double c11 = data[j][upper_right_index];
				double c12 = dy * diff_y(lower_right_index, j);
				double c13 = dy * diff_y(upper_right_index, j);
				double c20 = dx * diff_x(index, j);
				double c21 = dx * diff_x(upper_left_index, j);
				double c22 = dx * dy * diff_xy(index, j);
				double c23 = dx * dy * diff_xy(upper_left_index, j);
				double c30 = dx * diff_x(lower_right_index, j);
				double c31 = dx * diff_x(upper_right_index, j);
				double c32 = dx * dy * diff_xy(lower_right_index, j);
				double c33 = dx * dy * diff_xy(upper_right_index, j);

				double a00 = c00;
				double a01 = c02;
				double a02 = -3.0 * c00 + 3.0 * c01 - 2.0 * c02 - c03;
				double a03 = 2.0 * c00 - 2.0 * c01 + c02 + c03;
				double a10 = c20;
				double a11 = c22;
				double a12 = -3.0 * c20 + 3.0 * c21 - 2.0 * c22 - c23;
				double a13 = 2.0 * c20 - 2.0 * c21 + c22 + c23;
				double a20 = -3.0 * c00 + 3.0 * c10 - 2.0 * c20 - c30;
				double a21 = -3.0 * c02 + 3.0 * c12 - 2.0 * c22 - c32;
				double a22 = 9.0 * c00 - 9.0 * c01 + 6.0 * c02 + 3.0 * c03 - 9.0 * c10 + 9.0 * c11 - 6.0 * c12 - 3.0 * c13 + 6.0 * c20 - 6.0 * c21 + 4.0 * c22 + 2.0 * c23 + 3.0 * c30 - 3.0 * c31 + 2.0 * c32 + c33;
				double a23 = -6.0 * c00 + 6.0 * c01 - 3.0 * c02 - 3.0 * c03 + 6.0 * c10 - 6.0 * c11 + 3.0 * c12 + 3.0 * c13 - 4.0 * c20 + 4.0 * c21 - 2.0 * c22 - 2.0 * c23 - 2.0 * c30 + 2.0 * c31 - c32 - c33;
				double a30 = 2.0 * c00 - 2.0 * c10 + c20 + c30;
				double a31 = 2.0 * c02 - 2.0 * c12 + c22 + c32;
				double a32 = -6.0 * c00 + 6.0 * c01 - 4.0 * c02 - 2.0 * c03 + 6.0 * c10 - 6.0 * c11 + 4.0 * c12 + 2.0 * c13 - 3.0 * c20 + 3.0 * c21 - 2.0 * c22 - c23 - 3.0 * c30 + 3.0 * c31 - 2.0 * c32 - c33;
				double a33 = 4.0 * c00 - 4.0 * c01 + 2.0 * c02 + 2.0 * c03 - 4.0 * c10 + 4.0 * c11 - 2.0 * c12 - 2.0 * c13 + 2.0 * c20 - 2.0 * c21 + c22 + c23 + 2.0 * c30 - 2.0 * c31 + c32 + c33;

				coeff[j - 2][i][0] = a00;
				coeff[j - 2][i][1] = a01;
				coeff[j - 2][i][2] = a02;
				coeff[j - 2][i][3] = a03;
				coeff[j - 2][i][4] = a10;
				coeff[j - 2][i][5] = a11;
				coeff[j - 2][i][6] = a12;
				coeff[j - 2][i][7] = a13;
				coeff[j - 2][i][8] = a20;
				coeff[j - 2][i][9] = a21;
				coeff[j - 2][i][10] = a22;
				coeff[j - 2][i][11] = a23;
				coeff[j - 2][i][12] = a30;
				coeff[j - 2][i][13] = a31;
				coeff[j - 2][i][14] = a32;
				coeff[j - 2][i][15] = a33;
			}
		}
	}

	double bicubic_spline(double x, double y, int index, int index2)
	{	
		double dx = data[0][index+1] - data[0][index];
		double dy = data[1][index+ xsize] - data[1][index]; 
		double yd;
		double xd;
		double ybar;
		double xbar;

		if (M_PI_2 <= y && y <=M_PI){
			yd = M_PI - y;
		}

		else if (y<= 0 && y>= -M_PI_2){
			yd = fabs(y) - data[1][0];
		}
		else if (y>=M_PI && y<= 1.5*M_PI){
			yd = y - M_PI;
		}
		
		else {
			yd = y;
		}if (M_PI_2 < y && y <=M_PI){
			yd = M_PI - y;
		}

		else if (y<= 0 && y>= -M_PI_2){
			yd = fabs(y);
		}
		else if (y>=M_PI && y<= 1.5*M_PI){
			yd = y - M_PI;
		}

		else if (y<=-M_PI_2 && y>=-M_PI){
			yd = fabs(-M_PI -y);
		}

		else if (y>= 1.5*M_PI && y<= 2.0*M_PI){
			yd = 2.0*M_PI - y;
		}

		else if (y<=-M_PI && y>= -1.5*M_PI){
			yd = fabs(y +M_PI);
		}

		else if (y>=2.0*M_PI && y<= 2.5*M_PI){
			yd = y - 2.0*M_PI;
		}
		
		else {
			yd = y;
		}

		if ((y>=thm && y<=M_PI_2) || (yd>=thm && yd<=M_PI_2)){
			yd = thm-pow(10,-6);
		}

		xd = sqrt(x*x - rh*rh)/(sqrt(x*x - rh*rh) + 1);

		ybar = (yd - data[1][index])/dy;
		xbar = (xd - data[0][index]) / dx;
	

		return coeff[index2][index][0] + coeff[index2][index][1] * ybar + coeff[index2][index][2] * pow(ybar, 2) + coeff[index2][index][3] * pow(ybar, 3) + coeff[index2][index][4] * xbar + coeff[index2][index][5] * xbar * ybar + coeff[index2][index][6] * xbar * pow(ybar, 2) + coeff[index2][index][7] * xbar * pow(ybar, 3) + coeff[index2][index][8] * pow(xbar, 2) + coeff[index2][index][9] * pow(xbar, 2) * ybar + coeff[index2][index][10] * pow(xbar, 2) * pow(ybar, 2) + coeff[index2][index][11] * pow(xbar, 2) * pow(ybar, 3) + coeff[index2][index][12] * pow(xbar, 3) + coeff[index2][index][13] * pow(xbar, 3) * ybar + coeff[index2][index][14] * pow(xbar, 3) * pow(ybar, 2) + coeff[index2][index][15] * pow(xbar, 3) * pow(ybar, 3);
	}
	
	double bicubic_spline_diffx(double x, double y, int index, int index2)
	{	
		double dx = data[0][index+1] - data[0][index];
		double dy = data[1][index+ xsize] - data[1][index]; 
		double yd;
		double xd;
		double ybar;
		double xbar;

		if (M_PI_2 <= y && y <=M_PI){
			yd = M_PI - y;
		}

		else if (y<= 0 && y>= -M_PI_2){
			yd = fabs(y);
		}
		else if (y>=M_PI && y<= 1.5*M_PI){
			yd = y - M_PI;
		}

		else if (y<=-M_PI_2 && y>=-M_PI){
			yd = fabs(-M_PI -y);
		}

		else if (y>= 1.5*M_PI && y<= 2.0*M_PI){
			yd = 2.0*M_PI - y;
		}

		else if (y<=-M_PI && y>= -1.5*M_PI){
			yd = fabs(y +M_PI);
		}

		else if (y>=2.0*M_PI && y<= 2.5*M_PI){
			yd = y - 2.0*M_PI;
		}
		
		else {
			yd = y;
		}

		if ((y>=thm && y<=M_PI_2) || (yd>=thm && yd<=M_PI_2)){
			yd = thm-pow(10,-6);
		}

		xd = sqrt(x*x - rh*rh)/(sqrt(x*x - rh*rh) + 1);

		ybar = (yd - data[1][index])/dy;
		xbar = (xd - data[0][index]) / dx;
		
		double _x = coeff[index2][index][4] * (1.0 / dx) + coeff[index2][index][5] * (1.0 / dx) * ybar + coeff[index2][index][6] * (1.0 / dx) * pow(ybar, 2) + coeff[index2][index][7] * (1.0 / dx) * pow(ybar, 3) + coeff[index2][index][8] * (2.0 * xbar / dx) + coeff[index2][index][9] * (2.0 * xbar / dx) * ybar + coeff[index2][index][10] * (2.0 * xbar / dx) * pow(ybar, 2) + coeff[index2][index][11] * (2.0 * xbar / dx) * pow(ybar, 3) + coeff[index2][index][12] * (3.0 * xbar * xbar / dx) + coeff[index2][index][13] * (3.0 * xbar * xbar / dx) * ybar + coeff[index2][index][14] * (3.0 * xbar * xbar / dx) * pow(ybar, 2) + coeff[index2][index][15] * (3.0 * xbar * xbar / dx) * pow(ybar, 3);
		double x_r = x/(sqrt(x*x - rh*rh)*pow( sqrt(x*x - rh*rh)+ 1,2));
		return _x*x_r;
	}
	
	double bicubic_spline_diffy(double x, double y, int index, int index2)
	{	
		double dx = data[0][index+1] - data[0][index];
		double dy = data[1][index+ xsize] - data[1][index]; 
		double yd;
		double xd;
		double ybar;
		double xbar;

		if (M_PI_2 <= y && y <=M_PI){
			yd = M_PI - y;
		}

		else if (y<= 0 && y>= -M_PI_2){
			yd = fabs(y);
		}
		else if (y>=M_PI && y<= 1.5*M_PI){
			yd = y - M_PI;
		}

		else if (y<=-M_PI_2 && y>=-M_PI){
			yd = fabs(-M_PI -y);
		}

		else if (y>= 1.5*M_PI && y<= 2.0*M_PI){
			yd = 2.0*M_PI - y;
		}

		else if (y<=-M_PI && y>= -1.5*M_PI){
			yd = fabs(y +M_PI);
		}

		else if (y>=2.0*M_PI && y<= 2.5*M_PI){
			yd = y - 2.0*M_PI;
		}
		
		else {
			yd = y;
		}

		if ((y>=thm && y<=M_PI_2) || (yd>=thm && yd<=M_PI_2)){
			yd = thm-pow(10,-6);
		}

		xd = sqrt(x*x - rh*rh)/(sqrt(x*x - rh*rh) + 1);

		ybar = (yd - data[1][index])/dy;
		xbar = (xd - data[0][index]) / dx;
		

		return coeff[index2][index][1] * (1.0 / dy) + coeff[index2][index][2] * (2.0 * ybar / dy) + coeff[index2][index][3] * (3.0 * ybar * ybar / dy) + coeff[index2][index][5] * xbar * (1.0 / dy) + coeff[index2][index][6] * xbar * (2.0 * ybar / dy) + coeff[index2][index][7] * xbar * (3.0 * ybar * ybar / dy) + coeff[index2][index][9] * pow(xbar, 2) * (1.0 / dy) + coeff[index2][index][10] * pow(xbar, 2) * (2.0 * ybar / dy) + coeff[index2][index][11] * pow(xbar, 2) * (3.0 * ybar * ybar / dy) + coeff[index2][index][13] * pow(xbar, 3) * (1.0 / dy) + coeff[index2][index][14] * pow(xbar, 3) * (2.0 * ybar / dy) + coeff[index2][index][15] * pow(xbar, 3) * (3.0 * ybar * ybar / dy);
		
	}
	
};

class Stationary_axisymmetric_photon
{
public:
	double E;		  // E is the conserved energy set later on
	double L;		  // L is the conserved angular momentum set later on
	double r_min;	  // smallest r value in data set
	double tau_final; // final proper time after which RK-4 ray tracing ends
	double th_min;
	double th_max;
	double energy_checker;
	double energy_checker_0;
	bool inspiral;
	int size;
	double M;

	// Arrays to store points for RK-4 integration. Since its adaptive time step, using a vector
	double t;
	double r;
	double th;
	double phi;
	double v_t;
	double v_r;
	double v_th;
	double v_phi;

	double x;
	double y;
	double z;

	Stationary_axisymmetric_photon(double En, double Lz, double end, double mass, class Stationary_axisymmetric_data data)
	{
		E = En;
		L = Lz;
		r_min = data.rh;
		tau_final = end;
		energy_checker = 0.0;
		energy_checker_0 = 0.0;
		inspiral = false;
		size = 0;
		th_min =0.0;
		th_max = M_PI;
		M = mass;

		t = 0.0;
		r = 0.0;
		th = 0.0;
		phi = 0.0;
		v_t = 0.0;
		v_r = 0.0;
		v_th = 0.0;
		v_phi = 0.0;
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	// Slopes for t, phi and the effective_potential equation
	double t_slope(double r, double th, double f0, double f2, double w)
	{
		return (E * g_phiphi(r, th, f2) + L * g_tphi(r, th, f2, w)) / (pow(g_tphi(r, th, f2, w), 2.0) - g_tt(r, r_min, th, f0, f2, w) * g_phiphi(r, th, f2));
	}

	double phi_slope(double r, double th, double f0, double f2, double w)
	{
		return -(E * g_tphi(r, th, f2, w) + L * g_tt(r, r_min, th, f0, f2, w)) / (pow(g_tphi(r, th, f2, w), 2.0) - g_tt(r, r_min, th, f0, f2, w) * g_phiphi(r, th, f2));
	}

	double V_eff(double r, double th, double f0, double f1, double f2, double w)
	{
		return (1.0 / g_rr(r, r_min, f1)) * ((pow(E, 2.0) * g_phiphi(r, th, f2) + pow(L, 2.0) * g_tt(r, r_min, th, f0, f2, w) + 2.0 * E * L * g_tphi(r, th, f2, w)) / (g_tt(r, r_min, th, f0, f2, w) * g_phiphi(r, th, f2) - pow(g_tphi(r, th, f2, w), 2.0)));
	}

	// Slopes for v_r and v_th

	double v_r_slope(double r, double th, double v_t, double v_r, double v_th, double v_phi, double f0, double f1, double f2, double w, double f0_r, double f1_r, double f2_r, double w_r, double f1_th)
	{
		return (-1.0 / (2.0 * g_rr(r, r_min, f1))) * (2.0 * v_r * v_th * g_rr_th(r, r_min, f1, f1_th) + pow(v_r, 2) * g_rr_r(r, r_min, f1, f1_r) - pow(v_t, 2) * g_tt_r(r, r_min, th, f0, f2, w, f0_r, f2_r, w_r) - 2.0 * v_t * v_phi * g_tphi_r(r, th, f2, w, f2_r, w_r) - pow(v_th, 2) * g_thth_r(r, f1, f1_r) - pow(v_phi, 2) * g_phiphi_r(r, th, f2, f2_r));
	}

	double v_th_slope(double r, double th, double v_t, double v_r, double v_th, double v_phi, double f0, double f1, double f2, double w, double f1_r, double f0_th, double f1_th, double f2_th, double w_th)
	{
		return (-1.0 / (2.0 * g_thth(r, f1))) * (2.0 * v_r * v_th * g_thth_r(r, f1, f1_r) - pow(v_r, 2) * g_rr_th(r, r_min, f1, f1_th) - pow(v_t, 2) * g_tt_th(r, r_min, th, f0, f2, w, f0_th, f2_th, w_th) - 2.0 * v_t * v_phi * g_tphi_th(r, th, f2, w, f2_th, w_th) + pow(v_th, 2.0) * g_thth_th(r, f1, f1_th) - pow(v_phi, 2) * g_phiphi_th(r, th, f2, f2_th));
	}

	void set_initial_data(double t0, double r0, double th0, double phi0, class Stationary_axisymmetric_data data)
	{
		t = t0;
		r = r0;
		th = th0;
		phi = phi0;

		int index = data.sorter(r0, th0);
		double f0 = data.bicubic_spline(r0, th0, index, 0);
		double f1 = data.bicubic_spline(r0, th0, index, 1);
		double f2 = data.bicubic_spline(r0, th0, index, 2);
		double w = data.bicubic_spline(r0, th0, index, 3);

		v_t = t_slope(r0, th0, f0, f2, w);
		v_phi = phi_slope(r0, th0, f0, f2, w);
		
		x = r0 * sin(th0) * cos(phi0);
		y = r0 * sin(th0) * sin(phi0);
		z = r0 * cos(th0);

		size = 1;
	}

	void RKCK(class Stationary_axisymmetric_data data, double stop)
	{

		int index = 0;
		double tau = 0;
		double dtau = 0;

		double f0 = 0.0;
		double f1 = 0.0;
		double f2 = 0.0;
		double w = 0.0;
		double f0_r = 0.0;
		double f1_r = 0.0;
		double f2_r = 0.0;
		double w_r = 0.0;
		double f0_th = 0.0;
		double f1_th = 0.0;
		double f2_th = 0.0;
		double w_th = 0.0;

		double r_d = 0.0;
		double th_d = 0.0;
		double v_r_d = 0.0;
		double v_th_d = 0.0;

		double v_r_k1;
		double r_k1;
		double v_th_k1;
		double th_k1;
		double t_k1;
		double phi_k1;

		double v_r_k2;
		double r_k2;
		double v_th_k2;
		double th_k2;
		double t_k2;
		double phi_k2;

		double v_r_k3;
		double r_k3;
		double v_th_k3;
		double th_k3;
		double t_k3;
		double phi_k3;

		double v_r_k4;
		double r_k4;
		double v_th_k4;
		double th_k4;
		double t_k4;
		double phi_k4;

		double v_r_k5;
		double r_k5;
		double v_th_k5;
		double th_k5;
		double t_k5;
		double phi_k5;

		double v_r_k6;
		double r_k6;
		double v_th_k6;
		double th_k6;
		double t_k6;
		double phi_k6;

		while (tau > tau_final)
		{
			
			if (r <= 3.0 * M || (th>=-0.15 && th<=0.15 ) || (th>=M_PI-0.15 && th<=M_PI+0.15))
            {
                dtau = -pow(10, -3);
            }

            else if (r > 3.0 * M && r <= 10.0 * M)
            {
                dtau = -pow(10, -2);
            }

            else if (r > 10.0 * M && r <= 50.0 * M)
            {
                dtau = -5.0*pow(10, -2);
            }

            else
            {
                dtau = -pow(10, -1);
            }

			index = data.sorter(r, th);
			f0 = data.bicubic_spline(r, th, index, 0);
			f1 = data.bicubic_spline(r, th, index, 1);
			f2 = data.bicubic_spline(r, th, index, 2);
			w = data.bicubic_spline(r, th, index, 3);
			f0_r = data.bicubic_spline_diffx(r, th, index, 0);
			f1_r = data.bicubic_spline_diffx(r, th, index, 1);
			f2_r = data.bicubic_spline_diffx(r, th, index, 2);
			w_r = data.bicubic_spline_diffx(r, th, index, 3);
			f0_th = data.bicubic_spline_diffy(r, th, index, 0);
			f1_th = data.bicubic_spline_diffy(r, th, index, 1);
			f2_th = data.bicubic_spline_diffy(r, th, index, 2);
			w_th = data.bicubic_spline_diffy(r, th, index, 3);

			v_r_k1 = (dtau)*v_r_slope(r, th, v_t, v_r, v_th, v_phi, f0, f1, f2, w, f0_r, f1_r, f2_r, w_r, f1_th);
			r_k1 = (dtau)*v_r;
			v_th_k1 = (dtau)*v_th_slope(r, th, v_t, v_r, v_th, v_phi, f0, f1, f2, w, f1_r, f0_th, f1_th, f2_th, w_th);
			th_k1 = (dtau)*v_th;
			t_k1 = (dtau)*v_t;
			phi_k1 = (dtau)*v_phi;

			r_d = r + r_k1 / 5.0;
			th_d = th + th_k1 / 5.0;

			if (r_d <= r_min)
			{
				inspiral = true;
				break;
			}

			if (r_d >= stop)
			{
				break;
			}

			
			index = data.sorter(r_d, th_d);
			f0 = data.bicubic_spline(r_d, th_d, index, 0);
			f1 = data.bicubic_spline(r_d, th_d, index, 1);
			f2 = data.bicubic_spline(r_d, th_d, index, 2);
			w = data.bicubic_spline(r_d, th_d, index, 3);
			f0_r = data.bicubic_spline_diffx(r_d, th_d, index, 0);
			f1_r = data.bicubic_spline_diffx(r_d, th_d, index, 1);
			f2_r = data.bicubic_spline_diffx(r_d, th_d, index, 2);
			w_r = data.bicubic_spline_diffx(r_d, th_d, index, 3);
			f0_th = data.bicubic_spline_diffy(r_d, th_d, index, 0);
			f1_th = data.bicubic_spline_diffy(r_d, th_d, index, 1);
			f2_th = data.bicubic_spline_diffy(r_d, th_d, index, 2);
			w_th = data.bicubic_spline_diffy(r_d, th_d, index, 3);
			v_r_d = v_r + v_r_k1 / 5.0;
			v_th_d = v_th + v_th_k1 / 5.0;
			v_t = t_slope(r_d, th_d, f0, f2, w);
			v_phi = phi_slope(r_d, th_d, f0, f2, w);

			v_r_k2 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f0_r, f1_r, f2_r, w_r, f1_th);
			r_k2 = (dtau) * (v_r_d);
			v_th_k2 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f1_r, f0_th, f1_th, f2_th, w_th);
			th_k2 = (dtau) * (v_th_d);
			t_k2 = (dtau)*v_t;
			phi_k2 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 40.0) + r_k2 * (9.0 / 40.0);
			th_d = th + th_k1 * (3.0 / 40.0) + th_k2 * (9.0 / 40.0);

			if (r_d <= r_min)
			{
				inspiral = true;
				break;
			}

			if (r_d >= stop)
			{
				break;
			}
			
			index = data.sorter(r_d, th_d);
			f0 = data.bicubic_spline(r_d, th_d, index, 0);
			f1 = data.bicubic_spline(r_d, th_d, index, 1);
			f2 = data.bicubic_spline(r_d, th_d, index, 2);
			w = data.bicubic_spline(r_d, th_d, index, 3);
			f0_r = data.bicubic_spline_diffx(r_d, th_d, index, 0);
			f1_r = data.bicubic_spline_diffx(r_d, th_d, index, 1);
			f2_r = data.bicubic_spline_diffx(r_d, th_d, index, 2);
			w_r = data.bicubic_spline_diffx(r_d, th_d, index, 3);
			f0_th = data.bicubic_spline_diffy(r_d, th_d, index, 0);
			f1_th = data.bicubic_spline_diffy(r_d, th_d, index, 1);
			f2_th = data.bicubic_spline_diffy(r_d, th_d, index, 2);
			w_th = data.bicubic_spline_diffy(r_d, th_d, index, 3);
			v_r_d = v_r + v_r_k1 * (3.0 / 40.0) + v_r_k2 * (9.0 / 40.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 40.0) + v_th_k2 * (9.0 / 40.0);
			v_t = t_slope(r_d, th_d, f0, f2, w);
			v_phi = phi_slope(r_d, th_d, f0, f2, w);

			v_r_k3 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f0_r, f1_r, f2_r, w_r, f1_th);
			r_k3 = (dtau) * (v_r_d);
			v_th_k3 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f1_r, f0_th, f1_th, f2_th, w_th);
			th_k3 = (dtau) * (v_th_d);
			t_k3 = (dtau)*v_t;
			phi_k3 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 10.0) + r_k2 * (-9.0 / 10.0) + r_k3 * (6.0 / 5.0);
			th_d = th + th_k1 * (3.0 / 10.0) + th_k2 * (-9.0 / 10.0) + th_k3 * (6.0 / 5.0);

			if (r_d <= r_min)
			{
				inspiral = true;
				break;
			}

			
			if (r_d >= stop)
			{
				break;
			}
			
			index = data.sorter(r_d, th_d);
			f0 = data.bicubic_spline(r_d, th_d, index, 0);
			f1 = data.bicubic_spline(r_d, th_d, index, 1);
			f2 = data.bicubic_spline(r_d, th_d, index, 2);
			w = data.bicubic_spline(r_d, th_d, index, 3);
			f0_r = data.bicubic_spline_diffx(r_d, th_d, index, 0);
			f1_r = data.bicubic_spline_diffx(r_d, th_d, index, 1);
			f2_r = data.bicubic_spline_diffx(r_d, th_d, index, 2);
			w_r = data.bicubic_spline_diffx(r_d, th_d, index, 3);
			f0_th = data.bicubic_spline_diffy(r_d, th_d, index, 0);
			f1_th = data.bicubic_spline_diffy(r_d, th_d, index, 1);
			f2_th = data.bicubic_spline_diffy(r_d, th_d, index, 2);
			w_th = data.bicubic_spline_diffy(r_d, th_d, index, 3);
			v_r_d = v_r + v_r_k1 * (3.0 / 10.0) + v_r_k2 * (-9.0 / 10.0) + v_r_k3 * (6.0 / 5.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 10.0) + v_th_k2 * (-9.0 / 10.0) + v_th_k3 * (6.0 / 5.0);
			v_t = t_slope(r_d, th_d, f0, f2, w);
			v_phi = phi_slope(r_d, th_d, f0, f2, w);

			v_r_k4 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f0_r, f1_r, f2_r, w_r, f1_th);
			r_k4 = (dtau) * (v_r_d);
			v_th_k4 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f1_r, f0_th, f1_th, f2_th, w_th);
			th_k4 = (dtau) * (v_th_d);
			t_k4 = (dtau)*v_t;
			phi_k4 = (dtau)*v_phi;

			r_d = r + r_k1 * (-11.0 / 54.0) + r_k2 * (5.0 / 2.0) + r_k3 * (-70.0 / 27.0) + r_k4 * (35.0 / 27.0);
			th_d = th + th_k1 * (-11.0 / 54.0) + th_k2 * (5.0 / 2.0) + th_k3 * (-70.0 / 27.0) + th_k4 * (35.0 / 27.0);

			if (r_d <= r_min)
			{
				inspiral = true;
				break;
			}

			if (r_d >= stop)
			{
				break;
			}

			
			index = data.sorter(r_d, th_d);
			f0 = data.bicubic_spline(r_d, th_d, index, 0);
			f1 = data.bicubic_spline(r_d, th_d, index, 1);
			f2 = data.bicubic_spline(r_d, th_d, index, 2);
			w = data.bicubic_spline(r_d, th_d, index, 3);
			f0_r = data.bicubic_spline_diffx(r_d, th_d, index, 0);
			f1_r = data.bicubic_spline_diffx(r_d, th_d, index, 1);
			f2_r = data.bicubic_spline_diffx(r_d, th_d, index, 2);
			w_r = data.bicubic_spline_diffx(r_d, th_d, index, 3);
			f0_th = data.bicubic_spline_diffy(r_d, th_d, index, 0);
			f1_th = data.bicubic_spline_diffy(r_d, th_d, index, 1);
			f2_th = data.bicubic_spline_diffy(r_d, th_d, index, 2);
			w_th = data.bicubic_spline_diffy(r_d, th_d, index, 3);
			v_r_d = v_r + v_r_k1 * (-11.0 / 54.0) + v_r_k2 * (5.0 / 2.0) + v_r_k3 * (-70.0 / 27.0) + v_r_k4 * (35.0 / 27.0);
			v_th_d = v_th + v_th_k1 * (-11.0 / 54.0) + v_th_k2 * (5.0 / 2.0) + v_th_k3 * (-70.0 / 27.0) + v_th_k4 * (35.0 / 27.0);
			v_t = t_slope(r_d, th_d, f0, f2, w);
			v_phi = phi_slope(r_d, th_d, f0, f2, w);

			v_r_k5 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f0_r, f1_r, f2_r, w_r, f1_th);
			r_k5 = (dtau) * (v_r_d);
			v_th_k5 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f1_r, f0_th, f1_th, f2_th, w_th);
			th_k5 = (dtau) * (v_th_d);
			t_k5 = (dtau)*v_t;
			phi_k5 = (dtau)*v_phi;

			r_d = r + r_k1 * (1631.0 / 55296.0) + r_k2 * (175.0 / 512.0) + r_k3 * (575.0 / 13824.0) + r_k4 * (44275.0 / 110592.0) + r_k5 * (253.0 / 4096.0);
			th_d = th + th_k1 * (1631.0 / 55296.0) + th_k2 * (175.0 / 512.0) + th_k3 * (575.0 / 13824.0) + th_k4 * (44275.0 / 110592.0) + th_k5 * (253.0 / 4096.0);
			
			if (r_d <= r_min)
			{
				inspiral = true;
				break;
			}

			
			if (r_d >= stop)
			{
				break;
			}
			
			index = data.sorter(r_d, th_d);
			
			f0 = data.bicubic_spline(r_d, th_d, index, 0);
			
			f1 = data.bicubic_spline(r_d, th_d, index, 1);
			f2 = data.bicubic_spline(r_d, th_d, index, 2);
			w = data.bicubic_spline(r_d, th_d, index, 3);
			f0_r = data.bicubic_spline_diffx(r_d, th_d, index, 0);
			f1_r = data.bicubic_spline_diffx(r_d, th_d, index, 1);
			f2_r = data.bicubic_spline_diffx(r_d, th_d, index, 2);
			w_r = data.bicubic_spline_diffx(r_d, th_d, index, 3);
			f0_th = data.bicubic_spline_diffy(r_d, th_d, index, 0);
			f1_th = data.bicubic_spline_diffy(r_d, th_d, index, 1);
			f2_th = data.bicubic_spline_diffy(r_d, th_d, index, 2);
			w_th = data.bicubic_spline_diffy(r_d, th_d, index, 3);
			v_r_d = v_r + v_r_k1 * (1631.0 / 55296.0) + v_r_k2 * (175.0 / 512.0) + v_r_k3 * (575.0 / 13824.0) + v_r_k4 * (44275.0 / 110592.0) + v_r_k5 * (253.0 / 4096.0);
			v_th_d = v_th + v_th_k1 * (1631.0 / 55296.0) + v_th_k2 * (175.0 / 512.0) + v_th_k3 * (575.0 / 13824.0) + v_th_k4 * (44275.0 / 110592.0) + v_th_k5 * (253.0 / 4096.0);
			v_t = t_slope(r_d, th_d, f0, f2, w);
			v_phi = phi_slope(r_d, th_d, f0, f2, w);

			v_r_k6 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f0_r, f1_r, f2_r, w_r, f1_th);
			r_k6 = (dtau) * (v_r_d);
			v_th_k6 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, f0, f1, f2, w, f1_r, f0_th, f1_th, f2_th, w_th);
			th_k6 = (dtau) * (v_th_d);
			t_k6 = (dtau)*v_t;
			phi_k6 = (dtau)*v_phi;

			r_d = r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6;
			th_d = th + (37.0 / 378.0) * th_k1 + (250.0 / 621.0) * th_k3 + (125.0 / 594.0) * th_k4 + (512.0 / 1771.0) * th_k6;
			
			if (r_d <= r_min)
			{
				inspiral = true;
				break;
			}

			if (r_d >= stop)
			{
				break;
			}

			t = t + (37.0 / 378.0) * t_k1 + (250.0 / 621.0) * t_k3 + (125.0 / 594.0) * t_k4 + (512.0 / 1771.0) * t_k6;
			r = r_d;
			th = th_d;
			phi = phi + (37.0 / 378.0) * phi_k1 + (250.0 / 621.0) * phi_k3 + (125.0 / 594.0) * phi_k4 + (512.0 / 1771.0) * phi_k6;
			
			index = data.sorter(r, th);
			
			f0 = data.bicubic_spline(r, th, index, 0);
			f1 = data.bicubic_spline(r, th, index, 1);
			f2 = data.bicubic_spline(r, th, index, 2);
			w = data.bicubic_spline(r, th, index, 3);

			v_t = t_slope(r, th, f0, f2, w);
			v_phi = phi_slope(r, th, f0, f2, w);
			v_r = v_r + (37.0 / 378.0) * v_r_k1 + (250.0 / 621.0) * v_r_k3 + (125.0 / 594.0) * v_r_k4 + (512.0 / 1771.0) * v_r_k6;
			v_th = v_th + (37.0 / 378.0) * v_th_k1 + (250.0 / 621.0) * v_th_k3 + (125.0 / 594.0) * v_th_k4 + (512.0 / 1771.0) * v_th_k6;

			// converting to cartesian coordinates
			x = r * sin(th) * cos(phi);
			y = r * sin(th) * sin(phi);
			z = r * cos(th);

			energy_checker = energy_checker + fabs(pow(v_r, 2.0) + pow(v_th, 2.0) * g_thth(r, f1) / g_rr(r, r_min, f1) + V_eff(r, th, f0, f1, f2, w));
			tau = tau + dtau;
			size = size + 1;
			
		}
	}
	
	

	int colour_coder()
	{
		int colour;

		if (inspiral)
		{
			colour = 0;
		}

		else
		{
			if (y < 0 && z > 0)
			{
				colour = 1;
			}
			else if (y > 0 && z > 0)
			{
				colour = 2;
			}
			else if (y < 0 && z < 0)
			{
				colour = 3;
			}
			else if (y > 0 && z < 0)
			{
				colour = 4;
			}
		}
		return colour;
	}
};

int main()
{
	double rh = 0.0359508947168857;
    double M = 0.504708503201793;
    std::ofstream myFile("gc=0, q=0.txt");

	// make data class that will create/read data and the splines.
	Stationary_axisymmetric_data data(120, 30 ,rh);

	// make the grid of data
	data.make_data();

	std::cout << "Read the data file"
			  << "\n";

	// make coefficients for the spline
	data.bicubic_spline_coefficient();

	std::cout << "Made the splines"
			  << "\n";
	

	// make a photon of given E, L and integration limit+parameters
	Stationary_axisymmetric_photon photon(0.0, 0.0, -1000, M,  data);


	// Shadow formation

	// this is the new parametrization
	double *alpha = new double[125];
	double *beta = new double[250];

	double lim = 0.73;
	linspace(0, lim, 125, alpha);
	linspace(-lim, lim, 250, beta);


	double x;
	double y;
	double phi = 0.0;
	double th = 1.57;
	double E;
	double L;
	double v_r_0;
	double v_th_0;
	double r = 15.0*M;
	double mod_P = 1.0;
	double A_t;
	double P_r_obs;
	double P_th_obs;
	double P_phi_obs;

	double f0;
	double f1;
	double f2;
	double w;
	int index;

	index = data.sorter(r, th);
	f0 = data.bicubic_spline(r, th, index, 0);
	f1 = data.bicubic_spline(r, th, index, 1);
	f2 = data.bicubic_spline(r, th, index, 2);
	w = data.bicubic_spline(r, th, index, 3);

	std::cout<<"The value of gphiphi is: "<<g_phiphi(r,th,f2)<<"\n";

	for (int i = 0; i < 250; i++)
	{
		for (int j = 0; j < 125; j++)
		{
			x = -sin(beta[i]) * r;
			y = sin(alpha[j]) * r;
			P_r_obs = mod_P * cos(alpha[j]) * cos(beta[i]);
			P_th_obs = mod_P * sin(alpha[j]);
			P_phi_obs = mod_P * cos(alpha[j]) * sin(beta[i]);

			A_t = sqrt(g_phiphi(r, th, f2) / (pow(g_tphi(r, th, f2, w), 2) - g_tt(r, photon.r_min, th, f0, f2, w) * g_phiphi(r, th, f2)));

			E = (mod_P / A_t) - P_phi_obs * g_tphi(r, th, f2, w) / sqrt(g_phiphi(r, th, f2));
			L = P_phi_obs * sqrt(g_phiphi(r, th, f2));
			v_r_0 = P_r_obs / sqrt(g_rr(r, photon.r_min, f1));
			v_th_0 = P_th_obs / sqrt(g_thth(r, f1));

			photon.E = E;
			photon.L = L;
			photon.energy_checker = 0.0;
			photon.size = 0.0;
			photon.inspiral = false;

			// set initial data, check if its in allowed region
			photon.set_initial_data(0.0, r, th, phi, data);
			photon.v_r = v_r_0;
			photon.v_th = v_th_0;
			photon.energy_checker_0 = fabs(pow(photon.v_r, 2.0) + pow(photon.v_th, 2.0) * g_thth(r, f1) / g_rr(r, photon.r_min, f1) + photon.V_eff(r, th, f0, f1, f2, w));
			photon.energy_checker = photon.energy_checker + fabs(pow(photon.v_r, 2.0) + pow(photon.v_th, 2.0) * g_thth(r, f1) / g_rr(r, photon.r_min, f1) + photon.V_eff(r, th, f0, f1, f2, w));
			// do the RKCK and error analysis
			photon.RKCK(data, 15.0*M);
			myFile << x << "," << y << "," << photon.colour_coder() << "," << photon.r << "," << photon.energy_checker_0 << "," << photon.energy_checker / photon.size << ","<< photon.th<< ","<< photon.size<< ","<<photon.E<<","<< photon.L<<"\n";
			//std::cout<< j <<"\n";
		}
		std::cout << i << "\n";
	}
	myFile.close();

}