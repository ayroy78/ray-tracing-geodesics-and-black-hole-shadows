#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
//#include <vector>
#include <fstream>
#include <chrono>

//These are metric functions for Kerr
double sigma(double r, double th, double a) {
	return pow(r, 2) + pow(a * cos(th), 2);
}

double delta(double M, double r, double a) {
	return pow(r, 2) + pow(a, 2) - 2 * M * r;
}

double g_tt(double M, double r, double th, double a) {
	return -1 * (1 - (2 * M * r) / (sigma(r, th, a)));
}

double g_tphi(double M, double r, double th, double a) {
	return -(2 * M * r * a * pow(sin(th), 2.0)) / (sigma(r, th, a));
}

double g_rr(double M, double r, double th, double a) {
	return sigma(r, th, a)/delta(M, r, a);
}

double g_thth(double r, double th, double a) {
	return sigma(r, th, a);
}

double g_phiphi(double M, double r, double th, double a) {
	return (pow(a, 2) + pow(r, 2) + (2 * M * r * pow(a * sin(th), 2.0)) / (sigma(r, th, a))) * pow(sin(th), 2);
}

//derivatives for radial equation
double g_rr_th(double M, double r, double th, double a) {
	return -(sin(2 * th) * pow(a, 2))/delta(M,r,a);
}

double g_rr_r(double M, double r, double th, double a) {
	return 2.0 * (r * delta(M, r, a) - (r - M) * sigma(r, th, a))/pow(delta(M,r,a),2);

}

double g_tt_r(double M, double r, double th, double a) {
	return ((2 * M) / (pow(sigma(r, th, a), 2))) * (sigma(r, th, a) - 2 * pow(r, 2));
}

double g_tphi_r(double M, double r, double th, double a) {
	return (2 * M * a * pow(sin(th), 2) / pow(sigma(r, th, a), 2)) * (2 * pow(r, 2) - sigma(r, th, a));
}

double g_phiphi_r(double M, double r, double th, double a) {
	return (pow(sin(th), 2)) * (2 * r + ((2 * M * pow(a * sin(th), 2)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) - 2 * pow(r, 2)));
}

double g_thth_r(double r) {
	return 2 * r;
}

//now for the theta equation of motion

double g_tt_th(double M, double r, double th, double a) {
	return (2 * M * r * sin(2 * th) * pow(a, 2)) / pow(sigma(r, th, a), 2);
}

double g_tphi_th(double M, double r, double th, double a) {
	return -((2 * M * r * a * sin(2 * th)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) + pow(a * sin(th), 2));
}

double g_thth_th(double th, double a) {
	return -sin(2 * th) * pow(a, 2);
}

double g_phiphi_th(double M, double r, double th, double a) {
	return sin(2 * th)*(pow(r,2) + pow(a,2) +2.0*M*r*pow(a*sin(th),2)/sigma(r,th,a)) + ((2 * M * r * sin(2 * th) * pow(a * sin(th), 2)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) + pow(a * sin(th), 2));
}

double V_eff(double E, double L, double M, double r, double th, double a) {
	return (pow(g_rr(M, r, th, a), -1)) * ((g_phiphi(M, r, th, a) * pow(E, 2) + g_tt(M, r, th, a) * pow(L, 2) + 2 * E * L * g_tphi(M, r, th, a)) / (-pow(g_tphi(M, r, th, a), 2) + g_tt(M, r, th, a) * g_phiphi(M, r, th, a)));
}

//this is a function to create linearly separated points. 
void linspace(double start, double end, double numpoints, double* points) {
	//double * points  = new double [numpoints];
	double dx = (end - start) / (numpoints - 1);
	for (int i = 0; i < numpoints; i++) {
		points[i] = (start + i * dx);
	}
}


class Stationary_axisymmetric_data {
public:
	double** data = new double* [17];
	double*** coeff = new double** [15];
	double M;
	double a;
	int xsize;
	int ysize;
	double dx;
	double dy;

	Stationary_axisymmetric_data(int x_size, int y_size, double x_sep, double y_sep, double mass, double angmom)
	{	
		M = mass;
		a = angmom;
		xsize = x_size;
		ysize = y_size;
		dx = x_sep;
		dy = y_sep;

		for (int i = 0; i < 17; i++) {
			data[i] = new double[xsize * ysize];
		}

		for (int i = 0; i < 15; i++) {
			coeff[i] = new double* [xsize * ysize];
		}

		for (int i = 0; i < 15; i++) {
			for (int j = 0; j < xsize * ysize; j++) {
				coeff[i][j] = new double[16];
			}
		}


	}

	void make_data(double* x, double* y, double M, double a) {
		int counter = 0;
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {

				data[0][counter] = x[i];
				data[1][counter] = y[j];
				data[2][counter] = g_tt(M, x[i], y[j], a);
				data[3][counter] = g_tphi(M, x[i], y[j], a);
				data[4][counter] = g_rr(M, x[i], y[j], a)*delta(M,x[i],a) - pow(x[i],2);
				data[5][counter] = g_thth(x[i], y[j], a) - pow(x[i],2);
				data[6][counter] = g_phiphi(M, x[i], y[j], a) - pow(x[i]*sin(y[j]),2);
				data[7][counter] = g_rr_th(M, x[i], y[j], a)*delta(M, x[i],a);
				data[8][counter] = g_rr_r(M, x[i], y[j], a)*pow(delta(M, x[i],a),2);
				data[9][counter] = (g_tt_r(M, x[i], y[j], a));
				data[10][counter] = (g_tphi_r(M, x[i], y[j], a));
				data[11][counter] = (g_phiphi_r(M, x[i], y[j], a));
				data[12][counter] = (g_thth_r(x[i]));
				data[13][counter] = (g_tt_th(M, x[i], y[j], a));
				data[14][counter] = (g_tphi_th(M, x[i], y[j], a));
				data[15][counter] = (g_thth_th(y[j], a));
				data[16][counter] = (g_phiphi_th(M, x[i], y[j], a)) - pow(x[i],2)*sin(2*y[j]);

				counter++;
			}
		}

	}

	int sorter(double x, double y) {

		int i = int((x - data[0][0]) / dx);
		int j = int((y - data[1][0]) / dy);
		return i * ysize + j;

	}


	double diff_x(int index, int index2) { //index is the point in data array, index2 is the metric term

		if (0 <= index && index <= ysize - 1) {
			return 0.5 * (-3.0 * data[index2][index] + 4.0 * data[index2][index + ysize] - data[index2][index + 2 * ysize]) / (data[0][index + ysize] - data[0][index]);

		}
		else if (xsize * ysize - ysize <= index && index <= xsize * ysize - 1) {
			return  0.5 * (+3.0 * data[index2][index] - 4.0 * data[index2][index - ysize] + data[index2][index - 2 * ysize]) / (data[0][index] - data[0][index - ysize]);
		}
		else {
			return (data[index2][index + ysize] - data[index2][index - ysize]) / (data[0][index + ysize] - data[0][index - ysize]);
		}

	}

	double diff_y(int index, int index2) { //index is the point in data array, index2 is the metric term

		if (index % ysize == 0) {
			return 0.5 * (-3.0 * data[index2][index] + 4.0 * data[index2][index + 1] - data[index2][index + 2]) / (data[1][index + 1] - data[1][index]);

		}
		else if ((index + 1) % ysize == 0) {
			return 0.5 * (3.0 * data[index2][index] - 4.0 * data[index2][index - 1] + data[index2][index - 2]) / (data[1][index] - data[1][index - 1]);

		}
		else {
			return (data[index2][index + 1] - data[index2][index - 1]) / (data[1][index + 1] - data[1][index - 1]);

		}

	}

	double diff_xy(int index, int index2) {

		if (index == 0) {
			return (1.0 / 2.0) * (0.5 * (-3.0 * diff_x(index, index2) + 4.0 * diff_x(index + 1, index2) - diff_x(index + 2, index2)) / (data[1][index + 1] - data[1][index]) + 0.5 * (-3.0 * diff_y(index, index2) + 4.0 * diff_y(index + ysize, index2) - diff_y(index + 2 * ysize, index2)) / (data[0][index + ysize] - data[0][index]));

		}
		else if (index == ysize - 1) {
			return (1.0 / 2.0) * (0.5 * (+3.0 * diff_x(index, index2) - 4.0 * diff_x(index - 1, index2) + diff_x(index - 2, index2)) / (data[1][index] - data[1][index - 1]) + 0.5 * (-3.0 * diff_y(index, index2) + 4.0 * diff_y(index + ysize, index2) - diff_y(index + 2 * ysize, index2)) / (data[0][index + ysize] - data[0][index]));

		}
		else if (index == xsize * ysize - ysize) {
			return  (1.0 / 2.0) * (0.5 * (-3.0 * diff_x(index, index2) + 4.0 * diff_x(index + 1, index2) - diff_x(index + 2, index2)) / (data[1][index + 1] - data[1][index]) + 0.5 * (3.0 * diff_y(index, index2) - 4.0 * diff_y(index - ysize, index2) + diff_y(index - 2 * ysize, index2)) / (data[0][index] - data[0][index - ysize]));

		}
		else if (index == xsize * ysize - 1) {
			return (1.0 / 2.0) * (0.5 * (3.0 * diff_x(index, index2) - 4.0 * diff_x(index - 1, index2) + diff_x(index - 2, index2)) / (data[1][index] - data[1][index - 1]) + 0.5 * (3.0 * diff_y(index, index2) - 4.0 * diff_y(index - ysize, index2) + diff_y(index - 2 * ysize, index2)) / (data[0][index] - data[0][index - ysize]));

		}

		else if (0 < index && index < ysize - 1) {
			return (1.0 / 2.0) * ((diff_x(index + 1, index2) - diff_x(index - 1, index2)) / (data[1][index + 1] - data[1][index - 1]) + 0.5 * (-3.0 * diff_y(index, index2) + 4.0 * diff_y(index + ysize, index2) - diff_y(index + 2 * ysize, index2)) / (data[0][index + ysize] - data[0][index]));

		}
		else if (xsize * ysize - ysize < index && index < xsize * ysize - 1) {

			return (1.0 / 2.0) * ((diff_x(index + 1, index2) - diff_x(index - 1, index2)) / (data[1][index + 1] - data[1][index - 1]) + 0.5 * (+3.0 * diff_y(index, index2) - 4.0 * diff_y(index - ysize, index2) + diff_y(index - 2 * ysize, index2)) / (data[0][index] - data[0][index - ysize]));

		}

		else if (index % ysize == 0) {

			return (1.0 / 2.0) * (0.5 * (-3.0 * diff_x(index, index2) + 4.0 * diff_x(index + 1, index2) - diff_x(index + 2, index2)) / (data[1][index + 1] - data[1][index]) + (diff_y(index + ysize, index2) - diff_y(index - ysize, index2)) / (data[0][index + ysize] - data[0][index - ysize]));

		}
		else if ((index + 1) % ysize == 0) {
			return (1.0 / 2.0) * (0.5 * (+3.0 * diff_x(index, index2) - 4.0 * diff_x(index - 1, index2) + diff_x(index - 2, index2)) / (data[1][index] - data[1][index - 1]) + (diff_y(index + ysize, index2) - diff_y(index - ysize, index2)) / (data[0][index + ysize] - data[0][index - ysize]));

		}

		else {
			return (1.0 / 2.0) * ((diff_x(index + 1, index2) - diff_x(index - 1, index2)) / (data[1][index + 1] - data[1][index - 1]) + (diff_y(index + ysize, index2) - diff_y(index - ysize, index2)) / (data[0][index + ysize] - data[0][index - ysize]));

		}

	}

	void bicubic_spline_coefficient() {
		int index;
		int lower_right_index;
		int upper_left_index;
		int upper_right_index;
		for (int i = 0; i < xsize * ysize; i++) {
			if ((i + 1) % ysize == 0 || (xsize * ysize - ysize <= i && i <= xsize * ysize - 1)) {
				continue;
			}
			index = i;
			lower_right_index = index + ysize;
			upper_left_index = index + 1;
			upper_right_index = index + ysize + 1;

			for (int j = 2; j < 17; j++) {

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

	double bicubic_spline(double x, double y, int index, int index2) {
		double ybar = (y - data[1][index]) / (dy);
		double xbar = (x - data[0][index]) / (dx);
		double splineval = coeff[index2][index][0] + coeff[index2][index][1] * ybar + coeff[index2][index][2] * pow(ybar, 2) + coeff[index2][index][3] * pow(ybar, 3) + coeff[index2][index][4] * xbar + coeff[index2][index][5] * xbar * ybar + coeff[index2][index][6] * xbar * pow(ybar, 2) + coeff[index2][index][7] * xbar * pow(ybar, 3) + coeff[index2][index][8] * pow(xbar, 2) + coeff[index2][index][9] * pow(xbar, 2) * ybar + coeff[index2][index][10] * pow(xbar, 2) * pow(ybar, 2) + coeff[index2][index][11] * pow(xbar, 2) * pow(ybar, 3) + coeff[index2][index][12] * pow(xbar, 3) + coeff[index2][index][13] * pow(xbar, 3) * ybar + coeff[index2][index][14] * pow(xbar, 3) * pow(ybar, 2) + coeff[index2][index][15] * pow(xbar, 3) * pow(ybar, 3);
		if (index2 == 2) {
			return (splineval + pow(x,2))/delta(M, x, a);
		}
		else if (index2 == 3){
			return splineval + pow(x,2);
		}
		else if (index2 ==4){
			return splineval + pow(x*sin(y),2);
		}
		else if (index2 == 5) {
			return splineval/delta(M, x, a);
		}
		else if (index2 == 6) {
			return splineval/pow(delta(M, x, a),2);
		}
		else if (index2 == 14) {
			return splineval + pow(x,2)*sin(2.0*y);
		}
		else {
			return splineval;
		}
		


	}


};



class Stationary_axisymmetric_photon {
public:
	double M;
	double a;
	double E; //E is the conserved energy set later on
	double L; //L is the conserved angular momentum set later on
	double Q; //carter constant set later on
	double r_min; //smallest r value in data set
	double r_max; //largest r value in data set, // imp: since we don't know mass and angmom, can't find event hor...
	double tau_final; // final proper time after which RK-4 ray tracing ends
	double energy_checker;
	double carter_constant;
	bool inspiral;
	int size;
	

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

	Stationary_axisymmetric_photon(double En, double Lz, double end, class Stationary_axisymmetric_data data, double mass, double angmom)
	{	
		M = mass;
		a = angmom;
		E = En;
		L = Lz;
		Q = 0.0;
		r_min = data.data[0][0];
		r_max = data.data[0][data.xsize * data.ysize - 1];
		tau_final = end;
		energy_checker = 0.0;
		carter_constant =0.0;
		inspiral=false;
		size = 0.0;
		
		t =0.0;
		r =0.0;
		th =0.0;
		phi =0.0;
		v_t = 0.0;
		v_r = 0.0;
		v_th = 0.0;
		v_phi = 0.0;

		x=0.0;
		y=0.0;
		z=0.0;

	}

	//Slopes for t, phi and the effective_potential equation

	double t_slope(double r, double th, int index, class Stationary_axisymmetric_data data) {
		return (E * data.bicubic_spline(r, th, index, 4) + L * data.bicubic_spline(r, th, index, 1)) / (pow(data.bicubic_spline(r, th, index, 1), 2) - data.bicubic_spline(r, th, index, 0) * data.bicubic_spline(r, th, index, 4));
	}

	double phi_slope(double r, double th, int index, class Stationary_axisymmetric_data data) {
		return -(E * data.bicubic_spline(r, th, index, 1) + L * data.bicubic_spline(r, th, index, 0)) / (pow(data.bicubic_spline(r, th, index, 1), 2) - data.bicubic_spline(r, th, index, 0) * data.bicubic_spline(r, th, index, 4));
	}

	double V_eff(double r, double th, int index, class Stationary_axisymmetric_data data) {
		return (pow(data.bicubic_spline(r, th, index, 2), -1)) * ((data.bicubic_spline(r, th, index, 4) * pow(E, 2) + data.bicubic_spline(r, th, index, 0) * pow(L, 2) + 2 * E * L * data.bicubic_spline(r, th, index, 1)) / (-pow(data.bicubic_spline(r, th, index, 1), 2) + data.bicubic_spline(r, th, index, 0) * data.bicubic_spline(r, th, index, 4)));
	}

	//Slopes for v_r and v_th

	double v_r_slope(double r, double th, double v_t, double v_r, double v_th, double v_phi, int index, class Stationary_axisymmetric_data data) {
		return  (-1 / (2 * data.bicubic_spline(r, th, index, 2))) * (2 * v_r * v_th * data.bicubic_spline(r, th, index, 5) + pow(v_r, 2) * data.bicubic_spline(r, th, index, 6) - pow(v_t, 2) * data.bicubic_spline(r, th, index, 7) - 2 * v_t * v_phi * data.bicubic_spline(r, th, index, 8) - pow(v_th, 2) * data.bicubic_spline(r, th, index, 10) - pow(v_phi, 2) * data.bicubic_spline(r, th, index, 9));
	}

	double v_th_slope(double r, double th, double v_t, double v_r, double v_th, double v_phi, int index, class Stationary_axisymmetric_data data) {
		return (-1 / (2 * data.bicubic_spline(r, th, index, 3))) * (2 * v_r * v_th * data.bicubic_spline(r, th, index, 10) - pow(v_r, 2) * data.bicubic_spline(r, th, index, 5) - pow(v_t, 2) * data.bicubic_spline(r, th, index, 11) - 2 * v_t * v_phi * data.bicubic_spline(r, th, index, 12) + pow(v_th, 2) * data.bicubic_spline(r, th, index, 13) - pow(v_phi, 2) * data.bicubic_spline(r, th, index, 14));
	}

	void set_initial_data(double t0, double r0, double th0, double phi0, class Stationary_axisymmetric_data data) {
		t = t0;
		r = r0;
		th = th0;
		phi = phi0;

		int index = data.sorter(r0, th0);
		v_t = t_slope(r0, th0, index, data);
		v_phi = phi_slope(r0, th0, index, data);

		
		x = (sqrt(pow(r0,2) + pow(a,2) )* sin(th0) * cos(phi0));
		y = (sqrt(pow(r0,2) + pow(a,2) )* sin(th0) * sin(phi0));
		z = (r0 * cos(th0));
		size = size +1;

		
	}


	void RKCK(class Stationary_axisymmetric_data data, double stop) {

		int index = 0;
		double tau = 0;
		double dtau = 0.0;

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

		while (tau > tau_final) {
			if (r <= 3.0*M || (th>=-0.15 && th<=0.15) || (th>=M_PI-0.15 && th<=M_PI+0.15)){
				dtau = - pow(10,-3);
			}
			else if (3.0*M < r  && r<=10.0*M){
				dtau = -pow(10,-2);
			}
			else if (10.0*M < r && r<=50.0*M){
				dtau = -5.0*pow(10,-2);
			}
			else{
				dtau = -pow(10,-1);
			}

			index = data.sorter(r, th);

			v_r_k1 = (dtau)*v_r_slope(r, th, v_t, v_r, v_th, v_phi, index, data);
			r_k1 = (dtau)*v_r;
			v_th_k1 = (dtau)*v_th_slope(r, th, v_t, v_r, v_th, v_phi, index, data);
			th_k1 = (dtau)*v_th;
			t_k1 = (dtau)*v_t;
			phi_k1 = (dtau)*v_phi;

			r_d = r + r_k1 / 5.0;


			if (r_d <=r_min){
				inspiral=true;
				break;
			}
			if (r_d >= stop){
				break;
			}
		
			th_d = th + th_k1 / 5.0;
			index = data.sorter(r_d, th_d);
			v_r_d = v_r + v_r_k1 / 5.0;
			v_th_d = v_th + v_th_k1 / 5.0;
			v_t = t_slope(r_d, th_d, index, data);
			v_phi = phi_slope(r_d, th_d, index, data);


			v_r_k2 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			r_k2 = (dtau) * (v_r_d);
			v_th_k2 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			th_k2 = (dtau) * (v_th_d);
			t_k2 = (dtau)*v_t;
			phi_k2 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 40.0) + r_k2 * (9.0 / 40.0);


			if (r_d <=r_min){
				inspiral=true;
				break;
			}
			if (r_d >= stop){
				break;
			}

			th_d = th + th_k1 * (3.0 / 40.0) + th_k2 * (9.0 / 40.0);
			index = data.sorter(r_d, th_d);
			v_r_d = v_r + v_r_k1 * (3.0 / 40.0) + v_r_k2 * (9.0 / 40.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 40.0) + v_th_k2 * (9.0 / 40.0);
			v_t = t_slope(r_d, th_d, index, data);
			v_phi = phi_slope(r_d, th_d, index, data);

			v_r_k3 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			r_k3 = (dtau) * (v_r_d);
			v_th_k3 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			th_k3 = (dtau) * (v_th_d);
			t_k3 = (dtau)*v_t;
			phi_k3 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 10.0) + r_k2 * (-9.0 / 10.0) + r_k3 * (6.0 / 5.0);


			if (r_d <=r_min){
				inspiral=true;
				break;
			}
			if (r_d >= stop){
				break;
			}

			th_d = th + th_k1 * (3.0 / 10.0) + th_k2 * (-9.0 / 10.0) + th_k3 * (6.0 / 5.0);
			index = data.sorter(r_d, th_d);
			v_r_d = v_r + v_r_k1 * (3.0 / 10.0) + v_r_k2 * (-9.0 / 10.0) + v_r_k3 * (6.0 / 5.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 10.0) + v_th_k2 * (-9.0 / 10.0) + v_th_k3 * (6.0 / 5.0);
			v_t = t_slope(r_d, th_d, index, data);
			v_phi = phi_slope(r_d, th_d, index, data);

			v_r_k4 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			r_k4 = (dtau) * (v_r_d);
			v_th_k4 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			th_k4 = (dtau) * (v_th_d);
			t_k4 = (dtau)*v_t;
			phi_k4 = (dtau)*v_phi;

			r_d = r + r_k1 * (-11.0 / 54.0) + r_k2 * (5.0 / 2.0) + r_k3 * (-70.0 / 27.0) + r_k4 * (35.0 / 27.0);


			if (r_d <=r_min){
				inspiral=true;
				break;
			}
			if (r_d >= stop){
				break;
			}

			th_d = th + th_k1 * (-11.0 / 54.0) + th_k2 * (5.0 / 2.0) + th_k3 * (-70.0 / 27.0) + th_k4 * (35.0 / 27.0);
			index = data.sorter(r_d, th_d);
			v_r_d = v_r + v_r_k1 * (-11.0 / 54.0) + v_r_k2 * (5.0 / 2.0) + v_r_k3 * (-70.0 / 27.0) + v_r_k4 * (35.0 / 27.0);
			v_th_d = v_th + v_th_k1 * (-11.0 / 54.0) + v_th_k2 * (5.0 / 2.0) + v_th_k3 * (-70.0 / 27.0) + v_th_k4 * (35.0 / 27.0);
			v_t = t_slope(r_d, th_d, index, data);
			v_phi = phi_slope(r_d, th_d, index, data);

			v_r_k5 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			r_k5 = (dtau) * (v_r_d);
			v_th_k5 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			th_k5 = (dtau) * (v_th_d);
			t_k5 = (dtau)*v_t;
			phi_k5 = (dtau)*v_phi;

			r_d = r + r_k1 * (1631.0 / 55296.0) + r_k2 * (175.0 / 512.0) + r_k3 * (575.0 / 13824.0) + r_k4 * (44275.0 / 110592.0) + r_k5 * (253.0 / 4096.0);


			if (r_d <=r_min){
				inspiral=true;
				break;
			}
			if (r_d >= stop){
				break;
			}

			th_d = th + th_k1 * (1631.0 / 55296.0) + th_k2 * (175.0 / 512.0) + th_k3 * (575.0 / 13824.0) + th_k4 * (44275.0 / 110592.0) + th_k5 * (253.0 / 4096.0);
			index = data.sorter(r_d, th_d);
			v_r_d = v_r + v_r_k1 * (1631.0 / 55296.0) + v_r_k2 * (175.0 / 512.0) + v_r_k3 * (575.0 / 13824.0) + v_r_k4 * (44275.0 / 110592.0) + v_r_k5 * (253.0 / 4096.0);
			v_th_d = v_th + v_th_k1 * (1631.0 / 55296.0) + v_th_k2 * (175.0 / 512.0) + v_th_k3 * (575.0 / 13824.0) + v_th_k4 * (44275.0 / 110592.0) + v_th_k5 * (253.0 / 4096.0);
			v_t = t_slope(r_d, th_d, index, data);
			v_phi = phi_slope(r_d, th_d, index, data);

			v_r_k6 = (dtau)*v_r_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			r_k6 = (dtau) * (v_r_d);
			v_th_k6 = (dtau)*v_th_slope(r_d, th_d, v_t, v_r_d, v_th_d, v_phi, index, data);
			th_k6 = (dtau) * (v_th_d);
			t_k6 = (dtau)*v_t;
			phi_k6 = (dtau)*v_phi;

			r_d = r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6;
			if (r_d <= r_min ) {
				inspiral = true;
				break;
			}

			if (r_d >= stop){
				break;
			}
			th_d = th + (37.0 / 378.0) * th_k1 + (250.0 / 621.0) * th_k3 + (125.0 / 594.0) * th_k4 + (512.0 / 1771.0) * th_k6;
			//Order 6 values/estimate of true value
			t = t + (37.0 / 378.0) * t_k1 + (250.0 / 621.0) * t_k3 + (125.0 / 594.0) * t_k4 + (512.0 / 1771.0) * t_k6;
			r = r_d;
			th = th_d;
			phi = phi + (37.0 / 378.0) * phi_k1 + (250.0 / 621.0) * phi_k3 + (125.0 / 594.0) * phi_k4 + (512.0 / 1771.0) * phi_k6;

			index = data.sorter(r, th);

			v_t = t_slope(r, th, index, data);
			v_phi = phi_slope(r, th, index, data);
			v_r = v_r + (37.0 / 378.0) * v_r_k1 + (250.0 / 621.0) * v_r_k3 + (125.0 / 594.0) * v_r_k4 + (512.0 / 1771.0) * v_r_k6;
			v_th  = v_th + (37.0 / 378.0) * v_th_k1 + (250.0 / 621.0) * v_th_k3 + (125.0 / 594.0) * v_th_k4 + (512.0 / 1771.0) * v_th_k6;

			//converting to cartesian coordinates
			x = (sqrt(pow(r,2) + pow(a,2))* sin(th) * cos(phi));
			y = (sqrt(pow(r,2) + pow(a,2))* sin(th) * sin(phi));
			z = (r * cos(th));

			tau = tau + dtau;
			size = size + 1;
			//energy error
			energy_checker = energy_checker + fabs(pow(v_r, 2.0) + pow(v_th, 2.0) * data.bicubic_spline(r, th, index, 3) / data.bicubic_spline(r, th, index, 2) + V_eff(r, th, index, data));
			carter_constant = carter_constant + fabs(pow(sigma(r,th,a)*v_th,2) + pow(cos(th),2)*(-pow(E*a,2) + pow(L/sin(th),2)));

		}
	}

	int colour_coder() {
        int colour;
		
		

		if (inspiral) {
			colour =0;
		}

		
		else {
			if (y < 0 && z>0) {
				colour=1;
			}
			else if (y > 0 && z > 0) {
				colour=2;
			}
			else if (y < 0 && z < 0) {
				colour=3;
			}
			else if (y > 0 && z < 0) {
				colour=4;
			}

		}
		
		
        return colour;

	}

};

int main() {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	double M = 1.0;
	double J = 0.99;
	double a = J/M;
	double R = M + sqrt(pow(M,2) - pow(a,2));


	int r_size = 1000;
	int theta_size = 100;

	double* r = new double[r_size];
	double* theta = new double[theta_size];

	linspace(R + 0.1, 100*M + 0.1, r_size, r);
	linspace(0.0, M_PI, theta_size, theta);

	//make data class that will create/read data and the splines.
	Stationary_axisymmetric_data data(r_size, theta_size, r[1] - r[0], theta[1] - theta[0], M , a);

	//make the grid of data
	data.make_data(r, theta, M, a);

	delete[] r;
	r = NULL;

	delete[] theta;
	theta = NULL;

	//calculate the spline coefficients
	data.bicubic_spline_coefficient();
    std::cout<<"Created bicubic spline coefficients"<<"\n";

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;
	
	

	//shadow formation
	double* alpha = new double[250];
	double* beta = new double[250];

	double lim =  0.73;
	linspace(-lim, lim, 250, alpha);
	linspace(-lim, lim, 250, beta);

	std::ofstream myFile("kerr_shadow_data_a=0.99_spline.txt");
	Stationary_axisymmetric_photon photon(0.0, 0.0, -1000, data, M, a);
	double x;
	double y;
	double phi = 0.0;
	double th = M_PI_2;
	double E;
	double L;
	double v_r_0;
	double v_th_0;
	double r0 = 15.0*M;
	double mod_P = 1.0;
	double A_t;
	double P_r_obs;
	double P_th_obs;
	double P_phi_obs;
	for (int i = 0; i < 250; i++) {
		for (int j = 0; j < 250; j++) {
			x = - sin(beta[i])*r0;
			y = sin(alpha[j])*r0;
			P_r_obs = mod_P*cos(alpha[j])*cos(beta[i]);
			P_th_obs = mod_P*sin(alpha[j]);
			P_phi_obs = mod_P*cos(alpha[j])*sin(beta[i]);
			A_t = sqrt(g_phiphi(M,r0,th,a)/(pow(g_tphi(M,r0,th,a),2) - g_tt(M,r0,th,a)*g_phiphi(M,r0,th,a)));

			E = (1.0/A_t) - P_phi_obs*g_tphi(M,r0,th,a)/sqrt(g_phiphi(M,r0,th,a));
			L = P_phi_obs*sqrt(g_phiphi(M,r0,th,a));
			v_r_0 = P_r_obs/sqrt(g_rr(M,r0,th,a));
			v_th_0 = P_th_obs/sqrt(g_thth(r0,th,a));
			
			photon.E = E;
			photon.L= L;
			photon.Q = 0.0;
			photon.carter_constant = 0.0;
			photon.energy_checker = 0.0;
			photon.size = 0;
			photon.inspiral = false;
			
			
			//set initial data, check if its in allowed region
			photon.set_initial_data(0.0, r0, th, phi, data);
			photon.v_r = v_r_0;
			photon.v_th = v_th_0;
			//do the RKCK and error 
			photon.Q = pow(sigma(r0,th,a)*v_th_0,2) + pow(cos(th),2)*(-pow(E*a,2) + pow(L/sin(th),2));
			photon.energy_checker = photon.energy_checker + fabs(pow(v_r_0, 2.0) + pow(v_th_0, 2.0) * g_thth(r0,th,a)/g_rr(M,r0,th,a) + V_eff(E,L,M,r0,th,a));
			photon.carter_constant = photon.carter_constant + fabs(photon.Q);
			photon.RKCK(data, 15.0*M);
			myFile << x << "," << y << "," << photon.colour_coder() << ","<<photon.r << ","<<photon.th <<","<<photon.energy_checker/photon.size << ","<< fabs(photon.Q - photon.carter_constant/photon.size)<<"\n";
			//std::cout<< j <<"\n";
		}
		std::cout << i << "\n";
	}
	myFile.close();



}