#define _USE_MATH_DEFINES
#include <cmath>
#include<iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>

//this is a function to create linearly separated points. 
void linspace(double start, double end, double numpoints, double* points) {
	double dx = (end - start) / (numpoints - 1);
	for (int i = 0; i < numpoints; i++) {
		points[i] = (start + i * dx);
	}
}

//these are metric functions for Kerr
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
	return sigma(r, th, a) / delta(M, r, a);
}

double g_thth(double r, double th, double a) {
	return sigma(r, th, a);
}

double g_phiphi(double M, double r, double th, double a) {
	return (pow(a, 2) + pow(r, 2) + (2 * M * r * pow(a * sin(th), 2.0)) / (sigma(r, th, a))) * pow(sin(th), 2);
}

//Slopes for t, phi and the effective_potential equation

double t_slope(double E, double L, double M, double r, double th, double a) {
	return (E * g_phiphi(M, r, th, a) + L * g_tphi(M, r, th, a)) / (pow(g_tphi(M, r, th, a), 2) - g_tt(M, r, th, a) * g_phiphi(M, r, th, a));
}

double phi_slope(double E, double L, double M, double r, double th, double a) {
	return -(E * g_tphi(M, r, th, a) + L * g_tt(M, r, th, a)) / (pow(g_tphi(M, r, th, a), 2) - g_tt(M, r, th, a) * g_phiphi(M, r, th, a));
}

double V_eff(double E, double L, double M, double r, double th, double a) {
	return (pow(g_rr(M, r, th, a), -1)) * ( 1.0+ (g_phiphi(M, r, th, a) * pow(E, 2) + g_tt(M, r, th, a) * pow(L, 2) + 2 * E * L * g_tphi(M, r, th, a)) / (-pow(g_tphi(M, r, th, a), 2) + g_tt(M, r, th, a) * g_phiphi(M, r, th, a)));
}

// all the derivatives of the metric required are also coded as functions

//first is for the radial equation of motion
double g_rr_th(double M, double r, double th, double a) {
	return -(sin(2 * th) * pow(a, 2)) / (delta(M, r, a));
}

double g_rr_r(double M, double r, double th, double a) {
	return (2 / pow(delta(M, r, a), 2)) * (r * delta(M, r, a) - (r - M) * sigma(r, th, a));

}

double g_tt_r(double M, double r, double th, double a) {
	return ((2 * M) / (pow(sigma(r, th, a), 2))) * (sigma(r, th, a) - 2 * pow(r, 2));
}

double g_tphi_r(double M, double r, double th, double a) {
	return (2 * M * a * pow(sin(th), 2) / pow(sigma(r, th, a), 2)) * (2 * pow(r, 2) - sigma(r, th, a));
}

double g_thth_r(double r) {
	return 2 * r;
}

double g_phiphi_r(double M, double r, double th, double a) {
	return (pow(sin(th), 2)) * (2 * r + ((2 * M * pow(a * sin(th), 2)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) - 2 * pow(r, 2)));
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
	return (sin(2 * th) / pow(sin(th), 2)) * g_phiphi(M, r, th, a) + ((2 * M * r * sin(2 * th) * pow(a * sin(th), 2)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) + pow(a * sin(th), 2));
}

double v_r_slope(double M, double r, double th, double a, double v_t, double v_r, double v_th, double v_phi) {
	return  (-1 / (2 * g_rr(M, r, th, a))) * (2 * v_r * v_th * g_rr_th(M, r, th, a) + pow(v_r, 2) * g_rr_r(M, r, th, a) - pow(v_t, 2) * g_tt_r(M, r, th, a) - 2 * v_t * v_phi * g_tphi_r(M, r, th, a) - pow(v_th, 2) * g_thth_r(r) - pow(v_phi, 2) * g_phiphi_r(M, r, th, a));
}

double v_th_slope(double M, double r, double th, double a, double v_t, double v_r, double v_th, double v_phi) {
	return (-1 / (2 * g_thth(r, th, a))) * (2 * v_r * v_th * g_thth_r(r) - pow(v_r, 2) * g_rr_th(M, r, th, a) - pow(v_t, 2) * g_tt_th(M, r, th, a) - 2 * v_t * v_phi * g_tphi_th(M, r, th, a) + pow(v_th, 2) * g_thth_th(th, a) - pow(v_phi, 2) * g_phiphi_th(M, r, th, a));
}

class Kerr_potential_test_photon {

public:
	double mass; //here mass is the mass of the black hole
	double angmom; //here angmom is the angmom per unit mass(a) of the black hole
	//the test particle has unit mass
	double E; //E is the conserved energy set later on
	double L; //L is the conserved angular momentum set later on
	double Q; //Q is the carter constant found by setting a valie pf v_th later on
	double r_event_min; //inner event horizon
	double r_event_max; //outer event horizon
	double taufinal;
	double energy_checker;
	double carter_constant;
	int size;
	
	//Arrays to store data points from RK-4 integration

	double t;
	double r;
	double th;
	double phi;
	double v_t;
	double v_r;
	double v_th;
	double v_phi;

	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;


	Kerr_potential_test_photon(double M, double a, double En, double Lz, double tauend)
	{
		mass = M;
		angmom = a;
		E = En;
		L = Lz;
		Q =0.0;
		r_event_min = mass - sqrt(pow(mass, 2) - pow(angmom, 2));
		r_event_max = mass + sqrt(pow(mass, 2) - pow(angmom, 2));
		taufinal = tauend;
		energy_checker = 0.0;
		carter_constant =0.0;
		size = 0;
		
	}

	void set_initial_data(double t0, double r0, double th0, double phi0) {
		t = (t0);
		r = (r0);
		th = (th0);
		phi = (phi0);

		v_t = t_slope(E, L, mass, r0, th0, angmom);
		v_phi = phi_slope(E,L,mass,r0,th0,angmom);
		//this gives an outward orbit
		//v_th = 0.0;
		//v_r = + sqrt(-V_eff(E,L,mass,r0,th0,angmom));
		//this gives an inward orbit
		v_th = 0.0;
		v_r = - sqrt(-V_eff(E,L,mass,r0,th0,angmom));

		Q = pow(sigma(r0,th0,angmom)*v_th,2) + pow(cos(th0),2)*(pow(angmom,2)*(1.0 - pow(E,2))+ pow(L/sin(th0),2));
		
		x.push_back(sqrt(pow(r0,2) + pow(angmom,2)) * sin(th0) * cos(phi0));
		y.push_back(sqrt(pow(r0,2) + pow(angmom,2)) * sin(th0) * sin(phi0));
		z.push_back(r0 * cos(th0));

		energy_checker = energy_checker + fabs(pow(v_r, 2.0) + pow(v_th, 2.0) * g_thth(r0,th0,angmom) / g_rr(mass,r0,th0,angmom) + V_eff(E,L,mass,r0,th0,angmom));
		carter_constant = carter_constant + fabs(Q);
		size =1;
	}

	void RKCK() {
		double tau = 0;
		double dtau = 0;

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
		
		
		while (tau < taufinal) {

			if (r<=3.0*mass || (th>=-0.15 && th<=0.15) || (th>=M_PI-0.15 && th<=M_PI+0.15)){
				dtau = pow(10,-3);
			}
			else if (r>3.0*mass && r<=10.0*mass){
                dtau = pow(10,-2);
            }
            else if (r>10.0*mass && r<=50*mass){
                dtau = 5.0*pow(10,-2);
            }
            else {
                dtau = pow(10,-1);
            }

			v_r_k1 = (dtau)*v_r_slope(mass,r, th, angmom, v_t, v_r, v_th, v_phi);
			r_k1 = (dtau)*v_r;
			v_th_k1 = (dtau)*v_th_slope(mass,r, th, angmom, v_t, v_r, v_th, v_phi);
			th_k1 = (dtau)*v_th;
			t_k1 = (dtau)*v_t;
			phi_k1 = (dtau)*v_phi;

			r_d = r + r_k1 / 5.0;


			if (r_d <= r_event_max) {
				break;
			}

			th_d = th + th_k1 / 5.0;
			
			v_r_d = v_r + v_r_k1 / 5.0;
			v_th_d = v_th + v_th_k1 / 5.0;
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d,angmom);


			v_r_k2 = (dtau)*v_r_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k2 = (dtau) * (v_r_d);
			v_th_k2 = (dtau)*v_th_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			th_k2 = (dtau) * (v_th_d);
			t_k2 = (dtau)*v_t;
			phi_k2 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 40.0) + r_k2 * (9.0 / 40.0);


			if (r_d <= r_event_max) {
				break;
			}

			th_d = th + th_k1 * (3.0 / 40.0) + th_k2 * (9.0 / 40.0);
			
			v_r_d = v_r + v_r_k1 * (3.0 / 40.0) + v_r_k2 * (9.0 / 40.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 40.0) + v_th_k2 * (9.0 / 40.0);
			v_t = t_slope(E,L,mass,r_d, th_d,angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k3 = (dtau)*v_r_slope(mass,r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k3 = (dtau) * (v_r_d);
			v_th_k3 = (dtau)*v_th_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			th_k3 = (dtau) * (v_th_d);
			t_k3 = (dtau)*v_t;
			phi_k3 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 10.0) + r_k2 * (-9.0 / 10.0) + r_k3 * (6.0 / 5.0);


			if (r_d <= r_event_max) {
				break;
			}

			th_d = th + th_k1 * (3.0 / 10.0) + th_k2 * (-9.0 / 10.0) + th_k3 * (6.0 / 5.0);
			v_r_d = v_r + v_r_k1 * (3.0 / 10.0) + v_r_k2 * (-9.0 / 10.0) + v_r_k3 * (6.0 / 5.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 10.0) + v_th_k2 * (-9.0 / 10.0) + v_th_k3 * (6.0 / 5.0);
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k4 = (dtau)*v_r_slope(mass,r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k4 = (dtau) * (v_r_d);
			v_th_k4 = (dtau)*v_th_slope(mass,r_d, th_d, angmom,v_t, v_r_d, v_th_d, v_phi);
			th_k4 = (dtau) * (v_th_d);
			t_k4 = (dtau)*v_t;
			phi_k4 = (dtau)*v_phi;

			r_d = r + r_k1 * (-11.0 / 54.0) + r_k2 * (5.0 / 2.0) + r_k3 * (-70.0 / 27.0) + r_k4 * (35.0 / 27.0);


			if (r_d <= r_event_max) {
				break;
			}

			th_d = th + th_k1 * (-11.0 / 54.0) + th_k2 * (5.0 / 2.0) + th_k3 * (-70.0 / 27.0) + th_k4 * (35.0 / 27.0);
			v_r_d = v_r + v_r_k1 * (-11.0 / 54.0) + v_r_k2 * (5.0 / 2.0) + v_r_k3 * (-70.0 / 27.0) + v_r_k4 * (35.0 / 27.0);
			v_th_d = v_th + v_th_k1 * (-11.0 / 54.0) + v_th_k2 * (5.0 / 2.0) + v_th_k3 * (-70.0 / 27.0) + v_th_k4 * (35.0 / 27.0);
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k5 = (dtau)*v_r_slope(mass, r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k5 = (dtau) * (v_r_d);
			v_th_k5 = (dtau)*v_th_slope(mass,r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			th_k5 = (dtau) * (v_th_d);
			t_k5 = (dtau)*v_t;
			phi_k5 = (dtau)*v_phi;

			r_d = r + r_k1 * (1631.0 / 55296.0) + r_k2 * (175.0 / 512.0) + r_k3 * (575.0 / 13824.0) + r_k4 * (44275.0 / 110592.0) + r_k5 * (253.0 / 4096.0);


			if (r_d <= r_event_max) {
				break;
			}

			th_d = th + th_k1 * (1631.0 / 55296.0) + th_k2 * (175.0 / 512.0) + th_k3 * (575.0 / 13824.0) + th_k4 * (44275.0 / 110592.0) + th_k5 * (253.0 / 4096.0);
			v_r_d = v_r + v_r_k1 * (1631.0 / 55296.0) + v_r_k2 * (175.0 / 512.0) + v_r_k3 * (575.0 / 13824.0) + v_r_k4 * (44275.0 / 110592.0) + v_r_k5 * (253.0 / 4096.0);
			v_th_d = v_th + v_th_k1 * (1631.0 / 55296.0) + v_th_k2 * (175.0 / 512.0) + v_th_k3 * (575.0 / 13824.0) + v_th_k4 * (44275.0 / 110592.0) + v_th_k5 * (253.0 / 4096.0);
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k6 = (dtau)*v_r_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k6 = (dtau) * (v_r_d);
			v_th_k6 = (dtau)*v_th_slope(mass,r_d, th_d, angmom,v_t, v_r_d, v_th_d, v_phi);
			th_k6 = (dtau) * (v_th_d);
			t_k6 = (dtau)*v_t;
			phi_k6 = (dtau)*v_phi;

	
			//Order 6 values
			
			r = r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6;
			th = th + (37.0 / 378.0) * th_k1 + (250.0 / 621.0) * th_k3 + (125.0 / 594.0) * th_k4 + (512.0 / 1771.0) * th_k6;
			t = t + (37.0 / 378.0) * t_k1 + (250.0 / 621.0) * t_k3 + (125.0 / 594.0) * t_k4 + (512.0 / 1771.0) * t_k6;
			phi = phi + (37.0 / 378.0) * phi_k1 + (250.0 / 621.0) * phi_k3 + (125.0 / 594.0) * phi_k4 + (512.0 / 1771.0) * phi_k6;
			
			v_t = t_slope(E,L,mass,r,th,angmom);
			v_phi = phi_slope(E,L,mass,r,th,angmom);
			v_r = v_r + (37.0 / 378.0) * v_r_k1 + (250.0 / 621.0) * v_r_k3 + (125.0 / 594.0) * v_r_k4 + (512.0 / 1771.0) * v_r_k6;
			v_th = v_th + (37.0 / 378.0) * v_th_k1 + (250.0 / 621.0) * v_th_k3 + (125.0 / 594.0) * v_th_k4 + (512.0 / 1771.0) * v_th_k6; 
			
			
			//converting to cartesian coordinates
			x.push_back(sqrt(pow(r,2) + pow(angmom,2))* sin(th) * cos(phi));
			y.push_back(sqrt(pow(r,2) + pow(angmom,2))* sin(th) * sin(phi));
			z.push_back(r * cos(th));

			tau = tau + dtau;
			size = size + 1;

			//energy error
			energy_checker = energy_checker + fabs(pow(v_r, 2.0) + pow(v_th, 2.0) * g_thth(r, th, angmom) / g_rr(mass, r, th, angmom) + V_eff(E, L, mass, r, th, angmom));
			carter_constant = carter_constant + fabs(pow(sigma(r,th,angmom)*v_th,2) + pow(cos(th),2)*(-pow(E*angmom,2) + pow(L/sin(th),2)));
			

		}
		
	}
		
	void save_to_file() {
		std::ofstream myFile("particle.txt");

		for (int i = 0; i < size; i++) {
			myFile << x[i] << "," << y[i] << "," << z[i] << "\n";
		}
		myFile.close();

	}
};

int main(){

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	

	double M = 1.0;
	double a = 0.5;
	
    //declare a kerr photon
	Kerr_potential_test_photon particle(M, a,1.0,-2.0, 100);
	
    //set initial data
    //this is an ingoing orbit
	particle.set_initial_data(0.0,10.00,M_PI_4, M_PI_4);
		
	//Do the RKCK
	particle.RKCK();
    std::cout<< "The average deviation from expected zero is: "<< particle.energy_checker/particle.size <<"\n";
	std::cout<< "The starting value of the Carter Constant was: "<< particle.Q <<"\n";
	std::cout<< "The numerical value of the average Carter Constant was: "<< particle.carter_constant/particle.size <<"\n";
	std::cout<< "The deviation in the Carter Constant was: "<< fabs(particle.Q - particle.carter_constant/particle.size) <<"\n";
	std::cout<< "The size of the iteration was: "<< particle.size << "\n";
    std::cout <<"The ending radial position is: "<<particle.r<<"\n";
	
	//save to file
	particle.save_to_file();
	
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;

}