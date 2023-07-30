#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

class Schwarz_potential_test_photon{
public:
    double M;
    double e;
    double l;
    double H;
    double r_circ;
    double H_circ;
    double taufinal;
    double energy_checker;
    int size;
        
    //  To store r,v,t and phi
    double r;
    double v;
    double t;
    double phi;
    //Arrays that will contain x and y positions
    std::vector <double> x;
    std::vector <double> y;

    Schwarz_potential_test_photon(double mass, double angmom, double tf)
	{
		M = mass;
		l = angmom;
        H = 0.0;
		e = 0.0;
		r_circ = 3.0*mass;
        H_circ = pow(l,2)/(54.00*pow(M,2));
        taufinal = tf;
        energy_checker =0.0;
        size =0;
	}
    double V_eff(double r){
        return pow(l,2)/(2.0*r*r) - M*l*l/(r*r*r);
    }
    
    double v_slope(double r){
        return pow(l,2)/pow(r,3)  - (3.0*M*l*l)/(pow(r,4));
    }

    double t_slope(double r){
        return e/(1 -2.0*M/r );
    }

    double phi_slope(double r){
        return l/(r*r);
    }

    void V_eff_data() {
		std::ofstream myFile("schwarz_V_eff.txt");
		for (int r = int(2.0*M); r < 101; r++) {
			myFile << r << "," << V_eff(r) << "\n";
		}
        myFile.close();
	}

    void set_initial_data(double hamiltonian, double t0, double r0, double phi0){
        H= hamiltonian;
        e = sqrt(2.0*hamiltonian);
        //for a circular orbit
        v = (0.0); 
        //for a photon going in first
        //v.push_back(-sqrt(2.0*(H - V_eff(r0))));
        //for a photon going out
        //v.push_back(+sqrt(2.0*(H - V_eff(r0)))); 
        t = (t0);
        r = (r0);
        phi = (phi0);
        x.push_back(r0*cos(phi0));
        y.push_back(r0*sin(phi0));

        energy_checker = energy_checker + fabs( H - (1.0/2.0)*v*v - V_eff(r0));
        size =1;
    }
   
    void RKCK(){
        double tau = 0.0;
        double dtau = 0.0;

        double v_k1;
        double r_k1;
        double t_k1;
        double phi_k1;

        double v_k2;
        double r_k2;
        double t_k2;
        double phi_k2;

        double v_k3;
        double r_k3;
        double t_k3;
        double phi_k3;

        double v_k4;
        double r_k4;
        double t_k4;
        double phi_k4;

        double v_k5;
        double r_k5;
        double t_k5;
        double phi_k5;

        double v_k6;
        double r_k6;
        double t_k6;
        double phi_k6;

        double r_d;
        double v_d;


        while(tau<taufinal){
            //start with (r[i],v[i]) and employ RKCK algorithm to find the 4 slopes

            //change the step size for the current step
            if (r<=3.0*M){
                dtau = pow(10,-3);
            }
            else if (r>3.0*M && r<=10.0*M){
                dtau = pow(10,-2);
            }
            else if (r>10.0*M && r<=50*M){
                dtau = 5.0*pow(10,-2);
            }
            else {
                dtau = pow(10,-1);
            }

            v_k1 = (dtau)*v_slope(r);
            r_k1 = (dtau)*v;
            t_k1 = (dtau)*t_slope(r);
            phi_k1 = (dtau)*phi_slope(r);

            r_d = r + r_k1/5.0;
            if (r_d<=2.0*M){
                break;
            }
            v_d = v + v_k1/5.0;

            v_k2 = (dtau)*v_slope(r_d);
            r_k2 = (dtau)*(v_d);
            t_k2 = (dtau)*t_slope(r_d);
            phi_k2 = (dtau)*phi_slope(r_d);

            r_d = r + r_k1 * (3.0 / 40.0) + r_k2 * (9.0 / 40.0);
            if (r_d<=2.0*M){
                break;
            }
            v_d = v + v_k1 * (3.0 / 40.0) + v_k2 * (9.0 / 40.0);

            v_k3 = (dtau)*v_slope(r_d);
            r_k3 = (dtau)*(v_d);
            t_k3 = (dtau)*t_slope(r_d);
            phi_k3 = (dtau)*phi_slope(r_d);

            r_d = r + r_k1 * (3.0 / 10.0) + r_k2 * (-9.0 / 10.0) + r_k3 * (6.0 / 5.0);
            if (r_d<=2.0*M){
                break;
            }
            v_d = v + v_k1 * (3.0 / 10.0) + v_k2 * (-9.0 / 10.0) + v_k3 * (6.0 / 5.0);

            v_k4 = (dtau)*v_slope(r_d);
            r_k4 = (dtau)*(v_d);
            t_k4 = (dtau)*t_slope(r_d);
            phi_k4 = (dtau)*phi_slope(r_d);

            r_d = r + r_k1 * (-11.0 / 54.0) + r_k2 * (5.0 / 2.0) + r_k3 * (-70.0 / 27.0) + r_k4 * (35.0 / 27.0);
			if (r_d <= 2.0*M) {
				break;
			}
            v_d = v + v_k1 * (-11.0 / 54.0) + v_k2 * (5.0 / 2.0) + v_k3 * (-70.0 / 27.0) + v_k4 * (35.0 / 27.0);

            v_k5 = (dtau)*v_slope(r_d);
            r_k5 = (dtau)*(v_d);
            t_k5 = (dtau)*t_slope(r_d);
            phi_k5 = (dtau)*phi_slope(r_d);

            r_d = r + r_k1 * (1631.0 / 55296.0) + r_k2 * (175.0 / 512.0) + r_k3 * (575.0 / 13824.0) + r_k4 * (44275.0 / 110592.0) + r_k5 * (253.0 / 4096.0);

			if (r_d <= 2.0*M) {
				break;
			}
            v_d = v + v_k1 * (1631.0 / 55296.0) + v_k2 * (175.0 / 512.0) + v_k3 * (575.0 / 13824.0) + v_k4 * (44275.0 / 110592.0) + v_k5 * (253.0 / 4096.0);
            
            v_k6 = (dtau)*v_slope(r_d);
            r_k6 = (dtau)*(v_d);
            t_k6 = (dtau)*t_slope(r_d);
            phi_k6 = (dtau)*phi_slope(r_d);
            
            r_d = r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6;

            if (r_d <= 2.0*M) {
				break;
			}

            t = (t + (37.0 / 378.0) * t_k1 + (250.0 / 621.0) * t_k3 + (125.0 / 594.0) * t_k4 + (512.0 / 1771.0) * t_k6);
            r = r_d;
            phi = (phi + (37.0 / 378.0) * phi_k1 + (250.0 / 621.0) * phi_k3 + (125.0 / 594.0) * phi_k4 + (512.0 / 1771.0) * phi_k6);
            v = (v + (37.0 / 378.0) * v_k1 + (250.0 / 621.0) * v_k3 + (125.0 / 594.0) * v_k4 + (512.0 / 1771.0) * v_k6);
            //Once order 6 confirms this time step as valid, increment data
            size = size+1;
            energy_checker = energy_checker + fabs( H- (1.0/2.0)*v*v - V_eff(r));
            x.push_back(r*cos(phi));
            y.push_back(r*sin(phi));
            tau = tau + dtau;
            
        }
     
    }

    void save_data() {
		std::ofstream myFile("photon.txt");
		for (int i = 0; i < size; i++) {
			myFile << x[i] << "," << y[i] << "\n";
		}
        myFile.close();
	}


};

int main() {

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	double M = 1.0;
    double l = 1.0;
	
    //Declare instance of a photon in the schwarzschild potential
	Schwarz_potential_test_photon photon(M, l, 500.0);
	
    //First plot the potential and analyze it
    //photon.V_eff_data();

    //Set the initial data for desired kind of orbit
    //this is for a circular orbit
    photon.set_initial_data(photon.H_circ,0.0,photon.r_circ,M_PI_4);

    //photon.set_initial_data(0.99*(1.0/54.0),0.0,100.0,atan2(sqrt(27.0/0.99), 100.0));

    //Do the RKCK
    photon.RKCK();
    std::cout << "The size of the iteration was: "<< photon.size<<"\n";
    std::cout << "The average deviation from zero was: "<< (photon.energy_checker/photon.size) <<"\n";

	//save trajectory
	photon.save_data();
	
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;

}