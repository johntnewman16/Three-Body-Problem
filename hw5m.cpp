/*
AEP 4380 Homework 5

Three Body Problem

Jack Newman October 14, 2015

Last Edited: October 23, 2015

Intel(R) Core(TM) i7-4810MQ CPU @ 2.80 GHz

MinGW C++ compiler, GCC 4.8.1

Note: Comments are either inline or below the lines they describe
*/

#include <cstdlib> // importing plain c
#include <cmath> // math library

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output

using namespace std;

double w = 1.0;
bool collision;
double const msd2 = 7.46496E9, msd = 86400;
double const M[] = {5.976E24, 0.0123*5.976E24, 0.0123*5.976E24*0.2};//0.0213*5.976E24*0.2}; // mass constants in kg
double const Tm = 648, G = 6.6726E-11; // orbital period moon, gravitational const
double const R[] = {6378000.0, 3476000.0, 0.5*3476000.0}; // radius constants in m
double const c2=0.2,c3=0.3,c4=0.8,c5=8.0/9.0,a21=0.2,a31=3.0/40.0,
a32=9.0/40.0,a41=44.0/45.0,a42=-56.0/15.0,a43=32.0/9.0,a51=19372.0/6561.0,
a52=-25360.0/2187.0,a53=64448.0/6561.0,a54=-212.0/729.0,a61=9017.0/3168.0,
a62=-355.0/33.0,a63=46732.0/5247.0,a64=49.0/176.0,a65=-5103.0/18656.0,
a71=35.0/384.0, a72=0.0, a73=500.0/1113.0,a74=125.0/192.0,a75=-2187.0/6784.0,
a76=11.0/84.0,e1=5179.0/57600.0,e2=0.0,e3=7571.0/16695.0,e4=393.0/640.0,
e5=-92097.0/339200.0,e6=187.0/2100.0,e7=1.0/40.0;
// establishing the constants and parameters used later in the program
const int N1 = 3; // number of objects
const int N = 2*2*N1;
//const int N = 2*2*N1; // number of ODE
//evaluate the ODE right hand side
//change this for each different ODE


void myrhs( double y[], double t, double f[]){
	int istep, j;//right hand side of the Rossier system.
	//bool collision = false;
	double d, dx, dy;
	//cout << "call myrhs()\n" << endl;
	for(istep=0;istep<N1; istep++){
		//y[i] x pos. of ith object 
		//y[i+N1] y pos. of ith object
		//y[i+2*N1] vx velocity of ith object vxo
		//y[i+3*N1] vy velocity of ith object
		f[istep] = y[istep+2*N1]; // dx/dt
		f[istep + N1] = y[istep+3*N1]; // dy/dt
		f[istep+2*N1] = 0;
		f[istep+3*N1] = 0;
		for(j=0; j<N1; j++){ // + dy^2
			if( j != istep ){
				dx = (y[j]-y[istep]);
				dy = (y[j+N1]-y[istep+N1]);
				d = sqrt(dx*dx + dy*dy);
				f[istep+2*N1] = f[istep+2*N1] + G*M[j]*dx/(d*d*d); // dvx/dt
				f[istep+3*N1] = f[istep+3*N1] + G*M[j]*dy/(d*d*d); //dvy/dt
				//cout << "y[4] =" << y[4] << endl;
				if (d < (R[istep] + R[j])){
					collision = true; //register if there is a collision
				}
			}
		}
	}
	return;
}

void myrhs1( double y[], double t, double f[]){
	//right hand side of the simple harmonic oscillator
	f[1] = -w*w*y[0]; // f[1] = dy_1/dt
	f[0] = y[1]; // f[0] = dy_o/dt
	return;
}



void ode45(double yold[], double ynew[], double ynew2[], double yerr[], double h, 
	 double t, int neqns, void myrhs(double[], double, double[])){
	// ODE solver using 5th/4th Order Runge-Kutta AutoStep Method
	// yold[] is the array which stores the initial y-values of the equation
	// ynew[] stores the new y-values once the code has completed
	// h is the step size
	// t is the initial t value put into the ODE solver
	// neqns is the number of equations that are needed to describe the system
	// myrhs are the equations which define the system

	// ODE solver using 4th Order Runge-Kutta
	int i; // represents a counter for the number of equations
	double *k1, *k2, *k3, *k4, *k5, *k6, *k7, *temp;
	//initializing the 8 arrays needed for the 5th order method

	k1 = new double[8*neqns]; 
	// dynamically allocate 8*neqns memory locations for further use
	if(NULL == k1){ //check that array is allocated
		cout << "Cannot allocate k1 in ode4" << endl;
		exit(0);
	}
	k2 = k1 + neqns;
	k3 = k2 + neqns;
	k4 = k3 + neqns;
	k5 = k4 + neqns;
	k6 = k5 + neqns;
	k7 = k6 + neqns;	
	temp = k7 + neqns;

	myrhs(yold, t, k1); 
	// first, evaluate the rhs equations at t, and then store the values in k1

	for(i = 0; i < neqns; i++){
		temp[i] = yold[i] + a21*h*k1[i];
	}
	myrhs(temp, t + c2*h, k2);
	for(i = 0; i < neqns; i++){
		temp[i] = yold[i] + h*(a31*k1[i] + a32*k2[i]);

	}
	myrhs(temp, t + c3*h, k3);
	for(i = 0; i < neqns; i++){
		temp[i] = yold[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
	}
	myrhs(temp, t + c4*h, k4);
	for(i = 0; i < neqns; i++){
		temp[i] = yold[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
	}
	myrhs(temp, t + c5*h, k5);
	for(i = 0; i < neqns; i++){
		temp[i] = yold[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
	}
	myrhs(temp, t + h, k6);
	for(i = 0; i < neqns; i++){
		ynew[i] = yold[i] + h*(a71*k1[i] + a72*k2[i] + a73*k3[i] + 
			a74*k4[i] + a75*k5[i] + a76*k6[i]);

		// Use the calculated k-values to find the new ynew values.
	}
	myrhs(ynew, t + h, k7);
	for(i = 0; i < neqns; i++){
		ynew2[i] = yold[i] + h*(e1*k1[i] + e2*k2[i] + e3*k3[i] + 
			e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k7[i]);
		//calculate ynew and k-values to find ynew2, the 4th order solution
	}
	
	for(i = 0; i < neqns; i++){
		yerr[i] = abs(ynew[i]-ynew2[i]);
		//find the difference between the 5th and 4th order solution
	}
	// The above steps represent the 5th/4th Order Runge-Kutta formulas as
	// described in Numerical Recipes

	delete[] k1; // clear the memory so that new values can be stored.
	return;
}//end ode45

double maxi(double *arr, int neqns){
	double scale[] = {6378, 3E8, 3E8, 6378, 3E8, 3E8, 10, 500, 500, 10, 500, 500};
	//scale used to create a reasonable relative error for each of the values
	int i;
	double max = arr[0]/scale[0];//scale the difference
	for(i = 1; i < neqns; i++){
        if(max < arr[i]/scale[i]){
        	max = arr[i]/scale[i];
        }//sort through the arrays and get the maximum error/scale
    }
    return max;
}



int main()
{	
	ofstream fp;

	fp.open( "example5odea.dat" ); //open new file for graph
	if( fp.fail() ) { // or fp.bad()
		cout << "cannot open file" << endl;
		// catches fp.open failure, and returns appropriate error
		return( EXIT_SUCCESS );
	}

	int i; // counter
	double yold[N], ynew[N], yerr[N], ynew2[N];
	double tmin = 0, h=1000.0, tol = .00000001, tmax=200.0*msd*5, t, max;
	//set the initial values of the equation

	yold[0] = 0.0; yold[1] = 0.0; yold[2] = -6.0E8;
	yold[3] = 0.0; yold[4] = 3.84E8; yold[5] = 4.05E8;
	yold[6] = -12.593; yold[7] = 1019.000; yold[8] = 800.0000;
	yold[9] = 0.0; yold[10] = 0.0; yold[11] = 500.0000;
	// initialize yold[]
	t = tmin;
	while(t<tmax){//until we reach tmax
		ode45( yold, ynew, ynew2, yerr, h, t, 
			N, myrhs);
		max = maxi(yerr, N);
		//cout << yold[2] << endl;
		//exit(0);
		cout << "max = " << max << endl;
		//implement the ODE solver at each step

		if(max < tol){ //step succeeded
			cout << "STEP SUCCESS" << endl;;
    		h = h*pow(tol/(max), 0.2);
    		//if the step works, adjust h according to the ratio of tol to max
			cout << "h =" << h << endl;
    		t = t + h;
			cout << "t =" << t << endl;
    		for(i = 0; i < N; i++){
				yold[i] = ynew[i]; 
				//Set values of yold to ynew, iterating to the next step

				fp << h << setw(15) << t/msd << setw(15) << ynew[0] << setw(15) << 
				ynew[1] << setw(15) << ynew[2] << setw(15) << ynew[3] << 
				setw(15) << ynew[4] << setw(15) << ynew[5] << setw(15) << 
				ynew[6] << setw(15) << ynew[7] << setw(15) << 
				ynew[8] << setw(15) << ynew[9] << setw(15) << 
				ynew[10] << setw(15) << ynew[11] <<
				endl;
				//store values in array for plotting and tables
			}
    	}

    	else{// if maximum error is greater than the tolerance
			cout << "STEP FAILURE, max = " << max << endl;
    		h /= 5.0; //reduce the step size of h

    		if(abs(h) < 1E-20){
    			cout << "h too small, abort program" << endl;
    			return( EXIT_SUCCESS );
    			//if h is trending to 0, abort the program, and return error message
			}
		}

		if(collision == true){//if the radius of the objects come within d of each other
			cout << "collision occurred at t =" << t/msd << endl;
			fp.close();
			exit(0);
			//print error message, indicating the time of the collision and exit.

		}
		cout << setw(15) << t/msd << setw(15) << ynew[9] << setw(15) << ynew[10] << 
		setw(15) << ynew[11] << endl;
	}
	fp.close();

	return(EXIT_SUCCESS);
}



