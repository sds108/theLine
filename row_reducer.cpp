#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>

using namespace std;

#define PI 3.14159265

// variables

// Distances
long double d0 = 1000; // radius to point of turning
long double d1 = 0.3;
long double d2 = 1.5;
long double d3 = 1.5;
long double d4 = 0;
long double d5 = 0;
long double d6 = 1.5;
long double d7 = 1.5;
long double d8 = 1;
long double d9 = 1;
long double d10 = 0;

// Radii
long double R_COM = 0;
long double R_COP = 0;
long double R_LF = 0;
long double R_RF = 0;
long double R_LR = 0;
long double R_RR = 0;
long double R_Front = 0;
long double R_Rear = 0;

long double R_RR_COM = 0;
long double R_RR_COP = 0;
long double R_RR_RF = 0;

long double R_LR_COM = 0;
long double R_LR_COP = 0;
long double R_LR_LF = 0;

// Angles
long double A_COM = 0;
long double A_COP = 0;
long double A_LF = 0;
long double A_RF = 0;
long double A_LR = 0;
long double A_RR = 0;
long double A_Front = 0;
long double A_Rear = 0;

long double A_LF_RR = 0;
long double A_RR_COM = 0;
long double A_RR_COP = 0;
long double A_RR_RF = 0;

long double A_RF_LR = 0;
long double A_LR_COM = 0;
long double A_LR_COP = 0;
long double A_LR_LF = 0;

// Unknowns
long double N_LF = 0;
long double N_RF = 0;
long double N_LR = 0;
long double N_RR = 0;
long double omega2 = 0;
long double omega = 0;
long double linear_velocity = 0;

// Inputs
long double s = 0.001;
long double other_velocity = 0;
long double m = 900;
long double g = 9.81;
long double mew = 3.5;
long double I_COM = 0;

// Aero
long double C_Lift = 1;
long double C_Drag = 0;
long double A_Lift = 8;
long double A_Drag = 6;
long double rho = 1.225;

// Coefficients
long double e0[6] = {0, 0, 0, 0, 0, 0};
long double e1[6] = {0, 0, 0, 0, 0, 0};
long double e2[6] = {0, 0, 0, 0, 0, 0};
long double e3[6] = {0, 0, 0, 0, 0, 0};
long double e4[6] = {0, 0, 0, 0, 0, 0};


int main() {
	// Variable Assigning
	R_COM = sqrt(pow(d0, 2) + pow(d3, 2));
	A_COM = atan2(d3, d0);
	
	R_COP = sqrt(pow(d0 - d10, 2) + pow(d3 + d4, 2));
	A_COP = atan2(d3 + d4, d0 - d10);
	
	R_LF = sqrt(pow(d0 + d9, 2) + pow(d3 + d2, 2));
	A_LF = atan2(d3 + d2, d0 + d9);
	
	R_RF = sqrt(pow(d0 - d8, 2) + pow(d3 + d2, 2));
	A_RF = atan2(d3 + d2, d0 - d8);
	
	R_LR = d0 + d7;
	
	R_RR = d0 - d6;
	
	R_Front = sqrt(pow(d2, 2) + pow(d1, 2));
	A_Front = atan2(d1, d2);
	
	R_Rear = sqrt(pow(d3, 2) + pow(d1, 2));
	A_Rear = atan2(d1, d3);
	
	A_LF_RR = atan2(d2 + d3, d6 + d9);
	A_RR_COM = atan2(d3, d6);
	A_RR_COP = atan2(d3 + d4, d6 - d10);
	A_RR_RF = atan2(d2 + d3, d6 - d8);
	
	R_RR_COM = sqrt(pow(d3, 2) + pow(d6, 2));
	R_RR_COP = sqrt(pow(d3 + d4, 2) + pow(d6 - d10, 2));
	R_RR_RF = sqrt(pow(d2 + d3, 2) + pow(d6 - d8, 2));
	
	A_RF_LR = atan2(d2 + d3, d7 + d8);
	A_LR_COM = atan2(d3, d7);
	A_LR_COP = atan2(d3 + d4, d7 + d10);
	A_LR_LF = atan2(d2 + d3, d6 - d8);
	
	R_LR_COM = sqrt(pow(d3, 2) + pow(d7, 2));
	R_LR_COP = sqrt(pow(d3 + d4, 2) + pow(d7 + d10, 2));
	R_LR_LF = sqrt(pow(d2 + d3, 2) + pow(d7 - d9, 2));
	
	
	// Gaussian Elimination
    int i,j,k,n;
    
    //cout << "\nEnter the no. of equations: ";        
    //cin >> n;
	
	n = 5;
    
    /* if no of equations are n then size of augmented matrix will be n*n+1. So here we are declaring 2d array 'mat' of size n+n+1 */
    long double mat[n][n+1];
	long double matOriginal[n][n+1];
    
    /* for n equations there will be n unknowns which will be stored in array 'res' */
    long double res[n];
   	
	/*
    cout << "\nEnter the elements of the augmented matrix: ";
    for(i=0;i<n;i++) {
     	for(j=0;j<n+1;j++) {
      		cin>>mat[i][j]; 
    	}    
    }
	*/
	
	// Coefficient Calculating (No Acceleration)
	
	/* Matrix Format
	{   0  ,  1  ,  2  ,  3  ,  4  ,   5   }
	{omega2, N_LF, N_RF, N_LR, N_RR, answer}
	
	// Equation 0, Sum of F_vertical = 0
	// (rho * C_Lift * A_Lift * R_COP^2 * cos(theta_COP)^2) / 2 * omega^2 + N_LF + N_RF + N_LR + N_RR = mg
	e0[0] = (rho * C_Lift * A_Lift * R_COP * R_COP * pow(cos(A_COP), 2)) / 2;
	e0[1] = 1;
	e0[2] = 1;
	e0[3] = 1;
	e0[4] = 1;
	e0[5] = m * g;
	
	// Equation 1, M at Front Wheels = 0
	// omega^2 * (rho * R_COP^2 * cos(theta_COP)^2 * (d_5 * C_Drag * A_Drag - (d_2 - d_4) * C_Lift * A_Lift) - d_1 * m * R_COM * sin(theta_COM)) / 2 - (d_2 + d_3) * (N_LR + N_RR) = -d_2 * mg
	e1[0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * ((d_5 * C_Drag * A_Drag) - ((d_2 - d_4) * C_Lift * A_Lift))) / 2) - d_1 * m * R_COM * sin(A_COM); 
	e1[1] = 0;
	e1[2] = 0;
	e1[3] = -(d_2 + d_3);
	e1[4] = -(d_2 + d_3);
	e1[5] = - d_2 * m * g;
	
	// Equation 2, M under COM at side = 0
	// omega^2 * (rho * R_COP^2 * cos(theta_COP)^2 * (d_5 * C_Drag * A_Drag + d_4 * C_Lift * A_Lift) - d_1 * m * R_COM * sin(theta_COM)) / 2 + d_2 * (N_LF + N_RF) - d_3 * (N_LR + N_RR) = 0
	e2[0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * (d_5 * C_Drag * A_Drag + d_4 * C_Lift * A_Lift)) / 2) - d_1 * m * R_COM * sin(A_COM);
	e2[1] = d_2;
	e2[2] = d_2;
	e2[3] = -d_3;
	e2[4] = -d_3;
	e2[5] = 0;
	
	// Equation 3, M under COM at front = 0
	e3[0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * d_10 * C_Lift * A_Lift) / 2) + d_1 * m * R_COM * cos(A_COM);
	e3[1] = -d_9;
	e3[2] = d_8;
	e3[3] = -d_7;
	e3[4] = d_6;
	e3[5] = 0;
	
	// Equation 4, F at front sideways = 0
	e3[0] = 0.25 * m * (R_LF * cos(A_LF) + R_RF * cos(A_RF));
	e3[1] = -mew * cos(A_LF);
	e3[2] = -mew * cos(A_RF);
	e3[3] = 0;
	e3[4] = 0;
	e3[5] = 0;
	*/
	
	// Equation 0, Sum of F_vertical = 0
	// (rho * C_Lift * A_Lift * R_COP^2 * cos(theta_COP)^2) / 2 * omega^2 + N_LF + N_RF + N_LR + N_RR = mg
	mat[0][0] = (rho * C_Lift * A_Lift * R_COP * R_COP * pow(cos(A_COP), 2)) / 2;
	mat[0][1] = 1;
	mat[0][2] = 1;
	mat[0][3] = 1;
	mat[0][4] = 1;
	mat[0][5] = m * g;
	
	/*
	// Equation 1, M at Front Wheels = 0
	// omega^2 * (rho * R_COP^2 * cos(theta_COP)^2 * (d_5 * C_Drag * A_Drag - (d_2 - d_4) * C_Lift * A_Lift) - d_1 * m * R_COM * sin(theta_COM)) / 2 - (d_2 + d_3) * (N_LR + N_RR) = -d_2 * mg
	mat[1][0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * ((d5 * C_Drag * A_Drag) - ((d2 - d4) * C_Lift * A_Lift))) / 2) - d1 * m * R_COM * sin(A_COM); 
	mat[1][1] = 0;
	mat[1][2] = 0;
	mat[1][3] = -(d2 + d3);
	mat[1][4] = -(d2 + d3);
	mat[1][5] = - d2 * m * g;
	*/
	
	/* Also Cancelling out
	// Equation 1, F at rear sideways = 0
	mat[1][0] = (m * (R_LR + R_RR)) / 4;
	mat[1][1] = 0;
	mat[1][2] = 0;
	mat[1][3] = -mew;
	mat[1][4] = -mew;
	mat[1][5] = 0;
	//*/
	
	// Equation 1, F in direction of Centrifugal Force = 0
	mat[1][0] = m * R_COM;
	mat[1][1] = -mew * cos(abs(A_COM - A_LF));
	mat[1][2] = -mew * cos(abs(A_COM - A_RF));
	mat[1][3] = -mew * cos(abs(A_COM));
	mat[1][4] = -mew * cos(abs(A_COM));
	mat[1][5] = 0;
	
	/* Also Cancelling out
	// Equation 1, F at front forwards = 0
	mat[1][0] = (m * (R_LF * sin(A_LF) + R_RF * sin(A_RF))) / 4 - ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * C_Drag * A_Drag) / 2);
	mat[1][1] = -mew * sin(A_LF);
	mat[1][2] = -mew * sin(A_RF);
	mat[1][3] = 0;
	mat[1][4] = 0;
	mat[1][5] = 0;
	//*/
	
	/*
	// Equation 2, M under COM at side = 0
	// omega^2 * (rho * R_COP^2 * cos(theta_COP)^2 * (d_5 * C_Drag * A_Drag + d_4 * C_Lift * A_Lift) - d_1 * m * R_COM * sin(theta_COM)) / 2 + d_2 * (N_LF + N_RF) - d_3 * (N_LR + N_RR) = 0
	mat[2][0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * (d5 * C_Drag * A_Drag + d4 * C_Lift * A_Lift)) / 2) - d1 * m * R_COM * sin(A_COM);
	mat[2][1] = d2;
	mat[2][2] = d2;
	mat[2][3] = -d3;
	mat[2][4] = -d3;
	mat[2][5] = 0;
	//*/
	
	// Equation 2, M around RF to LR = 0
	mat[2][0] = -((((rho *  R_COP * R_COP * pow(cos(A_COP), 2)) * (C_Lift * A_Lift * R_LR_COP * sin(A_RF_LR - A_LR_COP) - C_Drag * A_Drag * d5 * sin((PI / 2) - A_RF_LR))) / 2) - d1 * m * R_COM * sin((PI - A_COM) - A_RF_LR));
	mat[2][1] = R_LR_LF * sin(A_LR_LF - A_RF_LR);
	mat[2][2] = 0;
	mat[2][3] = 0;
	mat[2][4] = -(d7 + d6) * sin(A_RF_LR);
	mat[2][5] = - m * g * R_LR_COM * sin(A_RF_LR - A_LR_COM);
	
	/* Cancelling out
	// Equation 3, M under COM at front = 0
	mat[3][0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * d10 * C_Lift * A_Lift) / 2) + d1 * m * R_COM * cos(A_COM);
	mat[3][1] = -d9;
	mat[3][2] = d8;
	mat[3][3] = -d7;
	mat[3][4] = d6;
	mat[3][5] = 0;
	//*/
	
	///* Cancelling out
	// Equation 3, F at front sideways = 0
	mat[3][0] = (m * (R_LF * cos(A_LF) + R_RF * cos(A_RF))) / 4;
	mat[3][1] = -mew * cos(A_LF);
	mat[3][2] = -mew * cos(A_RF);
	mat[3][3] = 0;
	mat[3][4] = 0;
	mat[3][5] = 0;
	//*/
	
	/* Cancelling out
	// Equation 3, M RF at front = 0
	mat[3][0] = (((d8 - d10) * rho *  R_COP * R_COP * pow(cos(A_COP), 2) * C_Lift * A_Lift) / 2) - d1 * m * R_COM * cos(A_COM);
	mat[3][1] = d8 + d9;
	mat[3][2] = 0;
	mat[3][3] = d8 + d7;
	mat[3][4] = -(d6 - d8);
	mat[3][5] = d8 * m * g;
	*/
	
	/* Cancellig out
	// Equation 3, M LF at front = 0
	mat[3][0] = -(((d9 + d10) * rho *  R_COP * R_COP * pow(cos(A_COP), 2) * C_Lift * A_Lift) / 2) - d1 * m * R_COM * cos(A_COM);
	mat[3][1] = 0;
	mat[3][2] = -(d8 + d9);
	mat[3][3] = (d7 - d9);
	mat[3][4] = -(d6 + d9);
	mat[3][5] = -d9 * m * g;
	*/
	
	/*
	// Equation 4, F at front sideways = 0
	mat[4][0] = (m * (R_LF * cos(A_LF) + R_RF * cos(A_RF))) / 4;
	mat[4][1] = -mew * cos(A_LF);
	mat[4][2] = -mew * cos(A_RF);
	mat[4][3] = 0;
	mat[4][4] = 0;
	mat[4][5] = 0;
	//*/
	
	/*
	// Equation 4, F at front forwards = 0
	mat[4][0] = ((rho *  R_COP * R_COP * pow(cos(A_COP), 2) * C_Drag * A_Drag) / 2) + m * R_COM * sin(A_COM);
	mat[4][1] = -mew * sin(A_LF);
	mat[4][2] = -mew * sin(A_RF);
	mat[4][3] = 0;
	mat[4][4] = 0;
	mat[4][5] = 0;
	//*/
	
	/*
	// Equation 4, F sideways = 0
	mat[4][0] = m * R_COM * cos(A_COM);
	mat[4][1] = -mew * cos(A_LF);
	mat[4][2] = -mew * cos(A_RF);
	mat[4][3] = -mew;
	mat[4][4] = -mew;
	mat[4][5] = 0;
	//*/
	
	/*
	// Equation 4, F rear sideways = 0
	mat[4][0] = (m * (R_LR + R_RR)) / 4;
	mat[4][1] = 0;
	mat[4][2] = 0;
	mat[4][3] = -mew;
	mat[4][4] = -mew;
	mat[4][5] = 0;
	//*/
	
	/*
	// Equation 4, M LF at front = 0
	mat[4][0] = (((d8 - d10) * rho *  R_COP * R_COP * pow(cos(A_COP), 2) * C_Lift * A_Lift) / 2) - d1 * m * R_COM * cos(A_COM);
	mat[4][1] = d8 + d9;
	mat[4][2] = 0;
	mat[4][3] = d8 + d7;
	mat[4][4] = -(d6 - d8);
	mat[4][5] = d8 * m * g;
	//*/
	
	// Equation 4, M diagonal from LF to RR = 0
	mat[4][0] = (((rho *  R_COP * R_COP * pow(cos(A_COP), 2)) * (C_Lift * A_Lift * R_RR_COP * sin(A_LF_RR - A_RR_COP) - C_Drag * A_Drag * d5 * sin((PI / 2) - A_LF_RR))) / 2) - d1 * m * R_COM * sin(A_LF_RR - A_COM);
	mat[4][1] = 0;
	mat[4][2] = -R_RR_RF * sin(A_RR_RF - A_LF_RR);
	mat[4][3] = (d7 + d6) * sin(A_LF_RR);
	mat[4][4] = 0;
	mat[4][5] = R_RR_COM * m * g * sin(A_LF_RR - A_RR_COM);
	
	/*
	for (int coeff = 0; coeff < 6; coeff++) {
		cout << e0[coeff] << endl;
	} cout << endl;
	
	for (int coeff = 0; coeff < 6; coeff++) {
		cout << e1[coeff] << endl;
	} cout << endl;
	
	for (int coeff = 0; coeff < 6; coeff++) {
		cout << e2[coeff] << endl;
	} cout << endl;
	
	for (int coeff = 0; coeff < 6; coeff++) {
		cout << e3[coeff] << endl;
	} cout << endl;
	
	for (int coeff = 0; coeff < 6; coeff++) {
		cout << e4[coeff] << endl;
	} cout << endl;
	*/
	
	for (int ROW = 0; ROW < n; ROW++) {
		for (int COLUMN = 0; COLUMN < n + 1; COLUMN++) {
			cout << mat[ROW][COLUMN] << ", ";
			matOriginal[ROW][COLUMN] = mat[ROW][COLUMN];
		} cout << endl;
	} cout << endl;
  	
	for (int iter = 0; iter < 5; iter++) {
		for(i=0;i<n;i++) {                   
			for(j=i+1;j<n;j++) {
				if(abs(mat[i][i]) < abs(mat[j][i])) {
					for(k=0;k<n+1;k++) {
						// swapping mat[i][k] and mat[j][k]
						mat[i][k]=mat[i][k]+mat[j][k];
						mat[j][k]=mat[i][k]-mat[j][k];
						mat[i][k]=mat[i][k]-mat[j][k];
					}
				}
			}
		}
   
     	// performing Gaussian elimination
		for(i=0;i<n-1;i++) {
			for(j=i+1;j<n;j++) {
				float f=mat[j][i]/mat[i][i];
				for(k=0;k<n+1;k++) {
					mat[j][k]=mat[j][k]-f*mat[i][k];
				}
			}
		}
	}
	
	cout << "\n REF Form\n";
	for (int ROW = 0; ROW < n; ROW++) {
		for (int COLUMN = 0; COLUMN < n + 1; COLUMN++) {
			cout << fixed << setprecision(6) << mat[ROW][COLUMN] << ", ";
		} cout << endl;
	} cout << endl;
	
    /* Backward substitution for discovering values of unknowns */
    for(i=n-1;i>=0;i--) {                     
        res[i]=mat[i][n];
                    
        for(j=i+1;j<n;j++) {
         	if(i!=j) {
            	res[i]=res[i]-mat[i][j]*res[j];
    		}          
  		}
  		res[i]=res[i]/mat[i][i];  
    }
    
	
    cout<<"\nThe values of unknowns for the above equations=>\n";
    for(i=0;i<n;i++) {
    	cout<<res[i]<<"\n";
    }
	
	N_LF = res[1];
	N_RF = res[2];
	N_LR = res[3];
	N_RR = res[4];
	
	omega2 = res[0];
	omega = sqrt(omega2);
	linear_velocity = omega * R_COM;
	
	cout << "Omega = " << omega << " rad / s\n";
	cout << "Linear Velocity = " << linear_velocity << " m / s\n\n";
	cout << "radius of turning = " << d0 << "m\n";
	cout << "Speed = " << (linear_velocity * 1000) / 3600 << " km / hr\n\n";
	
	cout << "Answer Check\n";
	
	long double sum = 0;
	for (int eq = 0; eq < n; eq++) {
		sum = 0;
		for (int el = 0; el < n; el++) {
			sum += matOriginal[eq][el] * res[el];
		}
		sum -= matOriginal[eq][5];
		cout << "Equation " << eq << " = " << int(sum) << endl;
	}
	
	cout << "Sliding to the side? (Positive or 0 if not) = " << 0.5 * mew * ((abs(N_LF*cos(A_LF)) + N_LF*cos(A_LF)) + (abs(N_RF*cos(A_RF)) + N_RF*cos(A_RF)) + (abs(N_LR) + N_LR) + (abs(N_RR) + N_RR)) - abs(omega2 * m * R_COM * cos(A_COM)) << endl;
      
    return 0;
}