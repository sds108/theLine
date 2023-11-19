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

const int Unknowns = 4;
const int n = 4;

struct rref {
	
	long double mat[n][Unknowns+1];
	long double matOriginal[n][Unknowns+1];
	long double A_Temp[Unknowns+1];
	long double res[Unknowns];
	
	int iterations = 1;
	
	int i, j, k;
	
	long double divisor;
	
	void rearrange (int element) {
		int row = element;
		
		if (mat[row][element] == 0) {
			for (int temp_element = 0; temp_element < Unknowns + 1; temp_element++) {
				A_Temp[temp_element] = mat[row][temp_element];
			}
		} else return;
			
		for (row = element; row < n; row++) {
			if (mat[row][element] != 0) {
				for (int temp_element = 0; temp_element < Unknowns + 1; temp_element++) {
					mat[element][temp_element] = mat[row][temp_element];
					mat[row][temp_element] = A_Temp[temp_element];
				}
				return;
			}
		}
		
		cout << "Could not solve\n";
	}
				
	
	void solve () {
		/*
		cout << "\nEnter the elements of the augmented matrix: ";
		for(i=0;i<n;i++) {
			for(j=0;j<Unknowns+1;) {
				cin>>mat[i][j]; 
				matOriginal[i][j] = mat[i][j];
				j++;
			}    
		}
		*/
		
		mat[0][0] = 1;
		mat[0][1] = 1;
		mat[0][2] = 1;
		mat[0][3] = 1;
		mat[0][4] = 376329;
		
		mat[1][0] = 2.25773008412;
		mat[1][1] = 0;
		mat[1][2] = 0;
		mat[1][3] = -3.01030677883;
		mat[1][4] = -231513.851583;
		
		mat[2][0] = -1.5;
		mat[2][1] = 1.5;
		mat[2][2] = -2;
		mat[2][3] = 2;
		mat[2][4] = -45000;
		
		mat[3][0] = 0;
		mat[3][1] = 0;
		mat[3][2] = -4;
		mat[3][3] = -4;
		mat[3][4] = -945192;
		
		cout << fixed << setprecision(20);
		
		cout << "\n Original Form\n";
		for (int ROW = 0; ROW < n; ROW++) {
			for (int COLUMN = 0; COLUMN < n + 1; COLUMN++) {
				cout << mat[ROW][COLUMN] << ", ";
			} cout << endl;
		} cout << endl;
		
		// REF
		for (int iteration = 0; iteration < iterations; iteration++) {
			for (j = 0; j < Unknowns; j++) {
				rearrange(j);
				
				for (int row = j + 1; row < n; row++) {
					for (int element = 0; element < Unknowns + 1; element++) {
						A_Temp[element] = (mat[j][element] / mat[j][j]) * mat[row][j];
					}
					for (int element = 0; element < Unknowns + 1; element++) {
						mat[row][element] -= A_Temp[element];
					}
				}
			}
		}
		
		cout << "\n REF Form\n";
		for (int ROW = 0; ROW < n; ROW++) {
			for (int COLUMN = 0; COLUMN < n + 1; COLUMN++) {
				cout << mat[ROW][COLUMN] << ", ";
			} cout << endl;
		} cout << endl;
		
		// RREF
		for (int iteration = 0; iteration < iterations; iteration++) {
			for (j = Unknowns - 1; j > 0; j--) {
				
				for (int row = j - 1; row >= 0; row--) {
					for (int element = 0; element < Unknowns + 1; element++) {
						A_Temp[element] = (mat[j][element] / mat[j][j]) * mat[row][j];
					}
					for (int element = 0; element < Unknowns + 1; element++) {
						mat[row][element] -= A_Temp[element];
					}
				}
			}
		}
		
		cout << "\n RREF Form\n";
		for (int ROW = 0; ROW < n; ROW++) {
			for (int COLUMN = 0; COLUMN < n + 1; COLUMN++) {
				cout << mat[ROW][COLUMN] << ", ";
			} cout << endl;
		} cout << endl;
		
		long double multiplier = 0;
		
		// Solve
		for (int j = 0; j < Unknowns; j++) {
			res[j] = mat[j][Unknowns] / mat[j][j];
		}
		
		cout<<"\nThe values of unknowns for the above equations=>\n";
		for(i=0;i<n;i++) {
			cout<<res[i]<<"\n";
		}
		
		cout << "\nAnswer Check\n";
	
		long double sum = 0;
		for (int eq = 0; eq < n; eq++) {
			sum = 0;
			for (int el = 0; el < Unknowns; el++) {
				sum += matOriginal[eq][el] * res[el];
			}
			sum -= matOriginal[eq][n];
			
			cout << "Equation " << eq << " = " << int(sum) << endl;
		}
	}
};

int main () {
	rref RREF;
	
	RREF.solve();
	
	return 0;
}