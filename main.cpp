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

const double step = 0.01;
const double moveStep = 0.1;
const double interStep = 0.0001;
const double interStepPath = 0.01;
const double stepMultiplier = 10;
const double Vmax = 100; // m/s
const double k = 1;
const double nSlowDown = 0.0001;
const double ISNstep = 0.00001;
const int parts = 100;
vector<double> xinter;
vector<double> yinter;
vector<double> tinter;
const double mew = 1;
const double m = 1000;
const double g = 9.81;

const double curveMultiplier = 0;
const double lengthMultiplier = 1;
const double outsideMultiplier = 0;
const double distanceMultiplier = 0;
const double areaIntersectionMultiplier = 0;
		

double Bezier (double t, double a, double b, double c, double d) {
	return ((a * (-(t * t * t) + (3 * t * t) - (3 * t) + 1)) + (b * 3 * t * (1 - (2 * t) + (t * t))) + (c * 3 * t * t * (1 - t)) + (d * t * t * t));
	//return ((a*(1-t)*(1-t)*(1-t)) + (b*3*t*(1-t)*(1-t)) + (c*3*t*t*(1-t)) + (d*t*t*t)); // Slower
}

double dBezierdT (double t, double a, double b, double c, double d) {
	return ((Bezier(t + step, a, b, c, d) - Bezier(t - step, a, b, c, d)) / (2 * step));
}

double BezierTangentAngle (double t, double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	double DX = dBezierdT(t, ax, bx, cx, dx);
	double DY = dBezierdT(t, ay, by, cy, dy);
	
	if (DY == 0) {
		if (DX > 0) return 0;
		else if (DX < 0) return (PI);
		else return 0;
	} else if (DX == 0) {
		if (DY > 0) return (PI / 2);
		else if (DY < 0) return (-PI / 2);
		else return 0;
	} else return (atan2(DY, DX));
}

double BezierNormalAngle (double t, double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	double DX = dBezierdT(t, ax, bx, cx, dx);
	double DY = dBezierdT(t, ay, by, cy, dy);
	
	if (DY == 0) {
		if (DX > 0) return (-PI / 2);
		else if (DX < 0) return (PI / 2);
		else return 0;
	} else if (DX == 0) {
		if (DY > 0) return (-PI);
		else if (DY < 0) return (PI);
		else return 0;
	} else return (atan2(-DX, DY));
}

double xTransform (double t, double pointx, double pointy, double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	double normalAngle = BezierNormalAngle(t, ax, ay, bx, by, cx, cy, dx, dy);
	return (((pointx - Bezier(t, ax, bx, cx, dx)) * (cos(-normalAngle))) - ((pointy - Bezier(t, ay, by, cy, dy)) * (sin(-normalAngle))));
}

double distanceSimple (double t, double pointx, double pointy, double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	return (sqrt(pow((pointx - Bezier(t, ax, bx, cx, dx)), 2) + pow((pointy - Bezier(t, ay, by, cy, dy)), 2)));
}

double distanceSigned (double t, double pointx, double pointy, double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	double xnew = xTransform(t, pointx, pointy, ax, ay, bx, by, cx, cy, dx, dy);
	
	if (xnew > 0) {
		return (sqrt(pow((pointx - Bezier(t, ax, bx, cx, dx)), 2) + pow((pointy - Bezier(t, ay, by, cy, dy)), 2)));
	} else return (-sqrt(pow((pointx - Bezier(t, ax, bx, cx, dx)), 2) + pow((pointy - Bezier(t, ay, by, cy, dy)), 2)));
}

bool inside (double t, double pointx, double pointy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy, double lax, double lay, double lbx, double lby, double lcx, double lcy, double ldx, double ldy) {
	if (xTransform(t, pointx, pointy, rax, ray, rbx, rby, rcx, rcy, rdx, rdy) < 0 && xTransform(t, pointx, pointy, lax, lay, lbx, lby, lcx, lcy, ldx, ldy) > 0) return (1 - abs(((sqrt(pow((pointx - Bezier(t + step, rax, rbx, rcx, rdx)), 2) + pow((pointy - Bezier(t + step, ray, rby, rcy, rdy)), 2)) + sqrt(pow((pointx - Bezier(t + step, lax, lbx, lcx, ldx)), 2) + pow((pointy - Bezier(t + step, lay, lby, lcy, ldy)), 2))) - (sqrt(pow((pointx - Bezier(t - step, rax, rbx, rcx, rdx)), 2) + pow((pointy - Bezier(t - step, ray, rby, rcy, rdy)), 2)) + sqrt(pow((pointx - Bezier(t - step, lax, lbx, lcx, ldx)), 2) + pow((pointy - Bezier(t - step, lay, lby, lcy, ldy)), 2)))) / (2 * step))) >= 0.9;
	else return false;
}

bool insideSegment (double cax, double cay, double cbx, double cby, double ccx, double ccy, double cdx, double cdy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy, double lax, double lay, double lbx, double lby, double lcx, double lcy, double ldx, double ldy) {
	double BX = 0;
	double BY = 0;
	int sum = 0;
	int o = 0;
	for (int i = 0; i < 1 / (step*stepMultiplier) + 1; i++) {
		BX = Bezier(i*(step*stepMultiplier), cax, cbx, ccx, cdx);
		BY = Bezier(i*(step*stepMultiplier), cay, cby, ccy, cdy);
		sum = 0;
		for (int n = 0; n < 1 / step + 1; n++) {
			sum += int(inside(n*step, BX, BY, rax, ray, rbx, rby, rcx, rcy, rdx, rdy, lax, lay, lbx, lby, lcx, lcy, ldx, ldy));
		}
		if (sum == 0) {
			cout << i*(step*stepMultiplier) << endl << endl;
			return false;
		}
	}
	return true;
}

double Lambda (double Bx, double By, double x0, double x1, double y0, double y1) {
	return (Bx - (By * ((x1 - x0) / (y1 - y0))) - (x0 - ((y0 * (x1 - x0)) / (y1 - y0))));
}

double magicEquation (double t, double H, double pointx, double pointy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy, double lax, double lay, double lbx, double lby, double lcx, double lcy, double ldx, double ldy) {
	return sqrt(pow(Bezier(t + H, rax, rbx, rcx, rdx) - pointx, 2) + pow(Bezier(t + H, ray, rby, rcy, rdy) - pointy, 2)) + sqrt(pow(Bezier(t + H, lax, lbx, lcx, ldx) - pointx, 2) + pow(Bezier(t + H, lay, lby, lcy, ldy) - pointy, 2));
}

// return false if intersecting
bool lineSegmentTest (double cax, double cay, double cbx, double cby, double ccx, double ccy, double cdx, double cdy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy, double lax, double lay, double lbx, double lby, double lcx, double lcy, double ldx, double ldy) {
	double x0 = 0;
	double x1 = 0;
	
	double y0 = 0;
	double y1 = 0;
	
	double lambda = 0;
	
	const int iterations = 100;
	
	double previousGuess = 0;
	double nextGuess = 0;
	
	double firstGuessRight = 0;
	double secondGuessRight = 0;
	double thirdGuessRight = 0;
	
	double firstGuessLeft = 0;
	double secondGuessLeft = 0;
	double thirdGuessLeft = 0;
	
	for (int i = 1; i < 1 / (step); i++) {
		x0 = Bezier(i * (step), cax, cbx, ccx, cdx);
		x1 = Bezier((i + 1) * (step), cax, cbx, ccx, cdx);
		
		y0 = Bezier(i * (step), cay, cby, ccy, cdy);
		y1 = Bezier((i + 1) * (step), cay, cby, ccy, cdy);
		
		// Newtons from 0.0001 for right side
		nextGuess = 0.0001;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store first guess
		firstGuessRight = nextGuess;
		
		// Newtons from 0.5 for right side
		nextGuess = 0.5;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store second guess
		secondGuessRight = nextGuess;
		
		// Newtons from 0.9999 for right side
		nextGuess = 0.9999;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store third guess
		thirdGuessRight = nextGuess;
		
		
		
		// Newtons from 0.0001 for left side
		nextGuess = 0.0001;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, lax, lbx, lcx, ldx), Bezier(previousGuess, lay, lby, lcy, ldy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, lax, lbx, lcx, ldx), Bezier(previousGuess + step, lay, lby, lcy, ldy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, lax, lbx, lcx, ldx), Bezier(previousGuess - step, lay, lby, lcy, ldy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, lax, lbx, lcx, ldx) - (Bezier(previousGuess, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, lax, lbx, lcx, ldx) - (Bezier(previousGuess + step, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, lax, lbx, lcx, ldx) - (Bezier(previousGuess - step, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store first guess
		firstGuessLeft = nextGuess;
		
		// Newtons from 0.5 for left side
		nextGuess = 0.5;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, lax, lbx, lcx, ldx), Bezier(previousGuess, lay, lby, lcy, ldy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, lax, lbx, lcx, ldx), Bezier(previousGuess + step, lay, lby, lcy, ldy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, lax, lbx, lcx, ldx), Bezier(previousGuess - step, lay, lby, lcy, ldy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, lax, lbx, lcx, ldx) - (Bezier(previousGuess, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, lax, lbx, lcx, ldx) - (Bezier(previousGuess + step, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, lax, lbx, lcx, ldx) - (Bezier(previousGuess - step, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store second guess
		secondGuessLeft = nextGuess;
		
		// Newtons from 0.9999 for left side
		nextGuess = 0.9999;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, lax, lbx, lcx, ldx), Bezier(previousGuess, lay, lby, lcy, ldy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, lax, lbx, lcx, ldx), Bezier(previousGuess + step, lay, lby, lcy, ldy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, lax, lbx, lcx, ldx), Bezier(previousGuess - step, lay, lby, lcy, ldy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, lax, lbx, lcx, ldx) - (Bezier(previousGuess, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, lax, lbx, lcx, ldx) - (Bezier(previousGuess + step, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, lax, lbx, lcx, ldx) - (Bezier(previousGuess - step, lay, lby, lcy, ldy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store third guess
		thirdGuessLeft = nextGuess;
		
		
		// Test different Lambdas
		lambda = (Bezier(firstGuessRight, ray, rby, rcy, rdy) - y0) / (y1 - y0);
		//cout << firstGuessRight << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (Bezier(secondGuessRight, ray, rby, rcy, rdy) - y0) / (y1 - y0);
		//cout << secondGuessRight << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (Bezier(thirdGuessRight, ray, rby, rcy, rdy) - y0) / (y1 - y0);
		//cout << thirdGuessRight << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (Bezier(firstGuessLeft, lay, lby, lcy, ldy) - y0) / (y1 - y0);
		//cout << firstGuessLeft << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (Bezier(secondGuessLeft, lay, lby, lcy, ldy) - y0) / (y1 - y0);
		//cout << secondGuessLeft << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (Bezier(thirdGuessLeft, lay, lby, lcy, ldy) - y0) / (y1 - y0);
		//cout << thirdGuessLeft << endl;
		//cout << lambda << endl << endl;
		if (lambda > 0 && lambda < 1) return false;
	}
	
	return true;
}

double NormalAnglePixels (double x0, double y0, double x1, double y1) {
	double DX = x1 - x0;
	double DY = y1 - y0;
	
	if (DY == 0) {
		if (DX > 0) return (-PI / 2);
		else if (DX < 0) return (PI / 2);
		else return 0;
	} else if (DX == 0) {
		if (DY > 0) return (-PI);
		else if (DY < 0) return (PI);
		else return 0;
	} else return (atan2(-DX, DY));
}

bool xTransformPixelsOutside (bool leftSide, double pathX, double pathY, double borderX0, double borderY0, double borderX1, double borderY1) {
	double normalAngle = NormalAnglePixels(borderX0, borderY0, borderX1, borderY1);
	return ((((pathX - borderX0) * (cos(-normalAngle))) - ((pathY - borderY0) * (sin(-normalAngle)))) < 0 && leftSide) || ((((pathX - borderX0) * (cos(-normalAngle))) - ((pathY - borderY0) * (sin(-normalAngle)))) > 0 && !leftSide);
}

bool xTransformBezierOutside (bool leftSide, double t, double pointx, double pointy, double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy) {
	double normalAngle = BezierNormalAngle(t, ax, ay, bx, by, cx, cy, dx, dy);
	return (leftSide && (((pointx - Bezier(t, ax, bx, cx, dx)) * (cos(-normalAngle))) - ((pointy - Bezier(t, ay, by, cy, dy)) * (sin(-normalAngle)))) < 0) ||  (!leftSide && (((pointx - Bezier(t, ax, bx, cx, dx)) * (cos(-normalAngle))) - ((pointy - Bezier(t, ay, by, cy, dy)) * (sin(-normalAngle)))) > 0);
}

bool exitingBorder (bool leftSide, double pathX0, double pathY0, double pathX1, double pathY1, double borderX0, double borderY0, double borderX1, double borderY1) {
	double normalAngle = NormalAnglePixels(borderX0, borderY0, borderX1, borderY1);
	
	double pathX0Transformed = ((pathX0 - borderX0) * (cos(-normalAngle))) - ((pathY0 - borderY0) * (sin(-normalAngle)));
	double pathY0Transformed = ((pathX0 - borderX0) * (sin(-normalAngle))) + ((pathY0 - borderY0) * (cos(-normalAngle)));
	
	double pathX1Transformed = ((pathX1 - borderX0) * (cos(-normalAngle))) - ((pathY1 - borderY0) * (sin(-normalAngle)));
	double pathY1Transformed = ((pathX1 - borderX0) * (sin(-normalAngle))) + ((pathY1 - borderY0) * (cos(-normalAngle)));
	
	double DX = pathX1Transformed - pathX0Transformed;
	double DY = pathY1Transformed - pathY0Transformed;
	
	long double pathAngle = abs(atan2(DX, DY)) * 18000 / PI;
	
	cout << "Path Angle : " << pathAngle << endl;
	
	return (leftSide && pathAngle > 9000) || (!leftSide && pathAngle < 9000);
}

bool exitingBorderBezier (bool leftSide, double pathTime, double borderTime, double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy, double bax, double bay, double bbx, double bby, double bcx, double bcy, double bdx, double bdy) {
	
	double pathX0 = Bezier(pathTime, pax, pbx, pcx, pdx);
	double pathY0 = Bezier(pathTime, pay, pby, pcy, pdy);
	
	double pathX1 = Bezier(pathTime + step, pax, pbx, pcx, pdx);
	double pathY1 = Bezier(pathTime + step, pay, pby, pcy, pdy);
	
	double borderX0 = Bezier(borderTime, bax, bbx, bcx, bdx);
	double borderY0 = Bezier(borderTime, bay, bby, bcy, bdy);
	
	double normalAngle = BezierTangentAngle(borderTime, bax, bay, bbx, bby, bcx, bcy, bdx, bdy);
	
	double pathX0Transformed = ((pathX0 - borderX0) * (cos(-normalAngle))) - ((pathY0 - borderY0) * (sin(-normalAngle)));
	double pathY0Transformed = ((pathX0 - borderX0) * (sin(-normalAngle))) + ((pathY0 - borderY0) * (cos(-normalAngle)));
	
	double pathX1Transformed = ((pathX1 - borderX0) * (cos(-normalAngle))) - ((pathY1 - borderY0) * (sin(-normalAngle)));
	double pathY1Transformed = ((pathX1 - borderX0) * (sin(-normalAngle))) + ((pathY1 - borderY0) * (cos(-normalAngle)));
	
	double DX = pathX1Transformed - pathX0Transformed;
	double DY = pathY1Transformed - pathY0Transformed;
	
	long double pathAngle = abs(atan2(DX, DY)) * 1800000 / PI;
	
	return (leftSide && pathAngle < 900000) || (!leftSide && pathAngle > 900000);
}

struct line {
	vector <int> xCoord;
	vector <int> yCoord;
	
	line (int x0, int y0, int x1, int y1) {
		editLine(x0, y0, x1, y1);
	}
	
	void editLine (int x0, int y0, int x1, int y1) {
		int rise = y1 - y0;
		int run = x1 - x0;
		
		if (run == 0) {
			if (y1 < y0) {
				for (int pixel = y1; pixel <= y0; pixel++) {
					xCoord.push_back(x0);
					yCoord.push_back(pixel);
				}
			} else {
				for (int pixel = y0; pixel <= y1; pixel++) {
					xCoord.push_back(x0);
					yCoord.push_back(pixel);
				}
			}
		} else {
			long double slope = double(rise) / double(run);
			long double delta = abs(slope);
			int adjust = 0;
			long double offset = 0;
			double threshold = 0.5;
			
			if (slope >= 0) adjust = 1;
			else adjust = -1;
			
			if (slope <= 1 && slope >= -1) {
				int y = y0;
				
				if (x1 < x0) {
					y = y1;
					for (int pixel = x1; pixel <= x0; pixel++) {
						xCoord.push_back(pixel);
						yCoord.push_back(y);
						
						offset += delta;
						
						if (offset >= threshold) {
							y += adjust;
							threshold++;
						}
					}
				} else {
					y = y0;
					for (int pixel = x0; pixel <= x1; pixel++) {
						xCoord.push_back(pixel);
						yCoord.push_back(y);
						
						offset += delta;
						
						if (offset >= threshold) {
							y += adjust;
							threshold++;
						}
					}
				}
			} else {
				delta = abs(double(run) / double(rise));
				
				int x = x0;
				
				if (y1 < y0) {
					x = x1;
					for (int pixel = y1; pixel <= y0; pixel++) {
						xCoord.push_back(x);
						yCoord.push_back(pixel);
						
						if (offset >= threshold) {
							x += adjust;
							threshold++;
						}
					}
				} else {
					x = x0;
					
					for (int pixel = y1; pixel <= y0; pixel++) {
						xCoord.push_back(x);
						yCoord.push_back(pixel);
						
						if (offset >= threshold) {
							x += adjust;
							threshold++;
						}
					}
				}
			}
		}
	}
};

struct BezierSegment {
	
	BezierSegment () {}
	
	vector <int> xCoord;
	vector <int> yCoord;
	vector <long double> tCoord;
	
	void produceLines () {
		xCoord.clear();
		yCoord.clear();
		tCoord.clear();
		
		for (int i = 0; i < (1 / interStep); i++) {
			int x0 = Bezier(i * interStep, params[0], params[2], params[4], params[6]);
			int x1 = Bezier((i + 1) * interStep, params[0], params[2], params[4], params[6]);

			int y0 = Bezier(i * interStep, params[1], params[3], params[5], params[7]);
			int y1 = Bezier((i + 1) * interStep, params[1], params[3], params[5], params[7]);

			editLine(i * interStep, x0, y0, x1, y1);
		}
	}
	
	void editLine (double t, int x0, int y0, int x1, int y1) {
		int rise = y1 - y0;
		int run = x1 - x0;
		
		if (run == 0) {
			if (y1 < y0) {
				for (int pixel = y0; pixel > y1; pixel--) {
					xCoord.push_back(x0);
					yCoord.push_back(pixel);
					tCoord.push_back(t + ((interStep) * abs(double(pixel - y0) / double(rise))));
				}
			} else {
				for (int pixel = y0; pixel < y1; pixel++) {
					xCoord.push_back(x0);
					yCoord.push_back(pixel);
					tCoord.push_back(t + ((interStep) * abs(double(pixel - y0) / double(rise))));
				}
			}
		} else {
			double slope = double(rise) / double(run);
			double delta = abs(slope);
			int adjust = 0;
			long double offset = 0;
			double threshold = 0.5;
			
			if (slope >= 0) adjust = 1;
			else adjust = -1;
			
			if (slope <= 1 && slope >= -1) {
				int y = y0;
				
				if (x1 < x0) {
					y = y0;
					for (int pixel = x0; pixel > x1; pixel--) {
						xCoord.push_back(pixel);
						yCoord.push_back(y);
						tCoord.push_back(t + ((interStep) * (double(pixel - x1) / abs(run))));
						
						offset += delta;
						
						if (offset >= threshold) {
							y -= adjust;
							threshold++;
						}
					}
				} else {
					y = y0;
					for (int pixel = x0; pixel < x1; pixel++) {
						xCoord.push_back(pixel);
						yCoord.push_back(y);
						tCoord.push_back(t + ((interStep) * (double(pixel - x0) / abs(run))));
						
						offset += delta;
						
						if (offset >= threshold) {
							y += adjust;
							threshold++;
						}
					}
				}
			} else {
				delta = abs(double(run) / double(rise));
				
				int x = x0;
				
				if (y1 < y0) {
					x = x0;
					
					for (int pixel = y0; pixel > y1; pixel--) {
						xCoord.push_back(x);
						yCoord.push_back(pixel);
						tCoord.push_back(t + ((interStep) * (double(pixel - y1) / abs(rise))));
						
						offset += delta;
						
						if (offset >= threshold) {
							x -= adjust;
							threshold++;
						}
					}
				} else {
					x = x0;
					
					for (int pixel = y0; pixel < y1; pixel++) {
						xCoord.push_back(x);
						yCoord.push_back(pixel);
						tCoord.push_back(t + ((interStep) * (double(pixel - y0) / abs(rise))));
						
						offset += delta;
						
						if (offset >= threshold) {
							x += adjust;
							threshold++;
						}
					}
				}
			}
		}
	}
	
	double ax = 0;
	double ay = 0;
	double controlLength = 1;
	double bx = 0;
	double by = 0;
	double cx = 0;
	double cy = 0;
	double dx = 0;
	double dy = 0;
	double curvature = 0;
	double mx = 1;
	double my = 1;
	
	double params[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // ax, ay, bx, by, cx, cy, dx, dy, CL
	
	void setupFirst (double AX, double AY, double BX, double BY, double CX, double CY, double DX, double DY) {
		ax = AX;
		ay = AY;
		bx = BX;
		by = BY;
		cx = CX;
		cy = CY;
		dx = DX;
		dy = DY;
		controlLength = 1;
		
		params[0] = AX;
		params[1] = AY;
		params[2] = BX;
		params[3] = BY;
		params[4] = CX;
		params[5] = CY;
		params[6] = DX;
		params[7] = DY;
		params[8] = 1;
	}
	
	void setupSecondary (double AX, double AY, double BX, double BY, double CX, double CY, double DX, double DY, double CL) {
		ax = AX;
		ay = AY;
		bx = ax + ((ax - BX) * controlLength);
		by = ay + ((ay - BY) * controlLength);
		cx = CX;
		cy = CY;
		dx = DX;
		dy = DY;
		controlLength = CL;
		
		params[0] = AX;
		params[1] = AY;
		params[8] = CL;
		params[2] = params[0] + ((params[0] - BX) * params[8]);
		params[3] = params[1] + ((params[1] - BY) * params[8]);
		params[4] = CX;
		params[5] = CY;
		params[6] = DX;
		params[7] = DY;
	}
	
	void updateBX (double BX) {
		mx = params[0] - BX;
		params[2] = params[0] + (mx * params[8]);
	}
	
	void updateBY (double BY) {
		my = params[1] - BY;
		params[3] = params[1] + (my * params[8]);
	}
	
	void updateCL (double CL) {
		params[8] = CL;
		params[2] = params[0] + (mx * params[8]);
		params[3] = params[1] + (my * params[8]);
	}
	
	double returnBezierX (double t) {
		return Bezier (t, params[0], params[2], params[4], params[6]);
	}
	
	double returnBezierY (double t) {
		return Bezier (t, params[1], params[3], params[5], params[7]);
	}
	
	double BY0 = 0;
	double BY1 = 0;
	double BY2 = 0;
	double BX0 = 0;
	double BX1 = 0;
	double BX2 = 0;
	double angle0 = 0;
	double angle1 = 0;
	double length = 0;
	
	double Curvature () {
		curvature = 0;
		for (int i = 1; i < 1 / step; i++) {
			BX0 = returnBezierX(i * step - step);
			BX1 = returnBezierX(i * step);
			BX2 = returnBezierX(i * step + step);
			
			BY0 = returnBezierY(i * step - step);
			BY1 = returnBezierY(i * step);
			BY2 = returnBezierY(i * step + step);
			
			if (BX2 - BX1 == 0) angle1 = 90;
			if (BX1 - BX0 == 0) angle0 = 90;
			
			length = sqrt(pow(BX2 - BX0, 2) + pow(BY2 - BY0, 2));
			
			if (length > 0) curvature += abs((atan(abs((BY2 - BY1) / (BX2 - BX1))) - atan(abs((BY1 - BY0) / (BX1 - BX0)))) / length);
		}
		
		return curvature * 2 * step;
	}
	
	double CurvatureV2 (double t) {
		curvature = 0;
		BX0 = returnBezierX(t - step);
		BX1 = returnBezierX(t);
		BX2 = returnBezierX(t + step);

		BY0 = returnBezierY(t - step);
		BY1 = returnBezierY(t);
		BY2 = returnBezierY(t + step);

		if (BX2 - BX1 == 0) angle1 = 90;
		if (BX1 - BX0 == 0) angle0 = 90;

		length = sqrt(pow(BX2 - BX0, 2) + pow(BY2 - BY0, 2));

		if (length > 0) curvature = abs((atan(abs((BY2 - BY1) / (BX2 - BX1))) - atan(abs((BY1 - BY0) / (BX1 - BX0)))) / length);
		
		return curvature * 2 * step;
	}
	
	double sDiff (double t) {
		BX0 = returnBezierX(t - step);
		BX1 = returnBezierX(t);
		BX2 = returnBezierX(t + step);

		BY0 = returnBezierY(t - step);
		BY1 = returnBezierY(t);
		BY2 = returnBezierY(t + step);

		if (BX2 - BX1 == 0) angle1 = 90;
		if (BX1 - BX0 == 0) angle0 = 90;

		length = sqrt(pow(BX2 - BX0, 2) + pow(BY2 - BY0, 2));
		
		return length;
	}
};

struct BezierPiecewise {
	int type = 0; // 0 = border, 1 = path.
	
	bool leftSide = false;
	
	vector <BezierSegment*> segments;
	
	double curvature = 0;
	
	void firstSegment (double AX, double AY, double BX, double BY, double CX, double CY, double DX, double DY) {
		segments.push_back(new BezierSegment());
		segments[0]->setupFirst(AX, AY, BX, BY, CX, CY, DX, DY);
	}
	
	void addSegment (double CX, double CY, double DX, double DY, double CL) {
		segments.push_back(new BezierSegment());
		segments[segments.size() - 1]->setupSecondary(segments[segments.size() - 2]->dx, segments[segments.size() - 2]->dy, segments[segments.size() - 2]->cx, segments[segments.size() - 2]->cy, CX, CY, DX, DY, CL);
	}
	
	double totalCurvature () {
		curvature = 0;
		for (int i = 0; i < segments.size(); i++) {
			curvature += segments[i]->Curvature();
		}
		return curvature;
	}
	
	void produceLines () {
		for (int seg = 0; seg < segments.size(); seg++) {
			segments[seg]->produceLines();
		}
	}
};

// return false if intersecting
bool lineSegmentTestSingle (BezierSegment* path, BezierSegment* border) {
	double x0 = 0;
	double x1 = 0;
	
	double y0 = 0;
	double y1 = 0;
	
	double lambda = 0;
	
	const int iterations = 100;
	
	double previousGuess = 0;
	double nextGuess = 0;
	
	double firstGuess = 0;
	double secondGuess = 0;
	double thirdGuess = 0;
	
	for (int i = 1; i < 1 / (step); i++) {
		x0 = path->returnBezierX(i * (step));
		x1 = path->returnBezierX((i + 1) * (step));
		
		y0 = path->returnBezierY(i * (step));
		y1 = path->returnBezierY((i + 1) * (step));
		
		// Newtons from 0.0001
		nextGuess = 0.0001;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(border->returnBezierX(previousGuess), border->returnBezierY(previousGuess), x0, x1, y0, y1)) / (Lambda(border->returnBezierX(previousGuess + step), border->returnBezierY(previousGuess + step), x0, x1, y0, y1) - Lambda(border->returnBezierX(previousGuess + step), border->returnBezierY(previousGuess + step), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store first guess
		firstGuess = nextGuess;
		
		// Newtons from 0.5
		nextGuess = 0.5;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(border->returnBezierX(previousGuess), border->returnBezierY(previousGuess), x0, x1, y0, y1)) / (Lambda(border->returnBezierX(previousGuess + step), border->returnBezierY(previousGuess + step), x0, x1, y0, y1) - Lambda(border->returnBezierX(previousGuess + step), border->returnBezierY(previousGuess + step), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store second guess
		secondGuess = nextGuess;
		
		// Newtons from 0.9999
		nextGuess = 0.9999;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(border->returnBezierX(previousGuess), border->returnBezierY(previousGuess), x0, x1, y0, y1)) / (Lambda(border->returnBezierX(previousGuess + step), border->returnBezierY(previousGuess + step), x0, x1, y0, y1) - Lambda(border->returnBezierX(previousGuess + step), border->returnBezierY(previousGuess + step), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
		}
		
		// Store third guess
		thirdGuess = nextGuess;
		
		
		
		// Test different Lambdas
		lambda = (border->returnBezierY(firstGuess) - y0) / (y1 - y0);
		//cout << firstGuessRight << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (border->returnBezierY(secondGuess) - y0) / (y1 - y0);
		//cout << secondGuessRight << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
		
		lambda = (border->returnBezierY(thirdGuess) - y0) / (y1 - y0);
		//cout << thirdGuessRight << endl;
		//cout << lambda << endl;
		if (lambda > 0 && lambda < 1) return false;
	}
	
	return true;
}

// return true if intersecting
bool notOverlapping (double cax, double cay, double cbx, double cby, double ccx, double ccy, double cdx, double cdy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy) {
	double x0 = 0;
	double x1 = 0;
	
	double y0 = 0;
	double y1 = 0;
	
	double lambda = 0;
	
	const int iterations = 100;
	const double acceptableError = 0.05;
	const double minimumDifference = 0.001;
	
	double previousGuess = 0;
	double nextGuess = 0;
	
	double firstGuess = 0;
	double secondGuess = 0;
	double thirdGuess = 0;
	
	for (int i = 1; i < 1 / (step); i++) {
		x0 = Bezier(i * (step), cax, cbx, ccx, cdx);
		x1 = Bezier((i + 1) * (step), cax, cbx, ccx, cdx);
		
		y0 = Bezier(i * (step), cay, cby, ccy, cdy);
		y1 = Bezier((i + 1) * (step), cay, cby, ccy, cdy);
		
		// Newtons from 0.0001
		nextGuess = 0.0001;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store first guess
		firstGuess = nextGuess;
		
		// Newtons from 0.5
		nextGuess = 0.5;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store second guess
		secondGuess = nextGuess;
		
		// Newtons from 0.9999
		nextGuess = 0.9999;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store third guess
		thirdGuess = nextGuess;	
		
		
		// Test different Lambdas		
		if (firstGuess > 0 && firstGuess < 1) {
			lambda = (Bezier(firstGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			//cout << firstGuessRight << endl;
			//cout << lambda << endl;
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(firstGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(firstGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
					cout << firstGuess << "t, ";
					return false;
				}
			}
		}
		
		if (secondGuess > 0 && secondGuess < 1) {
			lambda = (Bezier(secondGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			//cout << secondGuessRight << endl;
			//cout << lambda << endl;
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(secondGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(secondGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
					cout << secondGuess << "t, ";
					return false;
				}
			}
		}
		
		if (thirdGuess > 0 && thirdGuess < 1) {
			lambda = (Bezier(thirdGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			//cout << thirdGuessRight << endl;
			//cout << lambda << endl;
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(thirdGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(thirdGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
					cout << thirdGuess << "t, ";
					return false;
				}
			}
		}
	}
	
	return true;
}
	
// return true if intersecting
bool Overlapping (double cax, double cay, double cbx, double cby, double ccx, double ccy, double cdx, double cdy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy) {
	double x0 = 0;
	double x1 = 0;
	
	double y0 = 0;
	double y1 = 0;
	
	double lambda = 0;
	
	const int iterations = 100;
	const double acceptableError = 0.05;
	const double minimumDifference = 0.001;
	
	double previousGuess = 0;
	double nextGuess = 0;
	
	double firstGuess = 0;
	double secondGuess = 0;
	double thirdGuess = 0;
	
	for (int i = 1; i < 1 / (step); i++) {
		x0 = Bezier(i * (step), cax, cbx, ccx, cdx);
		x1 = Bezier((i + 1) * (step), cax, cbx, ccx, cdx);
		
		y0 = Bezier(i * (step), cay, cby, ccy, cdy);
		y1 = Bezier((i + 1) * (step), cay, cby, ccy, cdy);
		
		// Newtons from 0.0001
		nextGuess = 0.0001;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store first guess
		firstGuess = nextGuess;
		
		// Newtons from 0.5
		nextGuess = 0.5;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store second guess
		secondGuess = nextGuess;
		
		// Newtons from 0.9999
		nextGuess = 0.9999;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store third guess
		thirdGuess = nextGuess;	
		
		
		// Test different Lambdas		
		if (firstGuess > 0 && firstGuess < 1) {
			lambda = (Bezier(firstGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			//cout << firstGuessRight << endl;
			//cout << lambda << endl;
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(firstGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(firstGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
					cout << firstGuess << "t, ";
					return true;
				}
			}
		}
		
		if (secondGuess > 0 && secondGuess < 1) {
			lambda = (Bezier(secondGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			//cout << secondGuessRight << endl;
			//cout << lambda << endl;
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(secondGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(secondGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
					cout << secondGuess << "t, ";
					return true;
				}
			}
		}
		
		if (thirdGuess > 0 && thirdGuess < 1) {
			lambda = (Bezier(thirdGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			//cout << thirdGuessRight << endl;
			//cout << lambda << endl;
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(thirdGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(thirdGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
					cout << thirdGuess << "t, ";
					return true;
				}
			}
		}
	}
	
	return false;
}

// return true if intersecting
bool OverlappingV2 (double cax, double cay, double cbx, double cby, double ccx, double ccy, double cdx, double cdy, double rax, double ray, double rbx, double rby, double rcx, double rcy, double rdx, double rdy) {
	double x0 = 0;
	double x1 = 0;
	
	double y0 = 0;
	double y1 = 0;
	
	double lambda = 0;
	
	const int iterations = 100;
	const double acceptableError = 0.1;
	const double minimumDifference = 0.001;
	
	double previousGuess = 0;
	double nextGuess = 0;
	
	double firstGuess = 0;
	double secondGuess = 0;
	double thirdGuess = 0;
	
	for (int i = 1; i < 1 / (step); i++) {
		x0 = Bezier(i * (step), cax, cbx, ccx, cdx);
		x1 = Bezier((i + 1) * (step), cax, cbx, ccx, cdx);
		
		y0 = Bezier(i * (step), cay, cby, ccy, cdy);
		y1 = Bezier((i + 1) * (step), cay, cby, ccy, cdy);
		
		// Newtons from 0.0001
		nextGuess = 0.0001;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store first guess
		firstGuess = nextGuess;
		
		// Newtons from 0.5
		nextGuess = 0.5;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store second guess
		secondGuess = nextGuess;
		
		// Newtons from 0.9999
		nextGuess = 0.9999;
		for (int n = 0; n < iterations; n++) {
			previousGuess = nextGuess;
			nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, rax, rbx, rcx, rdx), Bezier(previousGuess, ray, rby, rcy, rdy), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, rax, rbx, rcx, rdx), Bezier(previousGuess + step, ray, rby, rcy, rdy), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, rax, rbx, rcx, rdx), Bezier(previousGuess - step, ray, rby, rcy, rdy), x0, x1, y0, y1)));
			//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
			if (abs(previousGuess - nextGuess) < minimumDifference) break;
		}
		
		// Store third guess
		thirdGuess = nextGuess;	
		
		
		// Test different Lambdas		
		if (firstGuess > 0 && firstGuess < 1) {
			lambda = (Bezier(firstGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(firstGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(firstGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) return true;
			}
		}
		
		if (secondGuess > 0 && secondGuess < 1) {
			lambda = (Bezier(secondGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(secondGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(secondGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) return true;
			}
		}
		
		if (thirdGuess > 0 && thirdGuess < 1) {
			lambda = (Bezier(thirdGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
			if (lambda > 0 && lambda < 1) {
				if (abs(Bezier(thirdGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(thirdGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) return true;
			}
		}
	}
	
	return false;
}			


struct handler {
	vector<BezierPiecewise*> borders;
	
	BezierPiecewise path;
	
	BezierPiecewise pathTemp0;
	
	BezierPiecewise pathTemp1;
	
	BezierPiecewise pathNext;
	
	const double multiplier = 1;
	
	vector<double> outsidex;
	vector<double> outsidey;
	vector<double> outsidet;
	
	double previousCondition = 100000000000000;
	
	void addBorder () {
		borders.push_back(new BezierPiecewise());
	}

	int collisionTest () {
		int collisionSum = 0;
		for (int bord = 0; bord < borders.size(); bord++) {
			for (int seg = 0; seg < path.segments.size(); seg++) {
				for (int bseg = 0; bseg < borders[bord]->segments.size(); bseg++) {
					collisionSum += Overlapping(path.segments[seg]->ax, path.segments[seg]->ay, path.segments[seg]->bx, path.segments[seg]->by, path.segments[seg]->cx, path.segments[seg]->cy, path.segments[seg]->dx, path.segments[seg]->dy, borders[bord]->segments[bseg]->ax, borders[bord]->segments[bseg]->ay, borders[bord]->segments[bseg]->bx, borders[bord]->segments[bseg]->by, borders[bord]->segments[bseg]->cx, borders[bord]->segments[bseg]->cy, borders[bord]->segments[bseg]->dx, borders[bord]->segments[bseg]->dy);
				}
			}
		}
		return collisionSum;
	}
	
	int collisionTestV2 () {
		int collisionSum = 0;
		for (int bord = 0; bord < borders.size(); bord++) {
			for (int seg = 0; seg < path.segments.size(); seg++) {
				for (int bseg = 0; bseg < borders[bord]->segments.size(); bseg++) {
					if (Overlapping(path.segments[seg]->ax, path.segments[seg]->ay, path.segments[seg]->bx, path.segments[seg]->by, path.segments[seg]->cx, path.segments[seg]->cy, path.segments[seg]->dx, path.segments[seg]->dy, borders[bord]->segments[bseg]->ax, borders[bord]->segments[bseg]->ay, borders[bord]->segments[bseg]->bx, borders[bord]->segments[bseg]->by, borders[bord]->segments[bseg]->cx, borders[bord]->segments[bseg]->cy, borders[bord]->segments[bseg]->dx, borders[bord]->segments[bseg]->dy)) {
						if (tinter.size() > 1) {
							double midx = 0;
							double midy = 0;
							double farx = 0;
							double fary = 0;
							midx = xinter[0] + 0.5 * (xinter[1] - xinter[0]);
							midy = yinter[0] + 0.5 * (yinter[1] - yinter[0]);
							farx = Bezier((tinter[0] + tinter[1]) / 2, path.segments[seg]->ax, path.segments[seg]->bx, path.segments[seg]->cx, path.segments[seg]->dx);
							fary = Bezier((tinter[0] + tinter[1]) / 2, path.segments[seg]->ay, path.segments[seg]->by, path.segments[seg]->cy, path.segments[seg]->dy);
							collisionSum += sqrt(pow(midx - farx, 2) + pow(midy - fary, 2));
						}
					} else return 0;
				}
			}
		}
		return collisionSum;
	}
	
	int collisionTestForward () {
		int collisionSum = 0;
		for (int bord = 0; bord < borders.size(); bord++) {
			for (int seg = 0; seg < pathTemp1.segments.size(); seg++) {
				for (int bseg = 0; bseg < borders[bord]->segments.size(); bseg++) {
					collisionSum += Overlapping(pathTemp1.segments[seg]->ax, pathTemp1.segments[seg]->ay, pathTemp1.segments[seg]->bx, pathTemp1.segments[seg]->by, pathTemp1.segments[seg]->cx, pathTemp1.segments[seg]->cy, pathTemp1.segments[seg]->dx, pathTemp1.segments[seg]->dy, borders[bord]->segments[bseg]->ax, borders[bord]->segments[bseg]->ay, borders[bord]->segments[bseg]->bx, borders[bord]->segments[bseg]->by, borders[bord]->segments[bseg]->cx, borders[bord]->segments[bseg]->cy, borders[bord]->segments[bseg]->dx, borders[bord]->segments[bseg]->dy);
				}
			}
		}
		return collisionSum;
	}
	
	int collisionTestForwardV2 () {
		int collisionSum = 0;
		for (int bord = 0; bord < borders.size(); bord++) {
			for (int seg = 0; seg < pathTemp1.segments.size(); seg++) {
				for (int bseg = 0; bseg < borders[bord]->segments.size(); bseg++) {
					if (Overlapping(pathTemp1.segments[seg]->ax, pathTemp1.segments[seg]->ay, pathTemp1.segments[seg]->bx, pathTemp1.segments[seg]->by, pathTemp1.segments[seg]->cx, pathTemp1.segments[seg]->cy, pathTemp1.segments[seg]->dx, pathTemp1.segments[seg]->dy, borders[bord]->segments[bseg]->ax, borders[bord]->segments[bseg]->ay, borders[bord]->segments[bseg]->bx, borders[bord]->segments[bseg]->by, borders[bord]->segments[bseg]->cx, borders[bord]->segments[bseg]->cy, borders[bord]->segments[bseg]->dx, borders[bord]->segments[bseg]->dy)) {
						if (tinter.size() > 1) {
							double midx = 0;
							double midy = 0;
							double farx = 0;
							double fary = 0;
							midx = xinter[0] + 0.5 * (xinter[1] - xinter[0]);
							midy = yinter[0] + 0.5 * (yinter[1] - yinter[0]);
							farx = Bezier((tinter[0] + tinter[1]) / 2, pathTemp1.segments[seg]->ax, pathTemp1.segments[seg]->bx, pathTemp1.segments[seg]->cx, pathTemp1.segments[seg]->dx);
							fary = Bezier((tinter[0] + tinter[1]) / 2, pathTemp1.segments[seg]->ay, pathTemp1.segments[seg]->by, pathTemp1.segments[seg]->cy, pathTemp1.segments[seg]->dy);
							collisionSum += sqrt(pow(midx - farx, 2) + pow(midy - fary, 2));
						}
					} else return 0;
				}
			}
		}
		return collisionSum;
	}
	
	int collisionTestBackward () {
		int collisionSum = 0;
		for (int bord = 0; bord < borders.size(); bord++) {
			for (int seg = 0; seg < pathTemp0.segments.size(); seg++) {
				for (int bseg = 0; bseg < borders[bord]->segments.size(); bseg++) {
					collisionSum += Overlapping(pathTemp0.segments[seg]->ax, pathTemp0.segments[seg]->ay, pathTemp0.segments[seg]->bx, pathTemp0.segments[seg]->by, pathTemp0.segments[seg]->cx, pathTemp0.segments[seg]->cy, pathTemp0.segments[seg]->dx, pathTemp0.segments[seg]->dy, borders[bord]->segments[bseg]->ax, borders[bord]->segments[bseg]->ay, borders[bord]->segments[bseg]->bx, borders[bord]->segments[bseg]->by, borders[bord]->segments[bseg]->cx, borders[bord]->segments[bseg]->cy, borders[bord]->segments[bseg]->dx, borders[bord]->segments[bseg]->dy);
				}
			}
		}
		return collisionSum;
	}
	
	int collisionTestBackwardV2 () {
		int collisionSum = 0;
		for (int bord = 0; bord < borders.size(); bord++) {
			for (int seg = 0; seg < pathTemp0.segments.size(); seg++) {
				for (int bseg = 0; bseg < borders[bord]->segments.size(); bseg++) {
					if (Overlapping(pathTemp0.segments[seg]->ax, pathTemp0.segments[seg]->ay, pathTemp0.segments[seg]->bx, pathTemp0.segments[seg]->by, pathTemp0.segments[seg]->cx, pathTemp0.segments[seg]->cy, pathTemp0.segments[seg]->dx, pathTemp0.segments[seg]->dy, borders[bord]->segments[bseg]->ax, borders[bord]->segments[bseg]->ay, borders[bord]->segments[bseg]->bx, borders[bord]->segments[bseg]->by, borders[bord]->segments[bseg]->cx, borders[bord]->segments[bseg]->cy, borders[bord]->segments[bseg]->dx, borders[bord]->segments[bseg]->dy)) {
						if (tinter.size() > 1) {
							double midx = 0;
							double midy = 0;
							double farx = 0;
							double fary = 0;
							midx = xinter[0] + 0.5 * (xinter[1] - xinter[0]);
							midy = yinter[0] + 0.5 * (yinter[1] - yinter[0]);
							farx = Bezier((tinter[0] + tinter[1]) / 2, pathTemp0.segments[seg]->ax, pathTemp0.segments[seg]->bx, pathTemp0.segments[seg]->cx, pathTemp0.segments[seg]->dx);
							fary = Bezier((tinter[0] + tinter[1]) / 2, pathTemp0.segments[seg]->ay, pathTemp0.segments[seg]->by, pathTemp0.segments[seg]->cy, pathTemp0.segments[seg]->dy);
							collisionSum += sqrt(pow(midx - farx, 2) + pow(midy - fary, 2));
						}
					} else return 0;
				}
			}
		}
		return collisionSum;
	}
	
	double fullEquation () {
		// This calculates the full equation at the current parameters
		return path.totalCurvature() * collisionTestV2() * multiplier;
	}

	double fullEquationDiff (int paramIndex, int segIndex) {
		// This calculates the derivative of the full equation at the current parameters
		if (segIndex < path.segments.size() - 1) {
			/*
			if (paramIndex == 0) {
				pathTemp0.segments[segIndex - 1]->params[4] = path.segments[segIndex - 1]->params[4] - step;
				//pathTemp0.segments[segIndex]->params[paramIndex] = path.segments[segIndex]->params[paramIndex] - step;
				pathTemp1.segments[segIndex - 1]->params[4] = path.segments[segIndex - 1]->params[4] + step;
			} else if (paramIndex == 1) {
				pathTemp0.segments[segIndex - 1]->params[5] = path.segments[segIndex - 1]->params[5] - step;
				pathTemp1.segments[segIndex - 1]->params[5] = path.segments[segIndex - 1]->params[5] + step;
			} else if (paramIndex == 2) {
				pathTemp0.segments[segIndex - 1]->params[6] = path.segments[segIndex - 1]->params[6] - step;
				pathTemp1.segments[segIndex - 1]->params[6] = path.segments[segIndex - 1]->params[6] + step;
			} else if (paramIndex == 3) {
				pathTemp0.segments[segIndex - 1]->params[7] = path.segments[segIndex - 1]->params[7] - step;
				pathTemp1.segments[segIndex - 1]->params[7] = path.segments[segIndex - 1]->params[7] + step;
			}
			*/
			
			if (paramIndex == 4) {
				pathTemp0.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] - step;
				pathTemp0.segments[segIndex + 1]->updateBX(path.segments[segIndex]->params[4] - step);
				pathTemp1.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] + step;
				pathTemp1.segments[segIndex + 1]->updateBX(path.segments[segIndex]->params[4] + step);
			} else if (paramIndex == 5) {
				pathTemp0.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] - step;
				pathTemp0.segments[segIndex + 1]->updateBY(path.segments[segIndex]->params[5] - step);
				pathTemp1.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] + step;
				pathTemp1.segments[segIndex + 1]->updateBY(path.segments[segIndex]->params[5] + step);
			} else if (paramIndex == 6) {
				pathTemp0.segments[segIndex + 1]->params[0] = path.segments[segIndex]->params[6] - step;
				pathTemp1.segments[segIndex + 1]->params[0] = path.segments[segIndex]->params[5] + step;
			} else if (paramIndex == 7) {
				pathTemp0.segments[segIndex + 1]->params[1] = path.segments[segIndex]->params[7] - step;
				pathTemp1.segments[segIndex + 1]->params[1] = path.segments[segIndex]->params[7] + step;
			} else if (paramIndex == 8) {
				pathTemp0.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] - step);
				pathTemp1.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] + step);
			}
		}
		
		else if (segIndex == path.segments.size() - 1) {
			if (paramIndex == 4) {
				pathTemp0.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] - step;
				pathTemp0.segments[0]->updateBX(path.segments[segIndex]->params[4] - step);
				pathTemp1.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] + step;
				pathTemp1.segments[0]->updateBX(path.segments[segIndex]->params[4] + step);
			} else if (paramIndex == 5) {
				pathTemp0.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] - step;
				pathTemp0.segments[0]->updateBY(path.segments[segIndex]->params[5] - step);
				pathTemp1.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] + step;
				pathTemp1.segments[0]->updateBY(path.segments[segIndex]->params[5] + step);
			} else if (paramIndex == 6) {
				pathTemp0.segments[0]->params[0] = path.segments[segIndex]->params[6] - step;
				pathTemp1.segments[0]->params[0] = path.segments[segIndex]->params[5] + step;
			} else if (paramIndex == 7) {
				pathTemp0.segments[0]->params[1] = path.segments[segIndex]->params[7] - step;
				pathTemp1.segments[0]->params[1] = path.segments[segIndex]->params[7] + step;
			} else if (paramIndex == 8) {
				pathTemp0.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] - step);
				pathTemp1.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] + step);
			}
		}
		
		else return 1;
		
		return (pathTemp1.totalCurvature() * exp(-collisionTestForwardV2())) - (pathTemp0.totalCurvature() * exp(-collisionTestBackwardV2()));
	}

	void newtonsParameterUpdate () {
		const int iterations = 10;
		
		for (int i = 0; i < iterations; i++) {
			cout << "Iteration: " << i + 1 << endl;
			for (int seg = 0; seg < path.segments.size(); seg++) {
				cout << (float(seg) / float(path.segments.size())) * 100 << "%\n";
				if (seg == path.segments.size() - 1) {
					for (int param = 4; param <= 8; param++) {
						path.segments[seg]->params[param] = path.segments[seg]->params[param] - ((fullEquation () * 2 * step) / fullEquationDiff(param, seg));
						if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
						else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
						else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
						else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
						else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
					}
				} else {
					for (int param = 4; param <= 8; param++) {
						path.segments[seg]->params[param] = path.segments[seg]->params[param] - ((fullEquation () * 2 * step) / fullEquationDiff(param, seg));
						if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
						else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
						else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
						else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
						else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
					}
				}
			}
		}
		
		cout << endl;
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param <= 8; param++) {
				cout << path.segments[seg]->params[param] << ", ";
			} cout << endl;
		}
	}
	
	void saveToFile () {
		// Save borders
		for (int border = 0; border < borders.size(); border++) {
			ofstream borderFile("borders/border" + to_string(border) + ".txt");
			for (int seg = 0; seg < borders[border]->segments.size(); seg++) {
				for (int param = 0; param < 8; param++) {
					borderFile << borders[border]->segments[seg]->params[param];
					if (param < 7) borderFile << " ";
				}
				if (seg < borders[border]->segments.size() - 1) borderFile << endl;
			}
			borderFile.close();
		}
		
		// Save path
		ofstream pathFile("paths/path" + to_string(0) + ".txt");
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param < 8; param++) {
				pathFile << path.segments[seg]->params[param];
				if (param < 7) pathFile << " ";
			}
			if (seg < path.segments.size() - 1) pathFile << endl;
		}
		pathFile.close();
		
		// Save outside points
		ofstream outFile("outside/outside.txt");
		for (int p = 0; p < outsidex.size(); p++) {
			outFile << outsidex[p] << " " << outsidey[p];
			if (p < outsidex.size() - 1) outFile << endl;
		}
		outFile.close();
		
		/*
		// Save outside point times
		ofstream outTFile("outside/outsideT.txt");
		for (int p = 0; p < outsidet.size(); p++) {
			outTFile << outsidet[p];
			if (p < outsidet.size() - 1) outTFile << endl;
		}
		outTFile.close();
		*/
		
		// Save pixels
		ofstream bPixelsFile("pixelsBorder/pixelsBorder.txt");
		
		for (int border = 0; border < borders.size(); border++) {
			for (int seg = 0; seg < borders[border]->segments.size(); seg++) {
				for (int pixel = 0; pixel < borders[border]->segments[seg]->tCoord.size(); pixel++) {
					bPixelsFile << borders[border]->segments[seg]->xCoord[pixel] << " " << borders[border]->segments[seg]->yCoord[pixel];
					if (pixel == borders[border]->segments[seg]->tCoord.size() - 1 && seg == borders[border]->segments.size() - 1 && border == borders.size() - 1) cout << ".";
					else bPixelsFile << endl;
				}
			}
		}
		bPixelsFile.close();
		
		ofstream pPixelsFile("pixelsPath/pixelsPath.txt");
		
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int pixel = 0; pixel < path.segments[seg]->tCoord.size(); pixel++) {
				pPixelsFile << path.segments[seg]->xCoord[pixel] << " " << path.segments[seg]->yCoord[pixel];
				if (pixel == path.segments[seg]->tCoord.size() - 1 && seg == path.segments.size() - 1) cout << ".";
				else pPixelsFile << endl;
			}
		}
		
		pPixelsFile.close();
		
		
	}
	
	/*
	// return true if intersecting
	bool OverlappingV3 (int pathSeg, int border, int borderSeg) {
		double x0 = 0;
		double x1 = 0;

		double y0 = 0;
		double y1 = 0;

		double lambda = 0;

		const int iterations = 100;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.001;

		double previousGuess = 0;
		double nextGuess = 0;

		double firstGuess = 0;
		double secondGuess = 0;
		double thirdGuess = 0;

		xinter.clear();
		yinter.clear();
		tinter.clear();

		for (int i = 1; i < 1 / (step); i++) {
			x0 = Bezier(i * (step), path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
			x1 = Bezier((i + 1) * (step), path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);

			y0 = Bezier(i * (step), path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
			y1 = Bezier((i + 1) * (step), path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);

			// Newtons from 0.0001
			nextGuess = 0.0001;
			for (int n = 0; n < iterations; n++) {
				previousGuess = nextGuess;
				nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
				//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
				if (abs(previousGuess - nextGuess) < minimumDifference) break;
			}

			// Store first guess
			firstGuess = nextGuess;

			// Newtons from 0.5
			nextGuess = 0.5;
			for (int n = 0; n < iterations; n++) {
				previousGuess = nextGuess;
				nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
				//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
				if (abs(previousGuess - nextGuess) < minimumDifference) break;
			}

			// Store second guess
			secondGuess = nextGuess;

			// Newtons from 0.9999
			nextGuess = 0.9999;
			for (int n = 0; n < iterations; n++) {
				previousGuess = nextGuess;
				nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
				//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
				if (abs(previousGuess - nextGuess) < minimumDifference) break;
			}

			// Store third guess
			thirdGuess = nextGuess;	


			// Test different Lambdas		
			if (firstGuess > 0 && firstGuess < 1) {
				lambda = (Bezier(firstGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
				//cout << firstGuessRight << endl;
				//cout << lambda << endl;
				if (lambda > 0 && lambda < 1) {
					if (abs(Bezier(firstGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(firstGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
						//cout << firstGuess << "t, ";
						tinter.push_back((i + 0.5) * step);
						xinter.push_back((x1 - x0) / 2);
						yinter.push_back((y1 - y0) / 2);
						return true;
					}
				}
			}

			if (secondGuess > 0 && secondGuess < 1) {
				lambda = (Bezier(secondGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
				//cout << secondGuessRight << endl;
				//cout << lambda << endl;
				if (lambda > 0 && lambda < 1) {
					if (abs(Bezier(secondGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(secondGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
						//cout << secondGuess << "t, ";
						tinter.push_back((i + 0.5) * step);
						xinter.push_back((x1 - x0) / 2);
						yinter.push_back((y1 - y0) / 2);
						return true;
					}
				}
			}

			if (thirdGuess > 0 && thirdGuess < 1) {
				lambda = (Bezier(thirdGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
				//cout << thirdGuessRight << endl;
				//cout << lambda << endl;
				if (lambda > 0 && lambda < 1) {
					if (abs(Bezier(thirdGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError && abs(Bezier(thirdGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
						//cout << thirdGuess << "t, ";
						tinter.push_back((i + 0.5) * step);
						xinter.push_back((x1 - x0) / 2);
						yinter.push_back((y1 - y0) / 2);
						return true;
					}
				}
			}
		}

		return false;
	}
	*/
	
	// return external area
	double OverlappingV3 () {
		outsidex.clear();
		outsidey.clear();
						
		double x0 = 0;
		double x1 = 0;

		double y0 = 0;
		double y1 = 0;

		double lambda = 0;

		const int iterations = 100;
		const double acceptableError = 0.1;
		const double minimumDifference = 0.0001;
		
		double artificialArea = 0;

		double previousGuess = 0;
		double nextGuess = 0;

		double firstGuess = 0;
		double secondGuess = 0;
		double thirdGuess = 0;

		double area = 0;
		
		int guessSections = 4;

		for (int border = 0; border < borders.size(); border++) {
			vector <double> uniqueGuessesPath;
			vector <double> uniqueGuessesBorder;

			for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
				for (int pathSeg = 0; pathSeg < path.segments.size(); pathSeg++) {
					for (int i = 0; i < 1 / (step); i++) {
						x0 = Bezier(i * (interStep), path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
						x1 = Bezier((i + 1) * (interStep), path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);

						y0 = Bezier(i * (interStep), path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
						y1 = Bezier((i + 1) * (interStep), path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
						
						vector <double> guesses;
						
						for (int g = 0; g <= guessSections; g++) {
							nextGuess = double(g) / double(guessSections);

							for (int n = 0; n < iterations; n++) {
								previousGuess = nextGuess;
								nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
								if (abs(previousGuess - nextGuess) < minimumDifference) break;
							}

							// Store the guess
							if (nextGuess > 0 && nextGuess < 1)	guesses.push_back(nextGuess);
						}
						
						// Sort the Guesses
						sort(guesses.begin(), guesses.end());
						
						// Filter for Valid Guesses
						vector <double> validGuessesPath;
						vector <double> validGuessesBorder;
						
						for (int guess = 0; guess < guesses.size(); guess++) {
							
							//cout << 0 << endl;
							
							// Test different Lambdas		
							if (guesses[guess] > 0 && guesses[guess] < 1) {
								
								//cout << 1 << endl;
								
								lambda = (Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
								if (lambda >= 0 && lambda < 1) {
									
									//cout << 2 << endl;
									
									if (abs(Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
										//cout << 3 << endl;
										
										//cout << borderSeg << endl;
										
										//cout << guess << endl;
										
										//cout << guesses[guess] << endl;
										
										validGuessesBorder.push_back(borderSeg + guesses[guess]);
										
										//cout << 4 << endl;
										
										validGuessesPath.push_back(pathSeg + ((i * step) + ((lambda * step) / sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)))));
										
										//cout << 5 << endl;
									}
								}
							}
						}

						// Filter for Unique Guesses
												
						int precision = 5;
						bool unique = true;

						for (int guess = 0; guess < validGuessesBorder.size(); guess++) {
							for (int uniqueGuess = 0; uniqueGuess < uniqueGuessesBorder.size(); uniqueGuess++) {
								if (ceil(validGuessesBorder[guess] * pow(10, precision)) == ceil(uniqueGuessesBorder[uniqueGuess] * pow(10, precision))) unique = false;
							}

							if (unique) {
								uniqueGuessesPath.push_back(validGuessesPath[guess]);
								uniqueGuessesBorder.push_back(validGuessesBorder[guess]);
							}
						}
					}
				}
			}
			
			cout << uniqueGuessesPath.size() << " unique guesses\n";
						
			// Find corresponding exits and re-entries		
			vector <double> validInterPath;
			vector <double> validInterBorder;
			
			bool lookingForExit = true;
			//double tempPathT = uniqueGuessesPath[0];
			double tempBorderT = uniqueGuessesBorder[0];
			
			double borderStep = 0;
			double pathStep = 0;
			double pathTempX = 0;
			double pathTempY = 0;
			double borderTempX = 0;
			double borderTempY = 0;
			int pathTempSeg = 0;
			int borderTempSeg = 0;
			double pathTempT = 0;
			double borderTempT = 0;
			int areaParts = 100;
			
			if (uniqueGuessesPath.size() > 0) {
				// Add the first intersection, take it as an exit
				validInterPath.push_back(uniqueGuessesPath[0]);
				validInterBorder.push_back(uniqueGuessesBorder[0]);
				
				for (int guess = 0; guess < uniqueGuessesBorder.size(); guess++) {
					pathTempSeg = floor(uniqueGuessesPath[guess]);
					pathTempT = uniqueGuessesPath[guess] - pathTempSeg;
					
					if (lookingForExit) {
						if (uniqueGuessesBorder[guess] > tempBorderT) {
							if (pointInsideSegmentNewtons(pathTempSeg, pathTempT + 0.1) || pointInsideSegmentNewtons(pathTempSeg, pathTempT + 0.05)) {
								if (!pointInsideSegmentNewtons(pathTempSeg, pathTempT - 0.1) || !pointInsideSegmentNewtons(pathTempSeg, pathTempT - 0.05)) {
									validInterPath.push_back(uniqueGuessesPath[guess]);
									validInterBorder.push_back(uniqueGuessesBorder[guess]);

									lookingForExit = false;
								}
							}
						}
					} else {
						if (!pointInsideSegmentNewtons(pathTempSeg, pathTempT + 0.1) || !pointInsideSegmentNewtons(pathTempSeg, pathTempT + 0.05)) {
							if (pointInsideSegmentNewtons(pathTempSeg, pathTempT - 0.1) || pointInsideSegmentNewtons(pathTempSeg, pathTempT - 0.05)) {
								validInterPath.push_back(uniqueGuessesPath[guess]);
								validInterBorder.push_back(uniqueGuessesBorder[guess]);
								tempBorderT = uniqueGuessesBorder[guess];
								lookingForExit = true;
							}
						}
					}
				}
			}
			
			// Find Area between exits and re-entries
			if (validInterPath.size() % 2 > 0) cout << "Odd amount of Intersections" << endl;

			for (int inter = 0; inter < validInterBorder.size() / 2; inter++) {
				borderStep = (validInterBorder[inter*2 + 1] - validInterBorder[inter*2]) / areaParts;
				pathStep = (validInterPath[inter*2 + 1] - validInterPath[inter*2]) / areaParts;

				for (int areaSeg = 1; areaSeg < areaParts; areaSeg++) {
					pathTempSeg = floor(validInterPath[inter*2] + (areaSeg * pathStep));
					borderTempSeg = floor(validInterBorder[inter*2] + (areaSeg * borderStep));
					pathTempT = (floor(validInterPath[inter*2]) + (areaSeg * pathStep)) - pathTempSeg;
					borderTempT = (floor(validInterBorder[inter*2]) + (areaSeg * borderStep)) - borderTempSeg;
					pathTempX = Bezier(pathTempT, path.segments[pathTempSeg]->params[0], path.segments[pathTempSeg]->params[2], path.segments[pathTempSeg]->params[4], path.segments[pathTempSeg]->params[6]);
					pathTempY = Bezier(pathTempT, path.segments[pathTempSeg]->params[1], path.segments[pathTempSeg]->params[3], path.segments[pathTempSeg]->params[5], path.segments[pathTempSeg]->params[7]);
					borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
					borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);
					artificialArea += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2));
					outsidex.push_back(pathTempX);
					outsidey.push_back(pathTempY);
				}
			}
		}
		
		return artificialArea;
	}
	
	bool pointInsideSegmentNewtons (int pathSeg, double t) {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		
		pointx = Bezier(t, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
		pointy = Bezier(t, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
		vector<double> distance0;
		vector<double> distance1;

		for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {

			for (int i = 1; i < parts; i++) {
				nextGuess = double(i) / double(parts);
				for (int n = 0; n < iterations; n++) {
					previousGuess = nextGuess;

					Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
					fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

					if (fx == 0) break;

					nextGuess = previousGuess - (Fx / fx);

					if (abs(previousGuess - nextGuess) < minimumDifference) break;
				}

				if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
					pointsOutside++;
					return true;
				}
			}
		}

		return false;
	}
	
	double Area (int pathIndex, int border) {
		double x0 = 0;
		double x1 = 0;

		double y0 = 0;
		double y1 = 0;

		double lambda = 0;

		const int iterations = 100;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.001;

		double previousGuess = 0;
		double nextGuess = 0;

		double firstGuess = 0;
		double secondGuess = 0;
		double thirdGuess = 0;
		
		vector<double> pathInterT;
		//vector<double> pathInterSeg;
		
		vector<double> borderInterT;
		//vector<double> borderInterSeg;
		
		//cout << "area" << endl;
		if (pathIndex == 0) {
			//cout << "path 0" << endl;
			for (int pathSeg = 0; pathSeg < path.segments.size(); pathSeg++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int i = 1; i < 1 / (interStep); i++) {
						x0 = Bezier(i * (interStep), path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
						x1 = Bezier((i + 1) * (interStep), path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);

						y0 = Bezier(i * (interStep), path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
						y1 = Bezier((i + 1) * (interStep), path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);

						// Newtons from 0.0001
						nextGuess = 0.0001;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store first guess
						firstGuess = nextGuess;

						// Newtons from 0.5
						nextGuess = 0.5;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store second guess
						secondGuess = nextGuess;

						// Newtons from 0.9999
						nextGuess = 0.9999;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store third guess
						thirdGuess = nextGuess;

						vector<double> uniqueGuesses;

						if (abs(firstGuess - secondGuess) < 2 * minimumDifference && abs(firstGuess - thirdGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference) {
							uniqueGuesses.push_back((firstGuess + secondGuess + thirdGuess) / 3.0);
						} else if (abs(firstGuess - secondGuess) < 2 * minimumDifference && !(abs(firstGuess - thirdGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
							uniqueGuesses.push_back(thirdGuess);
							/*
							if (firstGuess < thirdGuess || secondGuess < thirdGuess) {
								uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
								uniqueGuesses.push_back(thirdGuess);
							} else {
								uniqueGuesses.push_back(thirdGuess);
								uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
							}
							*/
						} else if (abs(firstGuess - thirdGuess) < 2 * minimumDifference && !(abs(firstGuess - secondGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
							uniqueGuesses.push_back(secondGuess);
							/*
							if (firstGuess < secondGuess || thirdGuess < secondGuess) {
								uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
								uniqueGuesses.push_back(secondGuess);
							} else {
								uniqueGuesses.push_back(secondGuess);
								uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
							}
							*/
						} else if (abs(secondGuess - thirdGuess) < 2 * minimumDifference && !(abs(secondGuess - firstGuess) < 2 * minimumDifference && abs(thirdGuess - firstGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
							uniqueGuesses.push_back(firstGuess);
							/*
							if (secondGuess < firstGuess || thirdGuess < firstGuess) {
								uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
								uniqueGuesses.push_back(firstGuess);
							} else {
								uniqueGuesses.push_back(firstGuess);
								uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
							}
							*/
						} else {
							uniqueGuesses.push_back(firstGuess);
							uniqueGuesses.push_back(secondGuess);
							uniqueGuesses.push_back(thirdGuess);
							/*
							if (firstGuess < secondGuess && secondGuess < thirdGuess) {
								uniqueGuesses.push_back(firstGuess);
								uniqueGuesses.push_back(secondGuess);
								uniqueGuesses.push_back(thirdGuess);
							else if (
							*/
						}

						// Sort the Guesses
						sort(uniqueGuesses.begin(), uniqueGuesses.end());

						if (uniqueGuesses.size() <= 0) cout << "Guesses Error" << endl;
						else {
							for (int guess = 0; guess < uniqueGuesses.size(); guess++) {
								// Test different Lambdas		
								if (uniqueGuesses[guess] > 0 && uniqueGuesses[guess] < 1) {
									lambda = (Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
									if (lambda > 0 && lambda < 1) {
										if (abs(Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
											borderInterT.push_back(borderSeg + uniqueGuesses[guess]);
											//borderInterSeg.push_back(borderSeg);

											pathInterT.push_back(pathSeg + ((i * step) + ((lambda * step) / sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)))));
											//pathInterSeg.push_back(pathSeg);
										}
									}
								}
							}
						}
					}
				}
			}

			// Find Area Outside Track
			//if (pathInterT.size() % 2 > 0) cout << "Odd amount of Intersections" << endl;
			//if (pathInterT.size() > 0) cout << pathInterT.size() << " intersections" << endl;

			double borderStep = 0;
			double pathStep = 0;
			double artificialArea = 0;
			double pathTempX = 0;
			double pathTempY = 0;
			double borderTempX = 0;
			double borderTempY = 0;
			int pathTempSeg = 0;
			int borderTempSeg = 0;
			double pathTempT = 0;
			double borderTempT = 0;

			for (int inter = 0; inter < borderInterT.size() / 2; inter++) {
				borderStep = (borderInterT[inter*2 + 1] - borderInterT[inter*2]) / parts;
				pathStep = (pathInterT[inter*2 + 1] - pathInterT[inter*2]) / parts;

				for (int areaSeg = 1; areaSeg < parts; areaSeg++) {
					pathTempSeg = floor(pathInterT[inter*2] + (areaSeg * pathStep));
					borderTempSeg = floor(borderInterT[inter*2] + (areaSeg * borderStep));

					pathTempT = (floor(pathInterT[inter*2]) + (areaSeg * pathStep)) - pathTempSeg;
					borderTempT = (floor(borderInterT[inter*2]) + (areaSeg * borderStep)) - borderTempSeg;

					pathTempX = Bezier(pathTempT, path.segments[pathTempSeg]->params[0], path.segments[pathTempSeg]->params[2], path.segments[pathTempSeg]->params[4], path.segments[pathTempSeg]->params[6]);
					pathTempY = Bezier(pathTempT, path.segments[pathTempSeg]->params[1], path.segments[pathTempSeg]->params[3], path.segments[pathTempSeg]->params[5], path.segments[pathTempSeg]->params[7]);

					borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
					borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);

					artificialArea += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2));
				}
			}

			return artificialArea;
		}
		
		else if (pathIndex == 1) {
			//cout << "path 1" << endl;
			for (int pathSeg = 0; pathSeg < pathTemp0.segments.size(); pathSeg++) {
				//cout << "pathSeg " << pathSeg << endl;
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					//cout << "borderSeg " << borderSeg << endl;
					for (int i = 1; i < 1 / (interStep); i++) {
						x0 = Bezier(i * (interStep), pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[6]);
						x1 = Bezier((i + 1) * (interStep), pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[6]);

						y0 = Bezier(i * (interStep), pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[7]);
						y1 = Bezier((i + 1) * (interStep), pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[7]);

						// Newtons from 0.0001
						nextGuess = 0.0001;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store first guess
						firstGuess = nextGuess;

						// Newtons from 0.5
						nextGuess = 0.5;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store second guess
						secondGuess = nextGuess;

						// Newtons from 0.9999
						nextGuess = 0.9999;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store third guess
						thirdGuess = nextGuess;

						vector<double> uniqueGuesses;

						if (abs(firstGuess - secondGuess) < 2 * minimumDifference && abs(firstGuess - thirdGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference) {
							uniqueGuesses.push_back((firstGuess + secondGuess + thirdGuess) / 3.0);
						} else if (abs(firstGuess - secondGuess) < 2 * minimumDifference && !(abs(firstGuess - thirdGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
							uniqueGuesses.push_back(thirdGuess);
							/*
							if (firstGuess < thirdGuess || secondGuess < thirdGuess) {
								uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
								uniqueGuesses.push_back(thirdGuess);
							} else {
								uniqueGuesses.push_back(thirdGuess);
								uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
							}
							*/
						} else if (abs(firstGuess - thirdGuess) < 2 * minimumDifference && !(abs(firstGuess - secondGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
							uniqueGuesses.push_back(secondGuess);
							/*
							if (firstGuess < secondGuess || thirdGuess < secondGuess) {
								uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
								uniqueGuesses.push_back(secondGuess);
							} else {
								uniqueGuesses.push_back(secondGuess);
								uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
							}
							*/
						} else if (abs(secondGuess - thirdGuess) < 2 * minimumDifference && !(abs(secondGuess - firstGuess) < 2 * minimumDifference && abs(thirdGuess - firstGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
							uniqueGuesses.push_back(firstGuess);
							/*
							if (secondGuess < firstGuess || thirdGuess < firstGuess) {
								uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
								uniqueGuesses.push_back(firstGuess);
							} else {
								uniqueGuesses.push_back(firstGuess);
								uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
							}
							*/
						} else {
							uniqueGuesses.push_back(firstGuess);
							uniqueGuesses.push_back(secondGuess);
							uniqueGuesses.push_back(thirdGuess);
							/*
							if (firstGuess < secondGuess && secondGuess < thirdGuess) {
								uniqueGuesses.push_back(firstGuess);
								uniqueGuesses.push_back(secondGuess);
								uniqueGuesses.push_back(thirdGuess);
							else if (
							*/
						}
						
						// Sort the Guesses
						sort(uniqueGuesses.begin(), uniqueGuesses.end());

						if (uniqueGuesses.size() <= 0) cout << "Guesses Error" << endl;
						else {
							for (int guess = 0; guess < uniqueGuesses.size(); guess++) {
								// Test different Lambdas		
								if (uniqueGuesses[guess] > 0 && uniqueGuesses[guess] < 1) {
									lambda = (Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
									if (lambda > 0 && lambda < 1) {
										if (abs(Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
											borderInterT.push_back(borderSeg + uniqueGuesses[guess]);
											//borderInterSeg.push_back(borderSeg);

											pathInterT.push_back(pathSeg + ((i * step) + ((lambda * step) / sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)))));
											//pathInterSeg.push_back(pathSeg);
										}
									}
								}
							}
						}
					}
				}
			}
			
			//cout << "Area Calc Time" << endl;

			// Find Area Outside Track
			if (pathInterT.size() % 2 > 0) cout << "Odd amount of Intersections" << endl;

			double borderStep = 0;
			double pathStep = 0;
			double artificialArea = 0;
			double pathTempX = 0;
			double pathTempY = 0;
			double borderTempX = 0;
			double borderTempY = 0;
			int pathTempSeg = 0;
			int borderTempSeg = 0;
			double pathTempT = 0;
			double borderTempT = 0;

			for (int inter = 0; inter < borderInterT.size() / 2; inter++) {
				borderStep = (borderInterT[inter*2 + 1] - borderInterT[inter*2]) / parts;
				pathStep = (pathInterT[inter*2 + 1] - pathInterT[inter*2]) / parts;

				for (int areaSeg = 1; areaSeg < parts; areaSeg++) {
					//cout << pathStep << endl;
					//cout << 1 << endl;
					pathTempSeg = floor(pathInterT[inter*2] + (areaSeg * pathStep));
					//cout << 2 << " PathTempSeg " << pathTempSeg << endl;
					borderTempSeg = floor(borderInterT[inter*2] + (areaSeg * borderStep));
					//cout << 3  << " BorderTempSeg " << borderTempSeg << endl;
					pathTempT = (floor(pathInterT[inter*2]) + (areaSeg * pathStep)) - pathTempSeg;
					//cout << 4 << " PathTempT " << pathTempT << endl;
					borderTempT = (floor(borderInterT[inter*2]) + (areaSeg * borderStep)) - borderTempSeg;
					//cout << 5 << " BorderTempT " << borderTempT << endl;
					pathTempX = Bezier(pathTempT, pathTemp1.segments[pathTempSeg]->params[0], pathTemp1.segments[pathTempSeg]->params[2], pathTemp1.segments[pathTempSeg]->params[4], pathTemp1.segments[pathTempSeg]->params[6]);
					//cout << 6 << endl;
					pathTempY = Bezier(pathTempT, pathTemp1.segments[pathTempSeg]->params[1], pathTemp1.segments[pathTempSeg]->params[3], pathTemp1.segments[pathTempSeg]->params[5], pathTemp1.segments[pathTempSeg]->params[7]);
					//cout << 7 << endl;
					borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
					//cout << 8 << endl;
					borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);
					//cout << 9 << endl;
					artificialArea += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2));
					//cout << 10 << endl;
				}
			}
			
			//cout << "Artificial Area : " << artificialArea << endl;

			return artificialArea;
		}
		
		else if (pathIndex == -1) {
			//cout << "Path -1" << endl;
			for (int pathSeg = 0; pathSeg < pathTemp0.segments.size(); pathSeg++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int i = 1; i < 1 / (interStep); i++) {
						x0 = Bezier(i * (interStep), pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[6]);
						x1 = Bezier((i + 1) * (interStep), pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[6]);

						y0 = Bezier(i * (interStep), pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[7]);
						y1 = Bezier((i + 1) * (interStep), pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[7]);

						// Newtons from 0.0001
						nextGuess = 0.0001;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store first guess
						firstGuess = nextGuess;

						// Newtons from 0.5
						nextGuess = 0.5;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store second guess
						secondGuess = nextGuess;

						// Newtons from 0.9999
						nextGuess = 0.9999;
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;
							nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
							//nextGuess = previousGuess - (((2 * step) * (Bezier(previousGuess, rax, rbx, rcx, rdx) - (Bezier(previousGuess, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))) / ((Bezier(previousGuess + step, rax, rbx, rcx, rdx) - (Bezier(previousGuess + step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0))))) - (Bezier(previousGuess - step, rax, rbx, rcx, rdx) - (Bezier(previousGuess - step, ray, rby, rcy, rdy) * ((x1 - x0) / (y1 - y0))) - (x0 - (y0 * ((x1 - x0) / (y1 - y0)))))));
							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}

						// Store third guess
						thirdGuess = nextGuess;

						vector<double> uniqueGuesses;

						if (abs(firstGuess - secondGuess) < 2 * minimumDifference && abs(firstGuess - thirdGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference) {
							uniqueGuesses.push_back((firstGuess + secondGuess + thirdGuess) / 3.0);
						} else if (abs(firstGuess - secondGuess) < 2 * minimumDifference && !(abs(firstGuess - thirdGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
							uniqueGuesses.push_back(thirdGuess);
							/*
							if (firstGuess < thirdGuess || secondGuess < thirdGuess) {
								uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
								uniqueGuesses.push_back(thirdGuess);
							} else {
								uniqueGuesses.push_back(thirdGuess);
								uniqueGuesses.push_back((firstGuess + secondGuess) / 2.0);
							}
							*/
						} else if (abs(firstGuess - thirdGuess) < 2 * minimumDifference && !(abs(firstGuess - secondGuess) < 2 * minimumDifference && abs(thirdGuess - secondGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
							uniqueGuesses.push_back(secondGuess);
							/*
							if (firstGuess < secondGuess || thirdGuess < secondGuess) {
								uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
								uniqueGuesses.push_back(secondGuess);
							} else {
								uniqueGuesses.push_back(secondGuess);
								uniqueGuesses.push_back((firstGuess + thirdGuess) / 2.0);
							}
							*/
						} else if (abs(secondGuess - thirdGuess) < 2 * minimumDifference && !(abs(secondGuess - firstGuess) < 2 * minimumDifference && abs(thirdGuess - firstGuess) < 2 * minimumDifference)) {
							uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
							uniqueGuesses.push_back(firstGuess);
							/*
							if (secondGuess < firstGuess || thirdGuess < firstGuess) {
								uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
								uniqueGuesses.push_back(firstGuess);
							} else {
								uniqueGuesses.push_back(firstGuess);
								uniqueGuesses.push_back((secondGuess + thirdGuess) / 2.0);
							}
							*/
						} else {
							uniqueGuesses.push_back(firstGuess);
							uniqueGuesses.push_back(secondGuess);
							uniqueGuesses.push_back(thirdGuess);
							/*
							if (firstGuess < secondGuess && secondGuess < thirdGuess) {
								uniqueGuesses.push_back(firstGuess);
								uniqueGuesses.push_back(secondGuess);
								uniqueGuesses.push_back(thirdGuess);
							else if (
							*/
						}

						// Sort the Guesses
						sort(uniqueGuesses.begin(), uniqueGuesses.end());

						if (uniqueGuesses.size() <= 0) cout << "Guesses Error" << endl;
						else {
							for (int guess = 0; guess < uniqueGuesses.size(); guess++) {
								// Test different Lambdas		
								if (uniqueGuesses[guess] > 0 && uniqueGuesses[guess] < 1) {
									lambda = (Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
									if (lambda > 0 && lambda < 1) {
										if (abs(Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(uniqueGuesses[guess], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
											borderInterT.push_back(borderSeg + uniqueGuesses[guess]);
											//borderInterSeg.push_back(borderSeg);

											pathInterT.push_back(pathSeg + ((i * step) + ((lambda * step) / sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)))));
											//pathInterSeg.push_back(pathSeg);
										}
									}
								}
							}
						}
					}
				}
			}
			
			//cout << "Valid Guesses number: " << pathInterT.size() << endl;

			// Find Area Outside Track
			//if (pathInterT.size() % 2 > 0) cout << "Odd amount of Intersections" << endl;

			double borderStep = 0;
			double pathStep = 0;
			double artificialArea = 0;
			double pathTempX = 0;
			double pathTempY = 0;
			double borderTempX = 0;
			double borderTempY = 0;
			int pathTempSeg = 0;
			int borderTempSeg = 0;
			double pathTempT = 0;
			double borderTempT = 0;

			for (int inter = 0; inter < borderInterT.size() / 2; inter++) {
				borderStep = (borderInterT[inter*2 + 1] - borderInterT[inter*2]) / parts;
				pathStep = (pathInterT[inter*2 + 1] - pathInterT[inter*2]) / parts;

				for (int areaSeg = 1; areaSeg < parts; areaSeg++) {
					pathTempSeg = floor(pathInterT[inter*2] + (areaSeg * pathStep));
					//cout << 2  << " pathTempSeg " << pathTempSeg << endl;
					borderTempSeg = floor(borderInterT[inter*2] + (areaSeg * borderStep));
					//cout << 3  << " BorderTempSeg " << borderTempSeg << endl;
					pathTempT = (floor(pathInterT[inter*2]) + (areaSeg * pathStep)) - pathTempSeg;
					borderTempT = (floor(borderInterT[inter*2]) + (areaSeg * borderStep)) - borderTempSeg;
				
					pathTempX = Bezier(pathTempT, pathTemp0.segments[pathTempSeg]->params[0], pathTemp0.segments[pathTempSeg]->params[2], pathTemp0.segments[pathTempSeg]->params[4], pathTemp0.segments[pathTempSeg]->params[6]);
					pathTempY = Bezier(pathTempT, pathTemp0.segments[pathTempSeg]->params[1], pathTemp0.segments[pathTempSeg]->params[3], pathTemp0.segments[pathTempSeg]->params[5], pathTemp0.segments[pathTempSeg]->params[7]);

					borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
					borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);

					artificialArea += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2));
				}
			}
			
			//cout << "Artificial Area : " << artificialArea << endl;

			return artificialArea;
		}
		
		return 0;
					/*
					// Test different Lambdas		
					if (firstGuess > 0 && firstGuess < 1) {
						lambda = (Bezier(firstGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
						//cout << firstGuessRight << endl;
						//cout << lambda << endl;
						if (lambda > 0 && lambda < 1) {
							if (abs(Bezier(firstGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(firstGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
								
							}
						}
					}

					if (secondGuess > 0 && secondGuess < 1) {
						lambda = (Bezier(secondGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
						//cout << secondGuessRight << endl;
						//cout << lambda << endl;
						if (lambda > 0 && lambda < 1) {
							if (abs(Bezier(secondGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(secondGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
								//cout << secondGuess << "t, ";
							}
						}
					}

					if (thirdGuess > 0 && thirdGuess < 1) {
						lambda = (Bezier(thirdGuess, ray, rby, rcy, rdy) - y0) / (y1 - y0);
						//cout << thirdGuessRight << endl;
						//cout << lambda << endl;
						if (lambda > 0 && lambda < 1) {
							if (abs(Bezier(thirdGuess, ray, rby, rcy, rdy) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(thirdGuess, rax, rbx, rcx, rdx) - (x0 + (lambda*(x1 - x0)))) < acceptableError * lambda) {
								//cout << thirdGuess << "t, ";
							}
						}
					}
					*/
	}
	
	double mainFunction () {
		double sum = 0;
		//double totalArea = 0;
		
		//cout << "starting mainFunction" << endl;
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int i = 0; i < 1 / step; i++) {
				sum += path.segments[seg]->sDiff(i * step) * exp(k * path.segments[seg]->CurvatureV2(i * step));
			}
		}
		
		/*
		for (int border = 0; border < borders.size(); border++) {
			totalArea += Area(0, border);
		}
		*/
		
		//cout << "mainFunction Complete" << endl;
		//cout << "Total Area : " << totalArea * multiplier << endl;
		
		cout << "MF : " << sum / Vmax << endl;
				
		return (sum / Vmax); // + (totalArea * multiplier);
		//return 0;
	}
	
	double mainFunctionBackward () {
		double sum = 0;
		//double totalAreaTemp0 = 0;
		
		//cout << "starting mainFunctionBackward" << endl;
		for (int seg = 0; seg < pathTemp0.segments.size(); seg++) {
			for (int i = 0; i < 1 / step; i++) {
				sum += pathTemp0.segments[seg]->sDiff(i * step) * exp(k * pathTemp0.segments[seg]->CurvatureV2(i * step));
			}
		}
		
		/*
		for (int border = 0; border < borders.size(); border++) {
			totalAreaTemp0 += Area(-1, border);
		}
		*/
		
		//cout << "mainFunctionBackward Complete" << endl;
		//cout << "Total Area Backward : " << totalAreaTemp0 * multiplier << endl;=
		
		cout << "MF Backward : " << sum / Vmax << endl;
				
		return (sum / Vmax); // + (totalAreaTemp0 * multiplier);
		//return 0;
	}
	
	double mainFunctionForward () {
		double sum = 0;
		//double totalAreaTemp1 = 0;
		
		//cout << "starting mainFunctionForward" << endl;
		for (int seg = 0; seg < pathTemp1.segments.size(); seg++) {
			for (int i = 0; i < 1 / step; i++) {
				sum += pathTemp1.segments[seg]->sDiff(i * step) * exp(k * pathTemp1.segments[seg]->CurvatureV2(i * step));
			}
		}
		
		//cout << "Sum: " << sum << endl;
		
		/*
		for (int border = 0; border < borders.size(); border++) {
			//cout << border << endl;
			totalAreaTemp1 += Area(1, border);
		}
		*/
		
		//cout << "mainFunctionForward Complete" << endl;
		//cout << "Total Area Forward : " << totalAreaTemp1 * multiplier << endl;
		
		cout << "MF Forward : " << sum / Vmax << endl;
				
		return (sum / Vmax); // + (totalAreaTemp1 * multiplier);
		//return 0;
	}
	
	void resetTempPaths () {
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param < 9; param++) {
				pathTemp0.segments[seg]->params[param] = path.segments[seg]->params[param];
				pathTemp1.segments[seg]->params[param] = path.segments[seg]->params[param];
				pathNext.segments[seg]->params[param] = path.segments[seg]->params[param];
			}
		}
		
		pathTemp0.produceLines();
		pathTemp1.produceLines();
		pathNext.produceLines();
	}
	
	void upgradePath () {
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param < 9; param++) {
				path.segments[seg]->params[param] = pathNext.segments[seg]->params[param];
			}
		}
		
		path.produceLines();
	}
	
	void upgradePathTo (BezierPiecewise& upgradeTo) {
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param < 9; param++) {
				path.segments[seg]->params[param] = upgradeTo.segments[seg]->params[param];
			}
		}
	}
	
	void newtonsParameterUpdateV2 (int iter) {
		const int iterations = iter;
		
		for (int i = 0; i < iterations; i++) {
			cout << "Iteration: " << i + 1 << endl;
			for (int seg = 0; seg < path.segments.size(); seg++) {
				cout << (float(seg) / float(path.segments.size())) * 100 << "%";
				if (seg == path.segments.size() - 1) {
					for (int param = 4; param <= 8; param++) {
						updateTempPaths(seg, param);
						
						//double topTemp = (((mainFunctionForward() + insideSegmentNewtonsForward()) - (mainFunctionBackward() + insideSegmentNewtonsBackward())) / (2 * step));
						//double bottomTemp = (((mainFunctionForward() + insideSegmentNewtonsForward()) - (2 * (mainFunction() + insideSegmentNewtons())) + (mainFunctionBackward() + insideSegmentNewtonsBackward())) / (step * step));
						
						//double topTemp = (((mainFunctionForward()) - (mainFunctionBackward())) / (2 * step));
						//double bottomTemp = (((mainFunctionForward()) - (2 * (mainFunction())) + (mainFunctionBackward())) / (step * step));
						//double topTemp = ((insideSegmentNewtonsForward() - insideSegmentNewtonsBackward()) / (2 * step));
						//double bottomTemp = (insideSegmentNewtonsForward() - (2 * insideSegmentNewtons()) + insideSegmentNewtonsBackward()) / (step * step);
						
						double topTemp = insideSegmentNewtons();
						double bottomTemp = ((insideSegmentNewtonsForward() - insideSegmentNewtonsBackward()) / (2 * step));
						//cout << "\n\n\n\nBottom" << bottomTemp << "\n\nTop" << topTemp << "\n\n\n\n";
						if (bottomTemp != 0) {
							path.segments[seg]->params[param] = path.segments[seg]->params[param] - (nSlowDown * (topTemp / bottomTemp));
							if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
							else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
						} else cout << ", skipping param " << param;
						
						resetTempPaths();
					}
				} else {
					for (int param = 4; param <= 8; param++) {
						updateTempPaths(seg, param);
						
						double topTemp = insideSegmentNewtons();
						double bottomTemp = ((insideSegmentNewtonsForward() - insideSegmentNewtonsBackward()) / (2 * step));
						//cout << "\n\n\n\nBottom" << bottomTemp << "\n\nTop" << topTemp << "\n\n\n\n";
						if (bottomTemp != 0) {
							path.segments[seg]->params[param] = path.segments[seg]->params[param] - (nSlowDown * (topTemp / bottomTemp));
							if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
							else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
						} else cout << ", skipping param " << param;
						
						resetTempPaths();
					}
				}
				cout << endl;
			}
		}
		
		cout << endl;
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param <= 8; param++) {
				cout << path.segments[seg]->params[param] << ", ";
			} cout << endl;
		}
	}
	
	void newtonSingleParameter (int iterations, int seg, int param) {
		resetTempPaths();
		const double minimumDifference = 0.0001;
		
		double nextGuess = 0;
		double previousGuess = 0;
		
		double previousConditionTemp = 0;
		
		if (seg == path.segments.size() - 1) {
			for (int iter = 0; iter < iterations; iter++) {
				updateTempPaths(seg, param);

				//double topTemp = (((mainFunctionForward() + insideSegmentNewtonsForward()) - (mainFunctionBackward() + insideSegmentNewtonsBackward())) / (2 * step));
				//double bottomTemp = (((mainFunctionForward() + insideSegmentNewtonsForward()) - (2 * (mainFunction() + insideSegmentNewtons())) + (mainFunctionBackward() + insideSegmentNewtonsBackward())) / (step * step));

				//double topTemp = (((mainFunctionForward()) - (mainFunctionBackward())) / (2 * step));
				//double bottomTemp = (((mainFunctionForward()) - (2 * (mainFunction())) + (mainFunctionBackward())) / (step * step));
				//double topTemp = ((insideSegmentNewtonsForward() - insideSegmentNewtonsBackward()) / (2 * step));
				//double bottomTemp = (insideSegmentNewtonsForward() - (2 * insideSegmentNewtons()) + insideSegmentNewtonsBackward()) / (step * step);
				
				/*
				double topTemp = insideSegmentNewtons();
				cout << "Current Condition ================================== " << topTemp << endl;
				double bottomTemp = ((insideSegmentNewtonsForward() - insideSegmentNewtonsBackward()) / (2 * step));
				*/
				
				/*
				double topTemp = insideSegmentNewtonsV3(path);
				cout << "Current Condition ================================== " << topTemp << endl;
				double forwardTemp = insideSegmentNewtonsV3(pathTemp1);
				double backwardTemp = insideSegmentNewtonsV3(pathTemp0);
				double bottomTemp = ((forwardTemp - backwardTemp) / (2 * step));
				*/
				
				//double topTemp = (fullEquationT(pathTemp1) - fullEquationT(pathTemp0)) / (2 * step);//ExternalAreaV2(path); //insideSegmentNewtonsV3(path)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(path)*pow(10, 1))));
				//cout << "Current Condition ================================== " << fullEquationT(path) << endl;
				//double bottomTemp = (fullEquationT(pathTemp1) - 2*fullEquationT(path) + fullEquationT(pathTemp0)) / (step * step);//(ExternalAreaV2(pathTemp1) - ExternalAreaV2(pathTemp0)) / (2 * step);//(((insideSegmentNewtonsV3(pathTemp1)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(pathTemp1)*pow(10, 1))))) -2*(insideSegmentNewtonsV3(path)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(path)*pow(10, 1))))) + (insideSegmentNewtonsV3(pathTemp0)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(pathTemp0)*pow(10, 1))))) )) / (step * step);
				
				double topTemp = fullEquationT(path);//ExternalAreaV2(path); //insideSegmentNewtonsV3(path)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(path)*pow(10, 1))));
				cout << "Current Condition ================================== " << fullEquationT(path) << endl;
				double bottomTemp = (fullEquationT(pathTemp1) - fullEquationT(pathTemp0)) / (2 * step);
				
				if (bottomTemp != 0) {
					previousGuess = path.segments[seg]->params[param];
					nextGuess = previousGuess - (nSlowDown * (topTemp / bottomTemp));
					pathNext.segments[seg]->params[param] = nextGuess;
					if (param == 4) pathNext.segments[0]->updateBX(pathNext.segments[seg]->params[4]);
					else if (param== 5) pathNext.segments[0]->updateBY(pathNext.segments[seg]->params[5]);
					else if (param == 6) {
						pathNext.segments[0]->params[0] = pathNext.segments[seg]->params[6];
						pathNext.segments[0]->updateBX(pathNext.segments[seg]->params[4]);
						pathNext.segments[0]->updateBY(pathNext.segments[seg]->params[5]);
					} else if (param == 7) {
						pathNext.segments[0]->params[1] = pathNext.segments[seg]->params[7];
						pathNext.segments[0]->updateBX(pathNext.segments[seg]->params[4]);
						pathNext.segments[0]->updateBY(pathNext.segments[seg]->params[5]);
					} else if (param == 8) pathNext.segments[seg]->updateCL(pathNext.segments[seg]->params[8]);
					
					//pathNext.produceLines();
					double nextTemp = fullEquationT(pathNext); //ExternalAreaV2(pathNext); //insideSegmentNewtonsV3(pathNext)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(pathNext)*pow(10, 1))));
					
					if (nextTemp != -1) {
						if (previousCondition > nextTemp) {
							upgradePath();
							/*
							path.segments[seg]->params[param] = nextGuess;
							if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
							else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
							*/
							previousCondition = (fullEquationT(pathTemp1) - fullEquationT(pathTemp0)) / (2 * step);
						} else {
							cout << "skipping" << endl;
							break;
						}
					} else {
						cout << "skipping" << endl;
						break;
					}
				}
				
				resetTempPaths();
				
				if (abs(nextGuess - previousGuess) < 2 * minimumDifference) break;
			}
		} else {
			for (int iter = 0; iter < iterations; iter++) {
				updateTempPaths(seg, param);

				double topTemp = fullEquationT(path);
				cout << "Current Condition ================================== " << fullEquationT(path) << endl;
				double bottomTemp = (fullEquationT(pathTemp1) - fullEquationT(pathTemp0)) / (2 * step);

				if (bottomTemp != 0) {
					previousGuess = path.segments[seg]->params[param];
					nextGuess = previousGuess - (nSlowDown * (topTemp / bottomTemp));
					pathNext.segments[seg]->params[param] = nextGuess;
					if (param == 4) pathNext.segments[seg + 1]->updateBX(pathNext.segments[seg]->params[4]);
					else if (param == 5) pathNext.segments[seg + 1]->updateBY(pathNext.segments[seg]->params[5]);
					else if (param == 6) {
						pathNext.segments[seg + 1]->params[0] = pathNext.segments[seg]->params[6];
						pathNext.segments[seg + 1]->updateBX(pathNext.segments[seg]->params[4]);
						pathNext.segments[seg + 1]->updateBY(pathNext.segments[seg]->params[5]);
					} else if (param == 7) {
						pathNext.segments[seg + 1]->params[1] = pathNext.segments[seg]->params[7];
						pathNext.segments[seg + 1]->updateBX(pathNext.segments[seg]->params[4]);
						pathNext.segments[seg + 1]->updateBY(pathNext.segments[seg]->params[5]);
					} else if (param == 8) pathNext.segments[seg]->updateCL(pathNext.segments[seg]->params[8]);
					
					pathNext.produceLines();
					double nextTemp = fullEquationT(pathNext); 
					
					if (nextTemp != -1) {
						if (previousCondition > nextTemp) {
							upgradePath();
							/*
							path.segments[seg]->params[param] = nextGuess;
							if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
							else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
							*/
							previousCondition = (fullEquationT(pathTemp1) - fullEquationT(pathTemp0)) / (2 * step);
						} else {
							cout << "skipping" << endl;
							break;
						}
					} else {
						cout << "skipping" << endl;
						break;
					}
				}

				resetTempPaths();
				
				if (abs(nextGuess - previousGuess) < 2 * minimumDifference) break;
			}
		}
	}
	
	void newtonsParameterUpdateV3 (int iter) {
		const int iterations = iter;
		
		double currentCondition = 0;
		double diffParams[5*path.segments.size()] = {};
		int maxDiffParam = 0;
		int derivedSeg = 0;
		int derivedParam = 0;
		
		for (int i = 0; i < iterations; i++) {
			cout << "Iteration: " << i + 1;
			currentCondition = insideSegmentNewtons();
			cout << ", Current Condition: " << currentCondition << endl;
			
			for (int seg = 0; seg < path.segments.size(); seg++) {
				for (int param = 4; param <= 8; param++) {
					cout << (float(seg * (param - 4)) / float(path.segments.size() + 5)) * 100 << "%";
					updateTempPaths(seg, param);
					diffParams[seg*5 + (param - 4)] = ((insideSegmentNewtonsForward() - insideSegmentNewtonsBackward()) / (2 * step));
					resetTempPaths();
				}
			}
				
			for (int param = 0; param < 5*path.segments.size(); param++) {
				if (abs(diffParams[param]) < abs(diffParams[maxDiffParam])) maxDiffParam = param;
			}
			
			if (diffParams[maxDiffParam] != 0) {
				derivedParam = maxDiffParam % 5;
				derivedSeg = (maxDiffParam - derivedParam) / path.segments.size();
				
				if (derivedSeg == path.segments.size() - 1) {
					path.segments[derivedSeg]->params[derivedParam + 4] = path.segments[derivedSeg]->params[derivedParam + 4] - (nSlowDown * (currentCondition / diffParams[maxDiffParam]));
					if (derivedParam + 4 == 4) path.segments[0]->updateBX(path.segments[derivedSeg]->params[4]);
					else if (derivedParam + 4 == 5) path.segments[0]->updateBY(path.segments[derivedSeg]->params[5]);
					else if (derivedParam + 4 == 6) path.segments[0]->params[0] = path.segments[derivedSeg]->params[6];
					else if (derivedParam + 4 == 7) path.segments[0]->params[1] = path.segments[derivedSeg]->params[7];
					else if (derivedParam + 4 == 8) path.segments[derivedSeg]->updateCL(path.segments[derivedSeg]->params[8]);
				} else {
					path.segments[derivedSeg]->params[derivedParam + 4] = path.segments[derivedSeg]->params[derivedParam + 4 ] - (nSlowDown * (currentCondition / diffParams[maxDiffParam]));
					if (derivedParam + 4 == 4) path.segments[derivedSeg + 1]->updateBX(path.segments[derivedSeg]->params[4]);
					else if (derivedParam + 4 == 5) path.segments[derivedSeg + 1]->updateBY(path.segments[derivedSeg]->params[5]);
					else if (derivedParam + 4 == 6) path.segments[derivedSeg + 1]->params[0] = path.segments[derivedSeg]->params[6];
					else if (derivedParam + 4 == 7) path.segments[derivedSeg + 1]->params[1] = path.segments[derivedSeg]->params[7];
					else if (derivedParam + 4 == 8) path.segments[derivedSeg]->updateCL(path.segments[derivedSeg]->params[8]);
				}
			}
		}
		
		cout << endl;
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param <= 8; param++) {
				cout << path.segments[seg]->params[param] << ", ";
			} cout << endl;
		}
	}
	
	void manualIteration () {
		const int iterations = 15;
		
		for (int i = 0; i < iterations; i++) {
			cout << "Iteration: " << i + 1 << endl;
			for (int seg = 0; seg < path.segments.size(); seg++) {
				cout << (float(seg) / float(path.segments.size())) * 100 << "%\n";
				if (seg == path.segments.size() - 1) {
					for (int param = 4; param <= 8; param++) {
						updateTempPathsManual(seg, param);
						
						double tempCentre = insideSegmentNewtons();
						double tempForward = insideSegmentNewtonsForward();
						double tempBackward = insideSegmentNewtonsBackward();
						
						if (tempForward < tempCentre) {
							cout << "Increasing" << endl;
							path.segments[seg]->params[param] += moveStep;
							if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
							else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
						} else if (tempBackward < tempCentre) {
							cout << "Reducing" << endl;
							path.segments[seg]->params[param] -= moveStep;
							if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
							else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
						} else cout << ", skipping param " << param;
						
						resetTempPaths();
					}
				} else {
					for (int param = 4; param <= 8; param++) {
						updateTempPathsManual(seg, param);
						
						double tempCentre = insideSegmentNewtons();
						double tempForward = insideSegmentNewtonsForward();
						double tempBackward = insideSegmentNewtonsBackward();
						
						if (tempForward < tempCentre) {
							cout << "Increasing" << endl;
							path.segments[seg]->params[param] += moveStep;
							if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
							else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
						} else if (tempBackward < tempCentre) {
							cout << "Reducing" << endl;
							path.segments[seg]->params[param] -= moveStep;
							if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
							else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
							else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
							else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
							else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
						} else cout << ", skipping param " << param;
						
						resetTempPaths();
					}
				}
				cout << endl;
			}
		}
		
		cout << endl;
		for (int seg = 0; seg < path.segments.size(); seg++) {
			for (int param = 0; param <= 8; param++) {
				cout << path.segments[seg]->params[param] << ", ";
			} cout << endl;
		}
	}
	
	void newtonSingleParamIterator (int iter, int paramIter) {
		resetTempPaths();
		
		// Set initial Previous Condition
		previousCondition = (fullEquationT(pathTemp1) - fullEquationT(pathTemp0)) / (2 * step); //insideSegmentNewtonsV3(path)*pow(10, 1) + ((mainFunctionRewrite(pathTemp1) - mainFunctionRewrite(pathTemp0)) / (2 * step * exp(insideSegmentNewtonsV3(path)*pow(10, 1))));
		
		for (int i = 0; i < iter; i++) {
			cout << "Iteration : " << i + 1 << endl;
			for (int seg = 0; seg < path.segments.size(); seg++) {
				cout << (float(seg) / float(path.segments.size())) * 100 << "%\n";
				for (int param = 4; param <= 8; param++) {
					newtonSingleParameter(paramIter, seg, param);
					//if (previousCondition <= 0) return;
				}
			}
		}
	
	}
	
	double lineIntersectionIterator (BezierPiecewise& currentPath) {
		
		double pointx = 0;
		double pointy = 0;
		double previousx = 0;
		double previousy = 0;
		
		double area = 0;
		int areaSegments = 10;
		
		for (int border = 0; border < borders.size(); border++) {
			vector <double> pathIntersections;
			vector <double> borderIntersections;
			pointx = Bezier(0, currentPath.segments[0]->params[0], currentPath.segments[0]->params[2], currentPath.segments[0]->params[4], currentPath.segments[0]->params[6]);
			pointy = Bezier(0, currentPath.segments[0]->params[1], currentPath.segments[0]->params[3], currentPath.segments[0]->params[5], currentPath.segments[0]->params[7]);
			
			for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
				for (int point = 0; point < 1 / interStep; point++) {
					previousx = pointx;
					previousy = pointy;
					pointx = Bezier(point * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
					pointy = Bezier(point * interStep, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
					
					double borderIntersection = lineIntersectionTest(border, previousx, previousy, pointx, pointy);
					
					if (borderIntersection != -1) {
						pathIntersections.push_back(point * interStep + pathSeg);
						borderIntersections.push_back(borderIntersection);
					}
				}
			}
			
			cout << pathIntersections.size() << " intersections with border " << border << endl;
			
			double pathIntersectionStart = pathIntersections[0];
			double pathIntersectionEnd = 0;
			double borderIntersectionStart = borderIntersections[0];
			double borderIntersectionEnd = 0;
			
			double lastBorderIntersection = borderIntersections[0];
			
			bool lookingForReEntry = true;
			
			if (pathIntersections.size() % 2 == 0) {
				for (int intersection = 0; intersection < pathIntersections.size(); intersection++) {
					if (lookingForReEntry) {
						if (borderIntersections[intersection] > lastBorderIntersection) {
							lookingForReEntry = false;
							
							pathIntersectionEnd = pathIntersections[intersection];
							borderIntersectionEnd = borderIntersections[intersection];
							
							double borderStep = 0;
							double pathStep = 0;
							double pathTempX = 0;
							double pathTempY = 0;
							double borderTempX = 0;
							double borderTempY = 0;
							int pathTempSeg = 0;
							int borderTempSeg = 0;
							double pathTempT = 0;
							double borderTempT = 0;

							borderStep = (borderIntersectionEnd - borderIntersectionStart) / areaSegments;
							pathStep = (pathIntersectionEnd - pathIntersectionStart) / areaSegments;

							for (int areaSeg = 1; areaSeg < areaSegments; areaSeg++) {
								pathTempSeg = floor(pathIntersectionStart + (areaSeg * pathStep));
								borderTempSeg = floor(borderIntersectionStart + (areaSeg * borderStep));

								pathTempT = (pathIntersectionStart + (areaSeg * pathStep)) - pathTempSeg;
								borderTempT = (borderIntersectionStart + (areaSeg * borderStep)) - borderTempSeg;

								pathTempX = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[0], currentPath.segments[pathTempSeg]->params[2], currentPath.segments[pathTempSeg]->params[4], currentPath.segments[pathTempSeg]->params[6]);
								pathTempY = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[1], currentPath.segments[pathTempSeg]->params[3], currentPath.segments[pathTempSeg]->params[5], currentPath.segments[pathTempSeg]->params[7]);

								borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
								borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);

								//outsidex.push_back(pathTempX);
								//outsidey.push_back(pathTempY);

								//outsidex.push_back(borderTempX);
								//outsidey.push_back(borderTempY);

								area += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2))*((pathStep + borderStep) / 2);
							}
						}
					} else {
						lookingForReEntry = true;
						
						pathIntersectionStart = pathIntersections[intersection];
						borderIntersectionStart = borderIntersections[intersection];
						
						lastBorderIntersection = borderIntersections[intersection];
					}
				}
			} else {
				cout << "Odd amount of intersections with border " << border << endl;
				return -1;
			}
		}
		
		return area;
	}
				
	void manualSingleParamIterator (int iter, int paramIter) {
		
		for (int i = 0; i < iter; i++) {
			cout << "Iteration : " << i + 1 << endl;
			for (int seg = 0; seg < path.segments.size(); seg++) {
				cout << (float(seg) / float(path.segments.size())) * 100 << "%\n";
				for (int param = 4; param <= 8; param++) {
					cout << "Param : " << param << endl;
					manualIterationSingle(paramIter, seg, param);
				}
			}
		}
	
	}
	
	void manualIterationSingle (int iterations, int seg, int param) {
		resetTempPaths();
		
		if (seg == path.segments.size() - 1) {
			for (int iter = 0; iter < iterations; iter++) {
				updateTempPathsManual(seg, param);
				
				double tempCentre = fullEquationT(path); //insideSegmentNewtons();
				cout << "Current Condition ================================== " << tempCentre << endl;
				double tempForward = fullEquationT(pathTemp1); //insideSegmentNewtonsForward();
				double tempBackward = fullEquationT(pathTemp0); //insideSegmentNewtonsBackward();
				
				if (tempForward < tempCentre) {
					cout << "Increasing" << endl;
					upgradePathTo(pathTemp1);
					/*
					path.segments[seg]->params[param] += moveStep;
					if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
					else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
					else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
					else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
					else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
					*/
				} else if (tempBackward < tempCentre) {
					cout << "Reducing" << endl;
					upgradePathTo(pathTemp0);
					/*
					path.segments[seg]->params[param] -= moveStep;
					if (param == 4) path.segments[0]->updateBX(path.segments[seg]->params[4]);
					else if (param== 5) path.segments[0]->updateBY(path.segments[seg]->params[5]);
					else if (param == 6) path.segments[0]->params[0] = path.segments[seg]->params[6];
					else if (param == 7) path.segments[0]->params[1] = path.segments[seg]->params[7];
					else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
					*/
				} else {
					cout << "skipping param" << param << endl;
					break;
				}
				
				resetTempPaths();
			}
		} else {
			for (int iter = 0; iter < iterations; iter++) {
				updateTempPathsManual(seg, param);

				double tempCentre = fullEquationT(path);
				cout << "Current Condition ================================== " << tempCentre << endl;
				double tempForward = fullEquationT(pathTemp1);
				double tempBackward = fullEquationT(pathTemp0);

				if (tempForward < tempCentre) {
					cout << "Increasing" << endl;
					upgradePathTo(pathTemp1);
					/*
					path.segments[seg]->params[param] += moveStep;
					if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
					else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
					else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
					else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
					else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
					*/
				} else if (tempBackward < tempCentre) {
					cout << "Reducing" << endl;
					upgradePathTo(pathTemp0);
					/*
					path.segments[seg]->params[param] -= moveStep;
					if (param == 4) path.segments[seg + 1]->updateBX(path.segments[seg]->params[4]);
					else if (param == 5) path.segments[seg + 1]->updateBY(path.segments[seg]->params[5]);
					else if (param == 6) path.segments[seg + 1]->params[0] = path.segments[seg]->params[6];
					else if (param == 7) path.segments[seg + 1]->params[1] = path.segments[seg]->params[7];
					else if (param == 8) path.segments[seg]->updateCL(path.segments[seg]->params[8]);
					*/
				} else {
					cout << "skipping param" << param << endl;
					break;
				}

				resetTempPaths();
			}
		}
	}
	
	void updateTempPaths (int segIndex, int paramIndex) {
		//cout << "Updating" << endl;
		if (segIndex < path.segments.size() - 1) {
			if (paramIndex == 4) {
				pathTemp0.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] - step;
				pathTemp0.segments[segIndex + 1]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp0.segments[segIndex + 1]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] + step;
				pathTemp1.segments[segIndex + 1]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp1.segments[segIndex + 1]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 5) {
				pathTemp0.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] - step;
				pathTemp0.segments[segIndex + 1]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp0.segments[segIndex + 1]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] + step;
				pathTemp1.segments[segIndex + 1]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp1.segments[segIndex + 1]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 6) {
				pathTemp0.segments[segIndex]->params[6] = path.segments[segIndex]->params[6] - step;
				pathTemp1.segments[segIndex]->params[6] = path.segments[segIndex]->params[6] + step;
				pathTemp0.segments[segIndex + 1]->params[0] = path.segments[segIndex]->params[6] - step;
				pathTemp1.segments[segIndex + 1]->params[0] = path.segments[segIndex]->params[6] + step;
				pathTemp0.segments[segIndex + 1]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp1.segments[segIndex + 1]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp0.segments[segIndex + 1]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[segIndex + 1]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 7) {
				pathTemp0.segments[segIndex]->params[7] = path.segments[segIndex]->params[7] - step;
				pathTemp1.segments[segIndex]->params[7] = path.segments[segIndex]->params[7] + step;
				pathTemp0.segments[segIndex + 1]->params[1] = path.segments[segIndex]->params[7] - step;
				pathTemp1.segments[segIndex + 1]->params[1] = path.segments[segIndex]->params[7] + step;
				pathTemp0.segments[segIndex + 1]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp1.segments[segIndex + 1]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp0.segments[segIndex + 1]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[segIndex + 1]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 8) {
				pathTemp0.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] - step);
				pathTemp1.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] + step);
			}
		}
		
		else if (segIndex == path.segments.size() - 1) {
			if (paramIndex == 4) {
				pathTemp0.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] - step;
				pathTemp0.segments[0]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp0.segments[0]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] + step;
				pathTemp1.segments[0]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp1.segments[0]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 5) {
				pathTemp0.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] - step;
				pathTemp0.segments[0]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp0.segments[0]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] + step;
				pathTemp1.segments[0]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp1.segments[0]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 6) {
				pathTemp0.segments[0]->params[0] = path.segments[segIndex]->params[6] - step;
				pathTemp1.segments[0]->params[0] = path.segments[segIndex]->params[5] + step;
				pathTemp0.segments[0]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp0.segments[0]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[0]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp1.segments[0]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 7) {
				pathTemp0.segments[0]->params[1] = path.segments[segIndex]->params[7] - step;
				pathTemp1.segments[0]->params[1] = path.segments[segIndex]->params[7] + step;
				pathTemp0.segments[0]->updateBX(pathTemp0.segments[segIndex]->params[4]);
				pathTemp0.segments[0]->updateBY(pathTemp0.segments[segIndex]->params[5]);
				pathTemp1.segments[0]->updateBX(pathTemp1.segments[segIndex]->params[4]);
				pathTemp1.segments[0]->updateBY(pathTemp1.segments[segIndex]->params[5]);
			} else if (paramIndex == 8) {
				pathTemp0.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] - step);
				pathTemp1.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] + step);
			}
		}
		
		pathTemp0.produceLines();
		pathTemp1.produceLines();
		//cout << "Finished Updating" << endl;
	}
	
	void updateTempPathsManual (int segIndex, int paramIndex) {
		//cout << "Updating" << endl;
		if (segIndex < path.segments.size() - 1) {
			if (paramIndex == 4) {
				pathTemp0.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] - moveStep;
				pathTemp0.segments[segIndex + 1]->updateBX(path.segments[segIndex]->params[4] - moveStep);
				pathTemp1.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] + moveStep;
				pathTemp1.segments[segIndex + 1]->updateBX(path.segments[segIndex]->params[4] + moveStep);
			} else if (paramIndex == 5) {
				pathTemp0.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] - moveStep;
				pathTemp0.segments[segIndex + 1]->updateBY(path.segments[segIndex]->params[5] - moveStep);
				pathTemp1.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] + moveStep;
				pathTemp1.segments[segIndex + 1]->updateBY(path.segments[segIndex]->params[5] + moveStep);
			} else if (paramIndex == 6) {
				pathTemp0.segments[segIndex + 1]->params[0] = path.segments[segIndex]->params[6] - moveStep;
				pathTemp1.segments[segIndex + 1]->params[0] = path.segments[segIndex]->params[5] + moveStep;
			} else if (paramIndex == 7) {
				pathTemp0.segments[segIndex + 1]->params[1] = path.segments[segIndex]->params[7] - moveStep;
				pathTemp1.segments[segIndex + 1]->params[1] = path.segments[segIndex]->params[7] + moveStep;
			} else if (paramIndex == 8) {
				pathTemp0.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] - moveStep);
				pathTemp1.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] + moveStep);
			}
		}
		
		else if (segIndex == path.segments.size() - 1) {
			if (paramIndex == 4) {
				pathTemp0.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] - moveStep;
				pathTemp0.segments[0]->updateBX(path.segments[segIndex]->params[4] - moveStep);
				pathTemp1.segments[segIndex]->params[4] = path.segments[segIndex]->params[4] + moveStep;
				pathTemp1.segments[0]->updateBX(path.segments[segIndex]->params[4] + moveStep);
			} else if (paramIndex == 5) {
				pathTemp0.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] - moveStep;
				pathTemp0.segments[0]->updateBY(path.segments[segIndex]->params[5] - moveStep);
				pathTemp1.segments[segIndex]->params[5] = path.segments[segIndex]->params[5] + moveStep;
				pathTemp1.segments[0]->updateBY(path.segments[segIndex]->params[5] + moveStep);
			} else if (paramIndex == 6) {
				pathTemp0.segments[0]->params[0] = path.segments[segIndex]->params[6] - moveStep;
				pathTemp1.segments[0]->params[0] = path.segments[segIndex]->params[5] + moveStep;
			} else if (paramIndex == 7) {
				pathTemp0.segments[0]->params[1] = path.segments[segIndex]->params[7] - moveStep;
				pathTemp1.segments[0]->params[1] = path.segments[segIndex]->params[7] + moveStep;
			} else if (paramIndex == 8) {
				pathTemp0.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] - moveStep);
				pathTemp1.segments[segIndex]->updateCL(path.segments[segIndex]->params[8] + moveStep);
			}
		}
		//cout << "Finished Updating" << endl;
	}
	
	double lineIntersectionTest (int border, double x0, double y0, double x1, double y1) {

		double lambda = 0;

		const int iterations = 5000;
		const double acceptableError = 0.0001;
		const double minimumDifference = 0.000000001;

		double previousGuess = 0;
		double nextGuess = 0;
		
		int guessSections = 100;
		
		/*
		double firstGuess = 0;
		double secondGuess = 0;
		double thirdGuess = 0;
		*/
		
		for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
			
			vector <double> guesses;
						
			for (int g = 0; g <= guessSections; g++) {
				nextGuess = double(g) / double(guessSections);

				for (int n = 0; n < iterations; n++) {
					previousGuess = nextGuess;
					nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
					if (abs(previousGuess - nextGuess) < minimumDifference) break;
				}

				// Store the guess
				if (nextGuess > 0 && nextGuess < 1)	guesses.push_back(nextGuess);
			}
			
			for (int guess = 0; guess < guesses.size(); guess++) {
				if (guesses[guess] > 0 && guesses[guess] < 1) {
					lambda = (Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
					if (lambda >= 0 && lambda < 1) {
						if (sqrt(pow(Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0))), 2) + pow(Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0))), 2)) < 0.01) {
						//if (abs(Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < lambda*acceptableError && abs(Bezier(guesses[guess], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < lambda*acceptableError) {
							outsidex.push_back(x0 + lambda*(x1 - x0));
							outsidey.push_back(y0 + lambda*(y1 - y0));
							return guesses[guess] + borderSeg;
						}
					}
				}
			}
		}
		
		return -1;
	}
			
			/*
			// Newtons from 0.0001
			nextGuess = 0.0001;
			for (int n = 0; n < iterations; n++) {
				previousGuess = nextGuess;
				nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
				if (abs(previousGuess - nextGuess) < minimumDifference) break;
			}

			// Store first guess
			firstGuess = nextGuess;

			// Newtons from 0.5
			nextGuess = 0.5;
			for (int n = 0; n < iterations; n++) {
				previousGuess = nextGuess;
				nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
				if (abs(previousGuess - nextGuess) < minimumDifference) break;
			}

			// Store second guess
			secondGuess = nextGuess;

			// Newtons from 0.9999
			nextGuess = 0.9999;
			for (int n = 0; n < iterations; n++) {
				previousGuess = nextGuess;
				nextGuess = previousGuess - (((2 * step) * Lambda(Bezier(previousGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)) / (Lambda(Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess + step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1) - Lambda(Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]), Bezier(previousGuess - step, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]), x0, x1, y0, y1)));
				if (abs(previousGuess - nextGuess) < minimumDifference) break;
			}

			// Store third guess
			thirdGuess = nextGuess;	


			// Test different Lambdas		
			if (firstGuess >= 0 && firstGuess < 1) {
				lambda = (Bezier(firstGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
				if (lambda >= 0 - acceptableError  && lambda < 1 + acceptableError) {
					if (abs(Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
						outsidex.push_back(x0 + lambda*(x1 - x0));
						outsidey.push_back(y0 + lambda*(y1 - y0));
						return true;
					}
				}
			}

			if (secondGuess >= 0 && secondGuess < 1) {
				lambda = (Bezier(secondGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
				if (lambda >= 0 - acceptableError  && lambda < 1 + acceptableError) {
					if (abs(Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
						outsidex.push_back(x0 + lambda*(x1 - x0));
						outsidey.push_back(y0 + lambda*(y1 - y0));
						return true;
					}
				}
			}

			if (thirdGuess >= 0 && thirdGuess < 1) {
				lambda = (Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - y0) / (y1 - y0);
				if (lambda >= 0 - acceptableError  && lambda < 1 + acceptableError) {
					if (abs(Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]) - (y0 + (lambda*(y1 - y0)))) < acceptableError * lambda && abs(Bezier(thirdGuess, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]) - (x0 + (lambda*(x1 - x0)))) < acceptableError) {
						outsidex.push_back(x0 + lambda*(x1 - x0));
						outsidey.push_back(y0 + lambda*(y1 - y0));
						return true;
					}
				}
			}
		}

		return false;
	}
	*/
		
	
	long double dualLineIntersectionTest (BezierPiecewise& currentPath) {
		//BezierPiecewise& currentPath = path;
		
		long double px0 = 0;
		long double py0 = 0;
		long double px1 = Bezier(0, currentPath.segments[0]->params[0], currentPath.segments[0]->params[2], currentPath.segments[0]->params[4], currentPath.segments[0]->params[6]);
		long double py1 = Bezier(0, currentPath.segments[0]->params[1], currentPath.segments[0]->params[3], currentPath.segments[0]->params[5], currentPath.segments[0]->params[7]);
		long long pSlope = 0;
		bool pUpright = false;
		long double pLambda = 0;
		
		
		long double bx0 = 0;
		long double by0 = 0;
		long double bx1 = 0;
		long double by1 = 0;
		long long bSlope = 0;
		bool bUpright = false;
		long double bLambda = 0;
		int bpoint = 1;
		
		long double area = 0;
		int areaSegments = 20;
		long double length = 0;
		double minimumDifference = 0.01;
		
		for (int border = 0; border < borders.size(); border++) {

			vector <long double> pathGuesses;
			vector <long double> borderGuesses;
			
			vector <long double> pathIntersections;
			vector <long double> borderIntersections;
			
			bx0 = Bezier(0, borders[border]->segments[0]->params[0], borders[border]->segments[0]->params[2], borders[border]->segments[0]->params[4], borders[border]->segments[0]->params[6]);
			by1 = Bezier(0, borders[border]->segments[0]->params[1], borders[border]->segments[0]->params[3], borders[border]->segments[0]->params[5], borders[border]->segments[0]->params[7]);

			for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int pathPoint = 1; pathPoint <= int(1 / interStepPath); pathPoint++) {
						px0 = px1;
						py0 = py1;

						px1 = Bezier(pathPoint * interStepPath, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
						py1 = Bezier(pathPoint * interStepPath, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
						
						if (px1 - px0 == 0) pUpright = true;
						else {
							pSlope = ceil(pow(10, 8) * ((py1 - py0) / (px1 - px0)));
							pUpright = false;
						}

						for (int borderPoint = 1; borderPoint <= int(1 / interStep); borderPoint++) {
							bx0 = bx1;
							by0 = by1;
							
							bx1 = Bezier(borderPoint * interStep, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]);
							by1 = Bezier(borderPoint * interStep, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]);
							
							//bpoint++;
							
							if (bx1 - bx0 == 0) bUpright = true;
							else {
								bSlope = ceil(pow(10, 8) * ((by1 - by0) / (bx1 - bx0)));
								bUpright = false;
							}
							
							// Start Comparing Lines
							if ((pUpright == true && bUpright == true) || pSlope == bSlope) break;
							
							else {
								//long double bLambda1 = ((px0 - bx0)*(py1 - py0) + (by0 - py0)*(px1 - px0)) / ((bx1 - bx0)*(py1 - py0) - (by1 - by0)*(px1 - px0));
								//long double pLambda1 = ((by1 - py0) + (((px0 - bx0)*(py1 - py0) + (by0 - py0)*(px1 - px0)) / ((bx1 - bx0)*(py1 - py0) - (by1 - by0)*(px1 - px0)))*(by1 - by0)) / (py1 - py0);
								//cout << ((px0 - px1)*(by0 - by1) - (py0 - py1)*(bx0 - bx1)) << endl;
								
								pLambda = ((px0 - bx0)*(by0 - by1) - (py0 - by0)*(bx0 - bx1)) / ((px0 - px1)*(by0 - by1) - (py0 - py1)*(bx0 - bx1));
								bLambda = ((px0 - bx0)*(py0 - py1) - (py0 - by0)*(px0 - px1)) / ((px0 - px1)*(by0 - by1) - (py0 - py1)*(bx0 - bx1));
								
							}
							
							if (pLambda > 0 && pLambda < 1 && bLambda > 0 && bLambda < 1) {
								/*
								//if (sqrt(pow((bx0 + bLambda*(bx1 - bx0)) - (px0 + pLambda*(px1 - px0)), 2) + pow((by0 + bLambda*(by1 - by0)) - (py0 + pLambda*(py1 - py0)), 2)) < 2) {
								if (
								long double borderTempTime = (borderPoint*interStep) - interStep*(1 - bLambda);
								long double pathTempTime = (pathPoint*interStepPath) - interStepPath*(1 - pLambda);

								long double borderTempX = Bezier(borderTempTime, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]);
								long double borderTempY = Bezier(borderTempTime, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]);

								long double pathTempX = Bezier(pathTempTime, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
								long double pathTempY = Bezier(pathTempTime, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
								
								//if (sqrt(pow(pathTempX - (px0 + pLambda*(px1 - px0)), 2) + pow(pathTempY - (py0 + pLambda*(py1 - py0)), 2)) < 1 && sqrt(pow(borderTempX - (bx0 + bLambda*(bx1 - bx0)), 2) + pow(borderTempY - (by0 + bLambda*(by1 - by0)), 2)) < 1) {
								
									outsidex.push_back(bx0 + bLambda*(bx1 - bx0));
									outsidey.push_back(by0 + bLambda*(by1 - by0));
									outsidex.push_back(px0 + pLambda*(px1 - px0));
									outsidey.push_back(py0 + pLambda*(py1 - py0));
									cout << "InterSection" << endl;
								*/
								long double pathIntersectionX = px0 + pLambda*(px1 - px0);
								long double pathIntersectionY = py0 + pLambda*(py1 - py0);
								long double borderIntersectionX = bx0 + bLambda*(bx1 - bx0);
								long double borderIntersectionY = by0 + bLambda*(by1 - by0);
								
								long double pathTempTime = (pathPoint - (1.0 - pLambda))*interStepPath;
								long double pathTempX = Bezier(pathTempTime, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
								long double pathTempY = Bezier(pathTempTime, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
								long double borderTempTime = (borderPoint - (1.0 - bLambda))*interStepPath;
								long double borderTempX = Bezier(borderTempTime, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]);
								long double borderTempY = Bezier(borderTempTime, borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]);

								double distanceToCurve = sqrt(pow(pathIntersectionX - pathTempX, 2) + pow(pathIntersectionY - pathTempY, 2));
								
								bool pointInside = true;

								if (abs(pathIntersectionX - borderIntersectionX) > minimumDifference || abs(pathIntersectionY - borderIntersectionY) > minimumDifference) pointInside = false;
								
								if (abs(pathTempX - borderTempX) > minimumDifference || abs(pathTempY - borderTempY) > minimumDifference) pointInside = false;

								if (px0 < px1){
									if (!(px0 <= pathIntersectionX) || !(pathIntersectionX <= px1)) pointInside = false;
								} else {
									if (!(px0 >= pathIntersectionX) || !(pathIntersectionX >= px1)) pointInside = false;
								}

								if (py0 < py1) {
									if (!(py0 <= pathIntersectionY) || !(pathIntersectionY <= py1)) pointInside = false;
								} else {
									if (!(py0 >= pathIntersectionY) || !(pathIntersectionY >= py1)) pointInside = false;
								}

								if (bx0 < bx1){
									if (!(bx0 <= borderIntersectionX) || !(borderIntersectionX <= bx1)) pointInside = false;
								} else {
									if (!(bx0 >= borderIntersectionX) || !(borderIntersectionX >= bx1)) pointInside = false;
								}

								if (by0 < by1) {
									if (!(by0 <= borderIntersectionY) || !(borderIntersectionY <= by1)) pointInside = false;
								} else {
									if (!(by0 >= borderIntersectionY) || !(borderIntersectionY >= by1)) pointInside = false;
								}


								if (px0 < px1){
									if (!(px0 <= pathTempX) || !(pathTempX <= px1)) pointInside = false;
								} else {
									if (!(px0 >= pathTempX) || !(pathTempX >= px1)) pointInside = false;
								}

								if (py0 < py1) {
									if (!(py0 <= pathTempY) || !(pathTempY <= py1)) pointInside = false;
								} else {
									if (!(py0 >= pathTempY) || !(pathTempY >= py1)) pointInside = false;
								}

								if (bx0 < bx1){
									if (!(bx0 <= borderTempX) || !(borderTempX <= bx1)) pointInside = false;
								} else {
									if (!(bx0 >= borderTempX) || !(borderTempX >= bx1)) pointInside = false;
								}

								if (by0 < by1) {
									if (!(by0 <= borderTempY) || !(borderTempY <= by1)) pointInside = false;
								} else {
									if (!(by0 >= borderTempY) || !(borderTempY >= by1)) pointInside = false;
								}

								if (distanceToCurve > minimumDifference) pointInside = false;
								//} else pointInside = false;
								
								if (pointInside) {
									outsidex.push_back(pathIntersectionX);
									outsidey.push_back(pathIntersectionY);
									//outsidex.push_back(bx1);
									//outsidey.push_back(by1);
									//outsidex.push_back(px0 + pLambda*(px1 - px0));
									//outsidey.push_back(py0 + pLambda*(py1 - py0));
									cout << "InterSection" << endl;
									
									pathGuesses.push_back(((double(pathPoint) - (1 - pLambda))*interStepPath) + pathSeg);
									borderGuesses.push_back(((double(borderPoint) - (1 - bLambda))*interStep) + borderSeg);
								}
							}
						}
					}
				}
			}
			
			bool duplicate = false;
			long double minimumDifference = 0.00001;
			
			// Filter for Duplicates
			for (int i = 0; i < pathGuesses.size(); i++) {
				duplicate = false;

				for (int j = 0; j < pathIntersections.size(); j++) {
					if (abs(pathIntersections[j] - pathGuesses[i]) < minimumDifference || abs(borderIntersections[j] - borderGuesses[i]) < minimumDifference) {
						duplicate = true;
						break;
					}
				}

				if (!duplicate) {
					pathIntersections.push_back(pathGuesses[i]);
					borderIntersections.push_back(borderGuesses[i]);

					//double borderTempX = Bezier(borderGuesses[i], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[6]);
					//double borderTempY = Bezier(borderGuesses[i], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[7]);


					//outsidex.push_back(borderGuesses[i]);
					//outsidey.push_back(borderGuesses[i]);
				}
			}
			
			if (pathIntersections.size() % 2 == 0 && pathIntersections.size() > 0) {
				long double pathIntersectionStart = pathIntersections[0];
				long double pathIntersectionEnd = 0;
				long double borderIntersectionStart = borderIntersections[0];
				long double borderIntersectionEnd = 0;

				long double lastBorderIntersection = borderIntersections[0];

				bool lookingForReEntry = true;

				long double borderStep = 0;
				long double pathStep = 0;
				long double pathTempX = 0;
				long double pathTempY = 0;
				long double borderTempX = 0;
				long double borderTempY = 0;
				int pathTempSeg = 0;
				int borderTempSeg = 0;
				long double pathTempT = 0;
				long double borderTempT = 0;
				
				for (int intersection = 1; intersection < pathIntersections.size(); intersection++) {
					if (lookingForReEntry) {
						if (borderIntersections[intersection] > lastBorderIntersection) {
							lookingForReEntry = false;
							
							pathIntersectionEnd = pathIntersections[intersection];
							borderIntersectionEnd = borderIntersections[intersection];

							borderStep = (borderIntersectionEnd - borderIntersectionStart) / areaSegments;
							pathStep = (pathIntersectionEnd - pathIntersectionStart) / areaSegments;
							
							for (int areaSeg = 1; areaSeg < areaSegments; areaSeg++) {
								
								pathTempSeg = floor(pathIntersectionStart + (areaSeg * borderStep));
								borderTempSeg = floor(borderIntersectionStart + (areaSeg * borderStep));
								
								pathTempT = (pathIntersectionStart + (areaSeg * pathStep)) - pathTempSeg;
								borderTempT = (borderIntersectionStart + (areaSeg * borderStep)) - borderTempSeg;
								pathTempX = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[0], currentPath.segments[pathTempSeg]->params[2], currentPath.segments[pathTempSeg]->params[4], currentPath.segments[pathTempSeg]->params[6]);
								pathTempY = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[1], currentPath.segments[pathTempSeg]->params[3], currentPath.segments[pathTempSeg]->params[5], currentPath.segments[pathTempSeg]->params[7]);
								

								borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
								borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);

								//outsidex.push_back(pathTempX);
								//outsidey.push_back(pathTempY);

								//outsidex.push_back(borderTempX);
								//outsidey.push_back(borderTempY);

								area += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2))*((abs(pathStep) + abs(borderStep)) / 2);
							}
						}
					} else {
						lookingForReEntry = true;
						pathIntersectionStart = pathIntersections[intersection];
						borderIntersectionStart = borderIntersections[intersection];
						lastBorderIntersection = borderIntersections[intersection];
					}
				}
			} else {
				cout << "Odd amount of intersections with border " << border << endl;
			}
		}
		
		return area;
	}
	
	double mainFunctionRewrite (BezierPiecewise& currentPath) {
		//BezierPiecewise& currentPath = path;
		
		double sum = 0;
		
		for (int seg = 0; seg < currentPath.segments.size(); seg++) {
			for (int i = 0; i < 1 / step; i++) {
				sum += currentPath.segments[seg]->sDiff(i * step) * exp(k * currentPath.segments[seg]->CurvatureV2(i * step));
			}
		}
		
		return sum / Vmax;
	}			
	
	double insideSegmentNewtonsV3 (BezierPiecewise& currentPath) {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		double distance = 0;
		bool previousInside = true;
		double previousx = 0;
		double previousy = 0;
		
		for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				pointx = Bezier(point * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
				vector<double> distance0;
				vector<double> distance1;
				insidePath = false;
				
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					//vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							} else if (nextGuess > 0 && nextGuess < 1) {
								distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
								distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
							}
						}
						*/
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						} else if (nextGuess > 0 && nextGuess < 1) {
							distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
							distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
						}
					}
				}
				
				if (!insidePath) {
					outsideTest++;
					int side = 0;
					double minDist = 9999999999999;
					
					for (int dist = 0; dist < distance0.size(); dist++) {
						if (distance0[dist] < minDist) minDist = distance0[dist];
						if (distance1[dist] < minDist) minDist = distance1[dist];
					}
					
					distance += minDist;
					
					curvature += abs(BezierTangentAngle((point + 1) * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]) - pointx, 2) + pow(Bezier((point + 1) * interStep, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]) - pointy, 2));
					
					if (previousInside) distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = false;
					
				} else if (!previousInside) {
					distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = true;
				}
				previousx = pointx;
				previousy = pointy;
			}
		}
		
		double outside = (currentPath.segments.size() / interStep) - pointsOutside;
		
		//double areaIntersection = lineIntersectionIterator(currentPath);
		
		//cout << "Points Outside : " << outside << endl;
		//cout << "Outside Real : Outside Test = " << double(outside) / double(outsideTest) << endl;
		//cout << "Total : Outside Test + Inside = " << (path.segments.size() / interStep) / double(outsideTest + pointsOutside) << endl;
		cout << "Inside : Total = " << double(pointsOutside) / double((currentPath.segments.size() / interStep)) << endl;
		cout << "Outside : Total = " << double(outsideTest) / double((currentPath.segments.size() / interStep)) << endl;
		
		//if (areaIntersection != -1){
		//	if (outside < 0) return 0;
		//	else return areaIntersectionMultiplier*areaIntersection+ outsideMultiplier*pow((double(outsideTest) / double((currentPath.segments.size() / interStep))), 4) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 2) + curveMultiplier*pow(curvature, 2); //+ outsideMultiplier*outside + distanceMultiplier*distance;
		//} else {
			if (outside < 0) return 0;
			else return outsideMultiplier*pow((double(outsideTest) / double((currentPath.segments.size() / interStep))), 1) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 1) + curveMultiplier*pow(curvature, 2); //+ outsideMultiplier*outside + distanceMultiplier*distance;
		//}
	}
	
	double insideSegmentNewtonsOriginal () {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.00001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool outsidePath = true;
		int outsideTest = 0;
		
		for (int pathSeg = 0; pathSeg < path.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				pointx = Bezier(point * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							}
						}
					}
				}
				
				if (outsidePath) {
					outsideTest++;
					curvature += abs(BezierTangentAngle((point + 1) * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[6], path.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[6], path.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]) - pointx, 2) + pow(Bezier((point + 1) * interStep, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]) - pointy, 2));
					
					//cout << "Segment " << pathSeg << ", Point " << point << " of " << 1 / interStep << " in segment is outside at time = " << point * interStep << endl;
				}
				outsidePath = true;
			}
		}
		
		double outside = (path.segments.size() / interStep) - pointsOutside;
		
		//cout << "Points Outside : " << outside << endl;
		//cout << "Outside Real : Outside Test = " << double(outside) / double(outsideTest) << endl;
		//cout << "Total : Outside Test + Inside = " << (path.segments.size() / interStep) / double(outsideTest + pointsOutside) << endl;
		cout << "Inside : Outside = " << double(pointsOutside) / double(outsideTest) << endl;
		
		if (outside < 0) return 0;
		else return lengthMultiplier*length + curveMultiplier*curvature + outsideMultiplier*outside;
	}
	
	double insideSegmentNewtons () {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		double distance = 0;
		bool previousInside = true;
		double previousx = 0;
		double previousy = 0;
		
		for (int pathSeg = 0; pathSeg < path.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				pointx = Bezier(point * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
				vector<double> distance0;
				vector<double> distance1;
				insidePath = false;
				
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					//vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							} else if (nextGuess > 0 && nextGuess < 1) {
								distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
								distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
							}
						}
						*/
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						} else if (nextGuess > 0 && nextGuess < 1) {
							distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
							distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
						}
					}
				}
				
				if (!insidePath) {
					outsideTest++;
					int side = 0;
					double minDist = 9999999999999;
					
					for (int dist = 0; dist < distance0.size(); dist++) {
						if (distance0[dist] < minDist) minDist = distance0[dist];
						if (distance1[dist] < minDist) minDist = distance1[dist];
					}
					
					//distance += minDist;
					
					curvature += abs(BezierTangentAngle((point + 1) * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[6], path.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[6], path.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]) - pointx, 2) + pow(Bezier((point + 1) * interStep, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]) - pointy, 2));
					
					if (previousInside) distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = false;
					
				} else if (!previousInside) {
					distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = true;
				}
				previousx = pointx;
				previousy = pointy;
			}
		}
		
		double outside = (path.segments.size() / interStep) - pointsOutside;
		
		//cout << "Points Outside : " << outside << endl;
		//cout << "Outside Real : Outside Test = " << double(outside) / double(outsideTest) << endl;
		//cout << "Total : Outside Test + Inside = " << (path.segments.size() / interStep) / double(outsideTest + pointsOutside) << endl;
		cout << "Inside : Total = " << double(pointsOutside) / double((path.segments.size() / interStep)) << endl;
		cout << "Outside : Total = " << double(outsideTest) / double((path.segments.size() / interStep)) << endl;
		
		if (outside < 0) return 0;
		else return outsideMultiplier*(double(outsideTest) / double((path.segments.size() / interStep))) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 1) + curveMultiplier*pow(curvature, 2); //+ outsideMultiplier*outside + distanceMultiplier*distance;
	}
	
	double insideSegmentNewtonsNext () {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		double distance = 0;
		bool previousInside = true;
		double previousx = 0;
		double previousy = 0;
		
		for (int pathSeg = 0; pathSeg < path.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				pointx = Bezier(point * interStep, pathNext.segments[pathSeg]->params[0], pathNext.segments[pathSeg]->params[2], pathNext.segments[pathSeg]->params[4], pathNext.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, pathNext.segments[pathSeg]->params[1], pathNext.segments[pathSeg]->params[3], pathNext.segments[pathSeg]->params[5], pathNext.segments[pathSeg]->params[7]);
				vector<double> distance0;
				vector<double> distance1;
				insidePath = false;
				
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					//vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							} else if (nextGuess > 0 && nextGuess < 1) {
								distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
								distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
							}
						}
						*/
						
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						} else if (nextGuess > 0 && nextGuess < 1) {
							distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
							distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
						}
					}
				}
				
				if (!insidePath) {
					outsideTest++;
					int side = 0;
					double minDist = 9999999999999;
					
					for (int dist = 0; dist < distance0.size(); dist++) {
						if (distance0[dist] < minDist) minDist = distance0[dist];
						if (distance1[dist] < minDist) minDist = distance1[dist];
					}
					
					//distance += minDist;
					
					curvature += abs(BezierTangentAngle((point + 1) * interStep, pathNext.segments[pathSeg]->params[0], pathNext.segments[pathSeg]->params[1], pathNext.segments[pathSeg]->params[2], pathNext.segments[pathSeg]->params[3], pathNext.segments[pathSeg]->params[4], pathNext.segments[pathSeg]->params[5], pathNext.segments[pathSeg]->params[6], pathNext.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, pathNext.segments[pathSeg]->params[0], pathNext.segments[pathSeg]->params[1], pathNext.segments[pathSeg]->params[2], pathNext.segments[pathSeg]->params[3], pathNext.segments[pathSeg]->params[4], pathNext.segments[pathSeg]->params[5], pathNext.segments[pathSeg]->params[6], pathNext.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, pathNext.segments[pathSeg]->params[0], pathNext.segments[pathSeg]->params[2], pathNext.segments[pathSeg]->params[4], pathNext.segments[pathSeg]->params[6]) - pointx, 2) + pow(Bezier((point + 1) * interStep, pathNext.segments[pathSeg]->params[1], pathNext.segments[pathSeg]->params[3], pathNext.segments[pathSeg]->params[5], pathNext.segments[pathSeg]->params[7]) - pointy, 2));
					
					if (previousInside) distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = false;
				} else if (!previousInside) {
					distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = true;
				}
				previousx = pointx;
				previousy = pointy;
			}
		}
		
		double outside = (pathNext.segments.size() / interStep) - pointsOutside;
		
		//cout << "Points Outside : " << outside << endl;
		//cout << "Outside Real : Outside Test = " << double(outside) / double(outsideTest) << endl;
		//cout << "Total : Outside Test + Inside = " << (path.segments.size() / interStep) / double(outsideTest + pointsOutside) << endl;
		cout << "Inside : Total = " << double(pointsOutside) / double((path.segments.size() / interStep)) << endl;
		cout << "Outside : Total = " << double(outsideTest) / double((path.segments.size() / interStep)) << endl;
		
		if (outside < 0) return 0;
		else return outsideMultiplier*(double(outsideTest) / double((pathNext.segments.size() / interStep))) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 1) + curveMultiplier*pow(curvature, 2); //+ outsideMultiplier*outside + distanceMultiplier*distance;
	}
	
	double insideSegmentNewtonsForward () {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		double distance = 0;
		bool previousInside = true;
		double previousx = 0;
		double previousy = 0;

		for (int pathSeg = 0; pathSeg < pathTemp1.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				pointx = Bezier(point * interStep, pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[7]);
				vector<double> distance0;
				vector<double> distance1;
				insidePath = false;
				
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					//vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							} else if (nextGuess > 0 && nextGuess < 1) {
								distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
								distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
							}
						}
						*/
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						} else if (nextGuess > 0 && nextGuess < 1) {
							distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
							distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
						}
					}
				}
				
				if (!insidePath) {
					outsideTest++;
					
					int side = 0;
					double minDist = 9999999999999;
					
					for (int dist = 0; dist < distance0.size(); dist++) {
						if (distance0[dist] < minDist) minDist = distance0[dist];
						if (distance1[dist] < minDist) minDist = distance1[dist];
					}
					
					//distance += minDist;
					
					curvature += abs(BezierTangentAngle((point + 1) * interStep, pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[6], pathTemp1.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[6], pathTemp1.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[6]) - Bezier(point * interStep, pathTemp1.segments[pathSeg]->params[0], pathTemp1.segments[pathSeg]->params[2], pathTemp1.segments[pathSeg]->params[4], pathTemp1.segments[pathSeg]->params[6]), 2) + pow(Bezier((point + 1) * interStep, pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[7]) - Bezier(point * interStep, pathTemp1.segments[pathSeg]->params[1], pathTemp1.segments[pathSeg]->params[3], pathTemp1.segments[pathSeg]->params[5], pathTemp1.segments[pathSeg]->params[7]), 2));
					
					if (previousInside) distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = false;
				} else if (!previousInside) {
					distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = true;
				}
				previousx = pointx;
				previousy = pointy;
			}
		}
		
		double outside = (pathTemp1.segments.size() / interStep) - pointsOutside;
		//cout << "Points Outside Forward : " << outside << endl;
		//cout << "Forward Outside Real : " << outside << ", Outside Test : " << outsideTest << endl;
		
		if (outside < 0) return 0;
		else return outsideMultiplier*(double(outsideTest) / double((path.segments.size() / interStep))) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 1) + curveMultiplier*pow(curvature, 2); //lengthMultiplier*length + curveMultiplier*curvature + outsideMultiplier*outside + distanceMultiplier*distance;
	}
	
	double insideSegmentNewtonsBackward () {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		double distance = 0;
		bool previousInside = true;
		double previousx = 0;
		double previousy = 0;

		for (int pathSeg = 0; pathSeg < pathTemp0.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				pointx = Bezier(point * interStep, pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[7]);
				vector<double> distance0;
				vector<double> distance1;
				insidePath = false;
				
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							} else if (nextGuess > 0 && nextGuess < 1) {
								distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
								distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
							}
						}
						*/
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						} else if (nextGuess > 0 && nextGuess < 1) {
							distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
							distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
						}
					}
				}
				
				if (!insidePath) {
					outsideTest++;
					int side = 0;
					double minDist = 9999999999999;
					
					for (int dist = 0; dist < distance0.size(); dist++) {
						if (distance0[dist] < minDist) minDist = distance0[dist];
						if (distance1[dist] < minDist) minDist = distance1[dist];
					}
					
					//distance += minDist;
					
					curvature += abs(BezierTangentAngle((point + 1) * interStep, pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[6], pathTemp0.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[6], pathTemp0.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[6]) - Bezier(point * interStep, pathTemp0.segments[pathSeg]->params[0], pathTemp0.segments[pathSeg]->params[2], pathTemp0.segments[pathSeg]->params[4], pathTemp0.segments[pathSeg]->params[6]), 2) + pow(Bezier((point + 1) * interStep, pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[7]) - Bezier(point * interStep, pathTemp0.segments[pathSeg]->params[1], pathTemp0.segments[pathSeg]->params[3], pathTemp0.segments[pathSeg]->params[5], pathTemp0.segments[pathSeg]->params[7]), 2));
					
					if (previousInside) distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = false;
				} else if (!previousInside) {
					distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
					previousInside = true;
				}
				previousx = pointx;
				previousy = pointy;
			}
		}
		double outside = (pathTemp0.segments.size() / interStep) - pointsOutside;
		//cout << "Points Outside Backward : " << outside << endl;
		
		//cout << "Backward Outside Real : " << outside << ", Outside Test : " << outsideTest << endl;
		
		if (outside < 0) return 0;
		else return outsideMultiplier*(double(outsideTest) / double((path.segments.size() / interStep))) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 1) + curveMultiplier*pow(curvature, 2); //lengthMultiplier*length + curveMultiplier*curvature + outsideMultiplier*outside + distanceMultiplier*distance;
	}
	
	void insideSegmentNewtonsRecord () {
		const int iterations = 1000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = 0;
		double pointy = 0;
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		
		for (int pathSeg = 0; pathSeg < path.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				insidePath = false;
				pointx = Bezier(point * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]);
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					//vector<double> previousGuesses;
					for (int i = 1; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 0 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								insidePath = true;
								break;
							}
						}
						*/
						
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						}
					}
					
					if (insidePath) break;
				}
				
				if (!insidePath) {
					outsideTest++;
					curvature += abs(BezierTangentAngle((point + 1) * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[6], path.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[6], path.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, path.segments[pathSeg]->params[0], path.segments[pathSeg]->params[2], path.segments[pathSeg]->params[4], path.segments[pathSeg]->params[6]) - pointx, 2) + pow(Bezier((point + 1) * interStep, path.segments[pathSeg]->params[1], path.segments[pathSeg]->params[3], path.segments[pathSeg]->params[5], path.segments[pathSeg]->params[7]) - pointy, 2));

					outsidex.push_back(pointx);
					outsidey.push_back(pointy);
				}
			}
		}
		
		double outside = (path.segments.size() / interStep) - pointsOutside;
		
		//cout << "Points Outside : " << outside << endl;
		//cout << "Outside Real : Outside Test = " << double(outside) / double(outsideTest) << endl;
		//cout << "Total : Outside Test + Inside = " << (path.segments.size() / interStep) / double(outsideTest + pointsOutside) << endl;
		cout << "Inside : Total = " << double(pointsOutside) / double(path.segments.size() / interStep) << endl;
		cout << "Inside + Outside : Total = " << double(pointsOutside + outsideTest) / double(path.segments.size() / interStep) << endl;
	}
	
	double insideSegmentNewtonsV2 (BezierPiecewise& currentPath) {
		const int iterations = 2000;
		const double acceptableError = 0.05;
		const double minimumDifference = 0.000001;
		double Fx = 0;
		double fx = 0;
		int pointsOutside = 0;
		double pointx = Bezier(0, currentPath.segments[0]->params[0], currentPath.segments[0]->params[2], currentPath.segments[0]->params[4], currentPath.segments[0]->params[6]);
		double pointy = Bezier(0, currentPath.segments[0]->params[1], currentPath.segments[0]->params[3], currentPath.segments[0]->params[5], currentPath.segments[0]->params[7]);
		double previousGuess = 0;
		double nextGuess = 0;
		bool uniqueGuess = true;
		double curvature = 0;
		double length = 0;
		bool insidePath = false;
		int outsideTest = 0;
		double distance = 0;
		bool previousInside = true;
		double previousx = 0;
		double previousy = 0;
		int currentBorder = 0;
		
		for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
			for (int point = 0; point < 1 / interStep; point++) {
				previousx = pointx;
				previousy = pointy;
				pointx = Bezier(point * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
				pointy = Bezier(point * interStep, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
				vector<double> distance0;
				vector<double> distance1;
				insidePath = false;
				
				for (int borderSeg = 0; borderSeg < borders[0]->segments.size(); borderSeg++) {
					//vector<double> previousGuesses;
					for (int i = 0; i < parts; i++) {
						nextGuess = double(i) / double(parts);
						for (int n = 0; n < iterations; n++) {
							previousGuess = nextGuess;

							Fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / (2 * ISNstep);
							fx = (magicEquation(previousGuess, ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) - (2 * magicEquation(previousGuess, 0, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) + magicEquation(previousGuess, -ISNstep, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7], borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7])) / pow(ISNstep, 2);

							if (fx == 0) {
								//cout << "fx = 0 error\n";
								break;
							}

							nextGuess = previousGuess - (Fx / fx);

							if (abs(previousGuess - nextGuess) < minimumDifference) break;
						}
						
						/*
						uniqueGuess = true;
						for (int guess = 0; guess < previousGuesses.size(); guess++) {
							if (abs(nextGuess - previousGuesses[guess]) < 2 * minimumDifference) uniqueGuess = false;
						}
						
						if (uniqueGuess) {
							previousGuesses.push_back(nextGuess);
							if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
								pointsOutside++;
								outsidePath = false;
							} else if (nextGuess > 0 && nextGuess < 1) {
								distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
								distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
							}
						}
						*/
						if (xTransform(nextGuess, pointx, pointy, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[6], borders[0]->segments[borderSeg]->params[7]) < 0 && xTransform(nextGuess, pointx, pointy, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[6], borders[1]->segments[borderSeg]->params[7]) > 0 && nextGuess > 0 && nextGuess < 1) {
							pointsOutside++;
							insidePath = true;
							break;
						} else if (nextGuess > 0 && nextGuess < 1) {
							distance0.push_back(sqrt(pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[0], borders[0]->segments[borderSeg]->params[2], borders[0]->segments[borderSeg]->params[4], borders[0]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[0]->segments[borderSeg]->params[1], borders[0]->segments[borderSeg]->params[3], borders[0]->segments[borderSeg]->params[5], borders[0]->segments[borderSeg]->params[7]) - pointy, 2)));
							distance1.push_back(sqrt(pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[0], borders[1]->segments[borderSeg]->params[2], borders[1]->segments[borderSeg]->params[4], borders[1]->segments[borderSeg]->params[6]) - pointx, 2) + pow(Bezier(nextGuess, borders[1]->segments[borderSeg]->params[1], borders[1]->segments[borderSeg]->params[3], borders[1]->segments[borderSeg]->params[5], borders[1]->segments[borderSeg]->params[7]) - pointy, 2)));
						}
					}
				}
				
				if (!insidePath) {
					outsideTest++;
					int side = 0;
					double minDist0 = 9999999999999;
					double minDist1 = 9999999999999;
					
					for (int dist = 0; dist < distance0.size(); dist++) {
						if (distance0[dist] < minDist0) minDist0 = distance0[dist];
						if (distance1[dist] < minDist1) minDist1 = distance1[dist];
					}
					
					curvature += abs(BezierTangentAngle((point + 1) * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7]) - BezierTangentAngle(point * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7])) / interStep;
					length += sqrt(pow(Bezier((point + 1) * interStep, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]) - pointx, 2) + pow(Bezier((point + 1) * interStep, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]) - pointy, 2));
					
					if (previousInside) {
						distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
						
						if (lineIntersectionTest(0, pointx, pointy, previousx, previousy)) {
							cout << "intersection with border 0";
							currentBorder = 0;
						}
						
						if (lineIntersectionTest(1, pointx, pointy, previousx, previousy)) {
							cout << "intersection with border 1";
							currentBorder = 1;
						}
					}
					
					if (currentBorder == 0) distance += minDist0;
					else if (currentBorder == 1) distance += minDist1;
					
					previousInside = false;
					
				} else {
					if (!previousInside) {
						distance += sqrt(pow(pointx - previousx, 2) + pow(pointy - previousy, 2));
						
						if (lineIntersectionTest(0, pointx, pointy, previousx, previousy)) {
							cout << "intersection with border 0";
							currentBorder = 0;
						}
						
						if (lineIntersectionTest(1, pointx, pointy, previousx, previousy)) {
							cout << "intersection with border 1";
							currentBorder = 1;
						}
					}
					previousInside = true;
				}
			}
		}
		
		double outside = (currentPath.segments.size() / interStep) - pointsOutside;
		
		//cout << "Points Outside : " << outside << endl;
		//cout << "Outside Real : Outside Test = " << double(outside) / double(outsideTest) << endl;
		//cout << "Total : Outside Test + Inside = " << (path.segments.size() / interStep) / double(outsideTest + pointsOutside) << endl;
		cout << "Inside : Total = " << double(pointsOutside) / double((currentPath.segments.size() / interStep)) << endl;
		cout << "Outside : Total = " << double(outsideTest) / double((currentPath.segments.size() / interStep)) << endl;
		
		if (outside < 0) return 0;
		else return outsideMultiplier*(double(outsideTest) / double((currentPath.segments.size() / interStep))) + distanceMultiplier*pow(distance, 2) + lengthMultiplier*pow(length, 1) + curveMultiplier*pow(curvature, 2); //+ outsideMultiplier*outside + distanceMultiplier*distance;
	}
	
	void findIntersections () { //(BezierPiecewise& currentPath) {
		BezierPiecewise& currentPath = path;
		
		int pixelRange = 1;
		
		int pathInterSeg = 0;
		int borderInterSeg = 0;
		
		double pathInterTime = 0;
		double borderInterTime = 0;
		int borderInter = 0;
		
		double lastTime = -1;
		
		bool onLine = false;
		bool previouslyExiting = false;
		
		for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
			for (int border = 0; border < borders.size(); border++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int pathPixel = 0; pathPixel < currentPath.segments[pathSeg]->tCoord.size(); pathPixel++) {
						int currentPathX = currentPath.segments[pathSeg]->xCoord[pathPixel];
						int currentPathY = currentPath.segments[pathSeg]->yCoord[pathPixel];
						
						int nextPathX;
						int nextPathY;
						
						double currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel];
						onLine = false;
						
						int currentBorderX;
						int currentBorderY;
						
						int nextBorderX;
						int nextBorderY;
						
						double borderT;
						
						for (int borderPixel = 0; borderPixel < borders[border]->segments[borderSeg]->tCoord.size(); borderPixel++) {
							
							currentBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel];
							currentBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel];
							
							if ((currentPathX == currentBorderX - pixelRange || currentPathX == currentBorderX + pixelRange || currentPathX == currentBorderX) && (currentPathY == currentBorderY - pixelRange || currentPathY == currentBorderY + pixelRange || currentPathY == currentBorderY)) {
								if (borderPixel == borders[border]->segments[borderSeg]->tCoord.size() - 1) {
									if (borderSeg == borders[border]->segments.size() - 1) {
										nextBorderX = borders[border]->segments[0]->xCoord[0];
										nextBorderY = borders[border]->segments[0]->yCoord[0];
									} else {
										nextBorderX = borders[border]->segments[borderSeg + 1]->xCoord[0];
										nextBorderY = borders[border]->segments[borderSeg + 1]->yCoord[0];
									}
								} else {
									nextBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel + 1];
									nextBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel + 1];
								}
								
								//if (xTransformPixelsOutside(borders[border]->leftSide, currentPathX, currentPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
								if (xTransformBezierOutside(borders[border]->leftSide, borders[border]->segments[borderSeg]->tCoord[borderPixel], currentPathX, currentPathY, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) {
									onLine = true;
									currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel];
									borderT = borders[border]->segments[borderSeg]->tCoord[borderPixel];
									break;
								}
							}
							
							/*
							if (currentPathX == currentBorderX && currentPathY == currentBorderY) {
								outsidex.push_back(currentPathX);
								outsidey.push_back(currentPathY);
							}
							*/
						}
						
						if (onLine) {
							if (pathPixel == currentPath.segments[pathSeg]->tCoord.size() - 1) {
								if (pathSeg == currentPath.segments.size() - 1) {
									nextPathX = currentPath.segments[0]->xCoord[0];
									nextPathY = currentPath.segments[0]->yCoord[0];
								} else {
									nextPathX = currentPath.segments[pathSeg + 1]->xCoord[0];
									nextPathY = currentPath.segments[pathSeg + 1]->yCoord[0];
								}
							} else {
								nextPathX = currentPath.segments[pathSeg]->xCoord[pathPixel + 1];
								nextPathY = currentPath.segments[pathSeg]->yCoord[pathPixel + 1];
							}
							
							//if (exitingBorder(borders[border]->leftSide, currentPathX, currentPathY, nextPathX, nextPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
							if (exitingBorderBezier(borders[border]->leftSide, currentPathT, borderT, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) {
								cout << "Exiting" << endl;
								cout << "X : " << currentPathX << ", Y : " << currentPathY << endl;
								if (!previouslyExiting) {
									cout << "Counted" << endl << endl;
									outsidex.push_back(currentPathX);
									outsidey.push_back(currentPathY);
									outsidet.push_back(currentPath.segments[pathSeg]->tCoord[pathPixel]);
									lastTime = currentPathT;
									previouslyExiting = true;
									//previousPrevious = previousPrevious;
									//previouslyOnLine = true;
								}
							} else {
								cout << "Entering" << endl;
								cout << "X : " << currentPathX << ", Y : " << currentPathY << endl;
								if (previouslyExiting) {
									cout << "Counted" << endl << endl;
									outsidex.push_back(currentPathX);
									outsidey.push_back(currentPathY);
									outsidet.push_back(currentPath.segments[pathSeg]->tCoord[pathPixel]);
									lastTime = currentPathT;
									previouslyExiting = false;
								}
							}
						}
					}
				}
			}
		}
	}
	
	double ExternalAreaV2 (BezierPiecewise& currentPath) {
		//BezierPiecewise& currentPath = path;
		
		int pixelRange = 1;
		
		int pathInterSeg = 0;
		int borderInterSeg = 0;
		
		double pathInterTime = 0;
		double borderInterTime = 0;
		int borderInter = 0;
		
		double lastTime = -1;
		
		bool onLine = false;
		bool previouslyExiting = false;
		
		vector <double> exitTPath;
		vector <double> reEntryTPath;
		
		vector <double> exitTBorder;
		vector <double> reEntryTBorder;
		
		double area = 0;
		int areaSegments = 5;
		
		for (int border = 0; border < borders.size(); border++) {
			for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int pathPixel = 0; pathPixel < currentPath.segments[pathSeg]->tCoord.size(); pathPixel++) {
						int currentPathX = currentPath.segments[pathSeg]->xCoord[pathPixel];
						int currentPathY = currentPath.segments[pathSeg]->yCoord[pathPixel];
						
						int nextPathX;
						int nextPathY;
						
						double currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel];
						onLine = false;
						
						int currentBorderX;
						int currentBorderY;
						
						int nextBorderX;
						int nextBorderY;
						
						double borderT;
						
						for (int borderPixel = 0; borderPixel < borders[border]->segments[borderSeg]->tCoord.size(); borderPixel++) {
							
							currentBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel];
							currentBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel];
							
							if ((currentPathX == currentBorderX - pixelRange || currentPathX == currentBorderX + pixelRange || currentPathX == currentBorderX) && (currentPathY == currentBorderY - pixelRange || currentPathY == currentBorderY + pixelRange || currentPathY == currentBorderY)) {
								if (borderPixel == borders[border]->segments[borderSeg]->tCoord.size() - 1) {
									if (borderSeg == borders[border]->segments.size() - 1) {
										nextBorderX = borders[border]->segments[0]->xCoord[0];
										nextBorderY = borders[border]->segments[0]->yCoord[0];
									} else {
										nextBorderX = borders[border]->segments[borderSeg + 1]->xCoord[0];
										nextBorderY = borders[border]->segments[borderSeg + 1]->yCoord[0];
									}
								} else {
									nextBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel + 1];
									nextBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel + 1];
								}
								
								//if (xTransformPixelsOutside(borders[border]->leftSide, currentPathX, currentPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
								if (xTransformBezierOutside(borders[border]->leftSide, borders[border]->segments[borderSeg]->tCoord[borderPixel], currentPathX, currentPathY, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) {
									onLine = true;
									currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel];
									borderT = borders[border]->segments[borderSeg]->tCoord[borderPixel];
									break;
								}
							}
							
							/*
							if (currentPathX == currentBorderX && currentPathY == currentBorderY) {
								outsidex.push_back(currentPathX);
								outsidey.push_back(currentPathY);
							}
							*/
						}
						
						if (onLine) {
							if (pathPixel == currentPath.segments[pathSeg]->tCoord.size() - 1) {
								if (pathSeg == currentPath.segments.size() - 1) {
									nextPathX = currentPath.segments[0]->xCoord[0];
									nextPathY = currentPath.segments[0]->yCoord[0];
								} else {
									nextPathX = currentPath.segments[pathSeg + 1]->xCoord[0];
									nextPathY = currentPath.segments[pathSeg + 1]->yCoord[0];
								}
							} else {
								nextPathX = currentPath.segments[pathSeg]->xCoord[pathPixel + 1];
								nextPathY = currentPath.segments[pathSeg]->yCoord[pathPixel + 1];
							}
							
							//if (exitingBorder(borders[border]->leftSide, currentPathX, currentPathY, nextPathX, nextPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
							if (exitingBorderBezier(borders[border]->leftSide, currentPathT, borderT, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) {
								cout << "Exiting" << endl;
								cout << "X : " << currentPathX << ", Y : " << currentPathY << endl;
								if (!previouslyExiting) {
									cout << "Counted" << endl << endl;
									outsidex.push_back(currentPathX);
									outsidey.push_back(currentPathY);
									outsidet.push_back(currentPath.segments[pathSeg]->tCoord[pathPixel]);
									lastTime = currentPathT;
									previouslyExiting = true;
									exitTPath.push_back(currentPathT + pathSeg);
									exitTBorder.push_back(borderT + borderSeg);
									//previousPrevious = previousPrevious;
									//previouslyOnLine = true;
								}
							} else {
								cout << "Entering" << endl;
								cout << "X : " << currentPathX << ", Y : " << currentPathY << endl;
								if (previouslyExiting) {
									cout << "Counted" << endl << endl;
									outsidex.push_back(currentPathX);
									outsidey.push_back(currentPathY);
									outsidet.push_back(currentPath.segments[pathSeg]->tCoord[pathPixel]);
									lastTime = currentPathT;
									previouslyExiting = false;
									reEntryTPath.push_back(currentPathT + pathSeg);
									reEntryTBorder.push_back(borderT + borderSeg);
								}
							}
						}
					}
				}
			}
			
			// Calculate External Area
			if (exitTPath.size() != reEntryTPath.size()) {
				cout << "Uneven Exits & Reentries" << endl;
			} else {

				long double pathIntersectionStart;
				long double pathIntersectionEnd;
				long double borderIntersectionStart;
				long double borderIntersectionEnd;
				long double borderStep = 0;
				long double pathStep = 0;
				long double pathTempX = 0;
				long double pathTempY = 0;
				long double borderTempX = 0;
				long double borderTempY = 0;
				int pathTempSeg = 0;
				int borderTempSeg = 0;
				long double pathTempT = 0;
				long double borderTempT = 0;

				for (int region = 0; region < exitTPath.size(); region++) {

					pathIntersectionStart = exitTPath[region];
					pathIntersectionEnd = reEntryTPath[region];

					borderIntersectionStart = exitTBorder[region];
					borderIntersectionEnd = reEntryTBorder[region];

					borderStep = (borderIntersectionEnd - borderIntersectionStart) / areaSegments;
					pathStep = (pathIntersectionEnd - pathIntersectionStart) / areaSegments;

					for (int areaSeg = 1; areaSeg < areaSegments; areaSeg++) {

						pathTempSeg = floor(pathIntersectionStart + (areaSeg * borderStep));
						borderTempSeg = floor(borderIntersectionStart + (areaSeg * borderStep));

						pathTempT = (pathIntersectionStart + (areaSeg * pathStep)) - pathTempSeg;
						borderTempT = (borderIntersectionStart + (areaSeg * borderStep)) - borderTempSeg;
						pathTempX = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[0], currentPath.segments[pathTempSeg]->params[2], currentPath.segments[pathTempSeg]->params[4], currentPath.segments[pathTempSeg]->params[6]);
						pathTempY = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[1], currentPath.segments[pathTempSeg]->params[3], currentPath.segments[pathTempSeg]->params[5], currentPath.segments[pathTempSeg]->params[7]);


						borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
						borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);

						area += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2))*((abs(pathStep) + abs(borderStep)) / 2);
					}
				}
			}
		}
		
		return area;
	}
	
	long double ExternalArea (BezierPiecewise& currentPath) {
		//BezierPiecewise& currentPath = path;
		
		cout << 2 << endl;
		
		int pixelRange = 2;
		
		int pathInterSeg = 0;
		int borderInterSeg = 0;
		
		double pathInterTime = 0;
		double borderInterTime = 0;
		int borderInter = 0;
		
		double lastTime = -1;
		
		bool onLine = false;
		bool previouslyExiting = false;
		
		vector <double> exitTPath;
		vector <double> reEntryTPath;
		
		vector <double> exitTBorder;
		vector <double> reEntryTBorder;
		
		long double area = 0;
		int areaSegments = 100;
		
		for (int border = 0; border < borders.size(); border++) {
			exitTPath.clear();
			reEntryTPath.clear();
			
			exitTBorder.clear();
			reEntryTBorder.clear();
			
			for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int pathPixel = 0; pathPixel < currentPath.segments[pathSeg]->tCoord.size(); pathPixel++) {
						int currentPathX = currentPath.segments[pathSeg]->xCoord[pathPixel];
						int currentPathY = currentPath.segments[pathSeg]->yCoord[pathPixel];
						
						int nextPathX;
						int nextPathY;
						
						double currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel];
						onLine = false;
						
						int currentBorderX;
						int currentBorderY;
						
						int nextBorderX;
						int nextBorderY;
						
						double borderT;
						
						for (int borderPixel = 0; borderPixel < borders[border]->segments[borderSeg]->tCoord.size(); borderPixel++) {
							
							currentBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel];
							currentBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel];
							
							if ((currentPathX == currentBorderX - pixelRange || currentPathX == currentBorderX + pixelRange || currentPathX == currentBorderX) && (currentPathY == currentBorderY - pixelRange || currentPathY == currentBorderY + pixelRange || currentPathY == currentBorderY)) {
								if (borderPixel == borders[border]->segments[borderSeg]->tCoord.size() - 1) {
									if (borderSeg == borders[border]->segments.size() - 1) {
										nextBorderX = borders[border]->segments[0]->xCoord[0];
										nextBorderY = borders[border]->segments[0]->yCoord[0];
									} else {
										nextBorderX = borders[border]->segments[borderSeg + 1]->xCoord[0];
										nextBorderY = borders[border]->segments[borderSeg + 1]->yCoord[0];
									}
								} else {
									nextBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel + 1];
									nextBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel + 1];
								}
								
								//if (xTransformPixelsOutside(borders[border]->leftSide, currentPathX, currentPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
								if (xTransformBezierOutside(borders[border]->leftSide, borders[border]->segments[borderSeg]->tCoord[borderPixel], currentPathX, currentPathY, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) {
									onLine = true;
									currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel] + pathSeg;
									borderT = borders[border]->segments[borderSeg]->tCoord[borderPixel] + borderSeg;
									break;
								}
							}
							
							/*
							if (currentPathX == currentBorderX && currentPathY == currentBorderY) {
								outsidex.push_back(currentPathX);
								outsidey.push_back(currentPathY);
							}
							*/
						}
						
						if (onLine) {
							if (pathPixel == currentPath.segments[pathSeg]->tCoord.size() - 1) {
								if (pathSeg == currentPath.segments.size() - 1) {
									nextPathX = currentPath.segments[0]->xCoord[0];
									nextPathY = currentPath.segments[0]->yCoord[0];
								} else {
									nextPathX = currentPath.segments[pathSeg + 1]->xCoord[0];
									nextPathY = currentPath.segments[pathSeg + 1]->yCoord[0];
								}
							} else {
								nextPathX = currentPath.segments[pathSeg]->xCoord[pathPixel + 1];
								nextPathY = currentPath.segments[pathSeg]->yCoord[pathPixel + 1];
							}
							
							//if (exitingBorder(borders[border]->leftSide, currentPathX, currentPathY, nextPathX, nextPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
							if (exitingBorderBezier(borders[border]->leftSide, currentPathT, borderT, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7], borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) {
								cout << "EXITING" << endl;
								if (!previouslyExiting) {
									cout << "X : " << currentPathX << ", Y : " << currentPathY << endl;
									outsidex.push_back(currentPathX);
									outsidey.push_back(currentPathY);
									outsidet.push_back(currentPath.segments[pathSeg]->tCoord[pathPixel]);
									lastTime = currentPathT;
									previouslyExiting = true;
									exitTPath.push_back(currentPathT);
									exitTBorder.push_back(borderT);
									
									//previousPrevious = previousPrevious;
									//previouslyOnLine = true;
								}
							} else {
								cout << "REENTERING" << endl;
								if (previouslyExiting) {
									cout << "X : " << currentPathX << ", Y : " << currentPathY << endl;
									outsidex.push_back(currentPathX);
									outsidey.push_back(currentPathY);
									outsidet.push_back(currentPath.segments[pathSeg]->tCoord[pathPixel]);
									lastTime = currentPathT;
									previouslyExiting = false;
									reEntryTPath.push_back(currentPathT);
									reEntryTBorder.push_back(borderT);
								}
							}
						}
					}
				}
			}
			
			// Calculate External Area
			if (exitTPath.size() != reEntryTPath.size()) {
				cout << "Uneven Exits & Reentries" << endl;
			} else {

				long double pathIntersectionStart;
				long double pathIntersectionEnd;
				long double borderIntersectionStart;
				long double borderIntersectionEnd;
				long double borderStep = 0;
				long double pathStep = 0;
				long double pathTempX = 0;
				long double pathTempY = 0;
				long double borderTempX = 0;
				long double borderTempY = 0;
				int pathTempSeg = 0;
				int borderTempSeg = 0;
				long double pathTempT = 0;
				long double borderTempT = 0;

				for (int region = 0; region < exitTPath.size(); region++) {

					pathIntersectionStart = exitTPath[region];
					pathIntersectionEnd = reEntryTPath[region];

					borderIntersectionStart = exitTBorder[region];
					borderIntersectionEnd = reEntryTBorder[region];

					borderStep = (borderIntersectionEnd - borderIntersectionStart) / areaSegments;
					pathStep = (pathIntersectionEnd - pathIntersectionStart) / areaSegments;

					for (int areaSeg = 1; areaSeg < areaSegments; areaSeg++) {

						pathTempSeg = floor(pathIntersectionStart + (areaSeg * borderStep));
						borderTempSeg = floor(borderIntersectionStart + (areaSeg * borderStep));

						pathTempT = (pathIntersectionStart + (areaSeg * pathStep)) - pathTempSeg;
						borderTempT = (borderIntersectionStart + (areaSeg * borderStep)) - borderTempSeg;
						pathTempX = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[0], currentPath.segments[pathTempSeg]->params[2], currentPath.segments[pathTempSeg]->params[4], currentPath.segments[pathTempSeg]->params[6]);
						pathTempY = Bezier(pathTempT, currentPath.segments[pathTempSeg]->params[1], currentPath.segments[pathTempSeg]->params[3], currentPath.segments[pathTempSeg]->params[5], currentPath.segments[pathTempSeg]->params[7]);


						borderTempX = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[0], borders[border]->segments[borderTempSeg]->params[2], borders[border]->segments[borderTempSeg]->params[4], borders[border]->segments[borderTempSeg]->params[6]);
						borderTempY = Bezier(borderTempT, borders[border]->segments[borderTempSeg]->params[1], borders[border]->segments[borderTempSeg]->params[3], borders[border]->segments[borderTempSeg]->params[5], borders[border]->segments[borderTempSeg]->params[7]);

						area += sqrt(pow(borderTempX - pathTempX, 2) + pow(borderTempY - pathTempY, 2))*((abs(pathStep) + abs(borderStep)) / 2);
					}
				}
			}
			
			cout << 1 << endl;
		}
		
		return area;
	}
	
	bool Intersecting (BezierPiecewise& currentPath) {
		//BezierPiecewise& currentPath = path;
		
		currentPath.produceLines();
		
		int pixelRange = 1;
		
		int pathInterSeg = 0;
		int borderInterSeg = 0;
		
		double pathInterTime = 0;
		double borderInterTime = 0;
		int borderInter = 0;
		
		double lastTime = -1;
		
		bool onLine = false;
		bool previouslyExiting = false;
		
		for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
			for (int border = 0; border < borders.size(); border++) {
				for (int borderSeg = 0; borderSeg < borders[border]->segments.size(); borderSeg++) {
					for (int pathPixel = 0; pathPixel < currentPath.segments[pathSeg]->tCoord.size(); pathPixel++) {
						int currentPathX = currentPath.segments[pathSeg]->xCoord[pathPixel];
						int currentPathY = currentPath.segments[pathSeg]->yCoord[pathPixel];
						
						int nextPathX;
						int nextPathY;
						
						double currentPathT = currentPath.segments[pathSeg]->tCoord[pathPixel];
						onLine = false;
						
						int currentBorderX;
						int currentBorderY;
						
						int nextBorderX;
						int nextBorderY;
						
						double borderT;
						
						for (int borderPixel = 0; borderPixel < borders[border]->segments[borderSeg]->tCoord.size(); borderPixel++) {
							
							currentBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel];
							currentBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel];
							
							if ((currentPathX == currentBorderX - pixelRange || currentPathX == currentBorderX + pixelRange || currentPathX == currentBorderX) && (currentPathY == currentBorderY - pixelRange || currentPathY == currentBorderY + pixelRange || currentPathY == currentBorderY)) {
								if (borderPixel == borders[border]->segments[borderSeg]->tCoord.size() - 1) {
									if (borderSeg == borders[border]->segments.size() - 1) {
										nextBorderX = borders[border]->segments[0]->xCoord[0];
										nextBorderY = borders[border]->segments[0]->yCoord[0];
									} else {
										nextBorderX = borders[border]->segments[borderSeg + 1]->xCoord[0];
										nextBorderY = borders[border]->segments[borderSeg + 1]->yCoord[0];
									}
								} else {
									nextBorderX = borders[border]->segments[borderSeg]->xCoord[borderPixel + 1];
									nextBorderY = borders[border]->segments[borderSeg]->yCoord[borderPixel + 1];
								}
								
								//if (xTransformPixelsOutside(borders[border]->leftSide, currentPathX, currentPathY, currentBorderX, currentBorderY, nextBorderX, nextBorderY)) {
								if (xTransformBezierOutside(borders[border]->leftSide, borders[border]->segments[borderSeg]->tCoord[borderPixel], currentPathX, currentPathY, borders[border]->segments[borderSeg]->params[0], borders[border]->segments[borderSeg]->params[1], borders[border]->segments[borderSeg]->params[2], borders[border]->segments[borderSeg]->params[3], borders[border]->segments[borderSeg]->params[4], borders[border]->segments[borderSeg]->params[5], borders[border]->segments[borderSeg]->params[6], borders[border]->segments[borderSeg]->params[7])) return true;
							}
						}
					}
				}
			}
		}
		
		return false;
	}
	
	double fullEquationT (BezierPiecewise& currentPath) {
		double t = 0;
		
		if (Intersecting(currentPath)) return 999999999999;
		
		for (int pathSeg = 0; pathSeg < currentPath.segments.size(); pathSeg++) {
			for (int i = 0; i < 1 / step; i++) {
				t += timeTakenCalc(currentPath, pathSeg, i * step);
			}
		}
		
		return t;
	}
		
	double timeTakenCalc (BezierPiecewise& currentPath, int pathSeg, double t) {
		long double currentAngle = BezierTangentAngle(t, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7]);
		long double nextAngle = BezierTangentAngle(t + step, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[6], currentPath.segments[pathSeg]->params[7]);

		double x0 = Bezier(t, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
		double y0 = Bezier(t, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
		
		double x1 = Bezier(t + step, currentPath.segments[pathSeg]->params[0], currentPath.segments[pathSeg]->params[2], currentPath.segments[pathSeg]->params[4], currentPath.segments[pathSeg]->params[6]);
		double y1 = Bezier(t + step, currentPath.segments[pathSeg]->params[1], currentPath.segments[pathSeg]->params[3], currentPath.segments[pathSeg]->params[5], currentPath.segments[pathSeg]->params[7]);
		
		double s = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
		
		double v = sqrt( (2 * s * mew * m * g) / (m * pow(nextAngle - currentAngle, 2)) );
		
		if (v <= Vmax) return s / v;
		else return s / Vmax;
	}
		

};
	
int main () {
	/*
	cout << setprecision(30);
	cout << Bezier(0.5, 7.56, -0.07, 9.32, 7.8) << endl;
	cout << dBezierdT(0.5, 7.56, -0.07, 9.32, 7.8) << endl;
	cout << BezierNormalAngle(0.75, 2, 0, 3, 0, 2, 2, 6, 0) << endl;
	cout << distanceSigned(0.7747, 3.97, 1.52, 2, 0, 3, 0, 2, 2, 6, 0) << endl;
	cout << distanceSigned(0.7747, 3.97, 1.52, 2, 1, 3, 1, 2, 3, 6, 0.5) << endl;
	cout << inside(0.7747, 3.89, 1.767, 2, 0, 3, 0, 2, 2, 6, 0, 2, 1, 3, 1, 2, 3, 6, 0.5) << endl;
	cout << (xTransform(0.7747, 3.89, 1.767, 2, 0, 3, 0, 2, 2, 6, 0) < 0 && xTransform(0.7725, 3.89, 1.767, 2, 1, 3, 1, 2, 3, 6, 0.5) > 0) << endl << endl;
	int sum = 0;
	for (int i = 1; i < 1 / step; i++) {
		if (inside(double(i*step), 3.89, 1.767, 2, 0, 3, 0, 2, 2, 6, 0, 2, 1, 3, 1, 2, 3, 6, 0.5)) cout << i*step << endl;
	}
	
	cout << endl;
	
	bool a;
	bool b;
	
	auto start = chrono::high_resolution_clock::now();
	//a = insideSegment(2.1, 0.5, 4, 2, 4, 0.47, 5.25, 0.48, 2, 0, 3, 0, 2, 2, 6, 0, 2, 1, 3, 1, 2, 3, 6, 0.5);
	auto stop = chrono::high_resolution_clock::now();
	
	auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
	cout << 0.000000001*duration.count() << endl;
	cout << a << endl << endl;
	
	start = chrono::high_resolution_clock::now();
	b = lineSegmentTest(2.1, 0.5, 4, 2, 4, 0.47, 5.25, 0.48, 2, 0, 3, 0, 2, 2, 6, 0, 2, 1, 3, 1, 2, 3, 6, 0.5);
	stop = chrono::high_resolution_clock::now();
	
	duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
	cout << 0.000000001*duration.count() << endl;
	cout << b << endl << endl;
	
	BezierSegment b1;
	b1.setupFirst(6.52, 5.14, 17.86, -1.7, 8.8, 15.63, 7.8, 9.14);
	cout << b1.Curvature() << endl;
	
	BezierPiecewise right;
	right.firstSegment(2, 0, 3, 0, 2, 2, 6, 0);
	right.addSegment(8, 2, 0, 4, 1);
	right.addSegment(-6, 0, 0, -5, 1);
	right.addSegment(9.74, 1.07, 11.037, 0.12, 1);
	right.addSegment(6, 5, -2, 10, 1);
	
	BezierPiecewise left;
	left.firstSegment(2, 1, 3, 1, 2, 3, 6, 0.5);
	left.addSegment(7, 2, 0, 3, 1);
	left.addSegment(-4, 1, 0, -2, 1);
	left.addSegment(11, 0, 10, 0.5, 1);
	left.addSegment(6, 2, -2, 5, 1);
	
	BezierPiecewise path;
	path.firstSegment(2, 0.5, 3, 2, 3.83, 1.3, 6.7, -0.17);
	path.addSegment(7.74, 1.85, -0.91, 3.56, 1);
	path.addSegment(-1.77, -1.89, 0.67, -3.93, 1);
	path.addSegment(10.91, -0.26, 10.4, 0.48, 1);
	path.addSegment(6, 2, -4.23, 8.59, 1);
	
	cout << "testing\n";
	
	for (int seg = 0; seg < path.segments.size(); seg++) {
		for (int bseg = 0; bseg < right.segments.size(); bseg++) {
			cout << notOverlapping(path.segments[seg]->ax, path.segments[seg]->ay, path.segments[seg]->bx, path.segments[seg]->by, path.segments[seg]->cx, path.segments[seg]->cy, path.segments[seg]->dx, path.segments[seg]->dy, right.segments[bseg]->ax, right.segments[bseg]->ay, right.segments[bseg]->bx, right.segments[bseg]->by, right.segments[bseg]->cx, right.segments[bseg]->cy, right.segments[bseg]->dx, right.segments[bseg]->dy) << ", ";
		}
		
		cout << endl; 
		
		for (int bseg = 0; bseg < left.segments.size(); bseg++) {
			cout << notOverlapping(path.segments[seg]->ax, path.segments[seg]->ay, path.segments[seg]->bx, path.segments[seg]->by, path.segments[seg]->cx, path.segments[seg]->cy, path.segments[seg]->dx, path.segments[seg]->dy, left.segments[bseg]->ax, left.segments[bseg]->ay, left.segments[bseg]->bx, left.segments[bseg]->by, left.segments[bseg]->cx, left.segments[bseg]->cy, left.segments[bseg]->dx, left.segments[bseg]->dy) << ", ";
		}
		
		cout << endl << endl;
	}
	
	// Minimise 
	for (int i = 0; i < 10; i++) {
	}
	*/
	
	handler h;
	h.addBorder();
	h.borders[0]->firstSegment(200, 0, 300, 0, 200, 200, 600, 0);
	h.borders[0]->addSegment(800, 200, 0, 400, 1);
	h.borders[0]->addSegment(-600, 0, 0, -500, 1);
	h.borders[0]->addSegment(1000, 0, 600, -50, 1);
	h.borders[0]->addSegment(150, 0, 200, 0, 1);
	h.borders[0]->produceLines();
	h.borders[0]->leftSide = false;
	
	h.addBorder();
	h.borders[1]->firstSegment(200, 100, 300, 100, 200, 300, 600, 50);
	h.borders[1]->addSegment(700, 200, 0, 300, 1);
	h.borders[1]->addSegment(-400, 100, 0, -200, 1);
	h.borders[1]->addSegment(1100, 0, 400, -100, 1);
	h.borders[1]->addSegment(100, 100, 200, 100, 1);
	h.borders[1]->produceLines();
	h.borders[1]->leftSide = true;
	
	/*
	h.path.firstSegment(30, 0, 30, 20, 40, 10, 67, -1.6);
	h.path.addSegment(77.5, 17.5, -18, 41.8, 1);
	h.path.addSegment(-20, -20, 10, -40, 1);
	h.path.addSegment(110, 0,77.3, -20.4, 1);
	h.path.addSegment(10, 0, 30, 0, 1);
	*/
	
	h.path.firstSegment(200, 50, 300, 200, 400, 100, 670, -16);
	h.path.addSegment(775, 175, -100, 350, 1);
	h.path.addSegment(-200, -200, 100, -400, 1);
	h.path.addSegment(1100, 0, 500, -75, 1);
	h.path.addSegment(100, 0, 200, 50, 1);
	h.path.produceLines();
	
	h.pathTemp0.firstSegment(20, 5, 30, 20, 40, 10, 67, -1.6);
	h.pathTemp0.addSegment(77.5, 17.5, 13, 25, 1);
	h.pathTemp0.addSegment(-20, -20, 10, -40, 1);
	h.pathTemp0.addSegment(110, 0, 77.3, -20.4, 1);
	h.pathTemp0.addSegment(10, 0, 20, 5, 1);
	
	h.pathTemp1.firstSegment(20, 5, 30, 20, 40, 10, 67, -1.6);
	h.pathTemp1.addSegment(77.5, 17.5, -10, 35, 1);
	h.pathTemp1.addSegment(-20, -20, 10, -40, 1);
	h.pathTemp1.addSegment(110, 0, 50, -7.5, 1);
	h.pathTemp1.addSegment(10, 0, 20, 5, 1);
	
	h.pathNext.firstSegment(20, 5, 30, 20, 40, 10, 67, -1.6);
	h.pathNext.addSegment(77.5, 17.5, -10, 35, 1);
	h.pathNext.addSegment(-20, -20, 10, -40, 1);
	h.pathNext.addSegment(110, 0, 50, -7.5, 1);
	h.pathNext.addSegment(10, 0, 20, 5, 1);
	
	//cout << "full equ " << h.fullEquation() << endl;
	
	//cout << "collision " << h.collisionTest() << endl;
	
	//cout << "collision f " << h.collisionTestForward() << endl;
	
	//cout << "collision b " << h.collisionTestBackward() << endl;
	
	//cout << "Area 0\n" << h.Area(0, 0) << endl;
	//cout << "Area 1\n" << h.Area(0, 1) << endl;
	
	//h.newtonsParameterUpdateV2();
	
	//h.manualIteration();
	
	//cout << "Points Outside : " << h.insideSegmentNewtons() << endl;
	
	//h.newtonsParameterUpdateV3(5);
	
	//h.newtonSingleParameter(50, 1, 4);
	
	//h.newtonSingleParamIterator(3, 1);
	
	h.manualSingleParamIterator(10, 10);
	
	//cout << h.OverlappingV3() << endl;
	
	//h.insideSegmentNewtonsRecord();
	
	//h.insideSegmentNewtonsV2();
	
	//cout << h.lineIntersectionIterator();
	
	//cout << h.dualLineIntersectionTest();
	
	//h.findIntersections();
	
	//cout << "External Area = " << h.ExternalAreaV2() << endl;
	
	h.saveToFile();
	
	return 0;
}