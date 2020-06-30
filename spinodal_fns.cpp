#include <iostream>
#include "spinodal.hpp"

// Global Variables

const double GM = 6.0;
const double Gh = 1;
const double Gdt = .001;
const double Geps2 = 0.5;
const double GFluct = 0.02;

double Gc[XSIZE][YSIZE];          /* concentration field */
struct Vec Gfield[XSIZE][YSIZE];  /* chemical potential field */

// Function definitions

//Checks whether we're at the boundary and loops back to the other side if we are for periodic boundary
int checkbc(int i, int lambda){

    if (i > (lambda - 1)){
        return i-lambda;
    }
    else if (i < 0){
        return i+lambda;
    }
    else {
        return i;
    }
}

//calculates the gradient of c (a vector)
struct Vec grad(double func[XSIZE][YSIZE], int i, int j) {
    struct Vec result;

    result.x = 0.0;
    result.y = 0.0;

    result.x = (func[checkbc(i+1, XSIZE)][j] - func[checkbc(i-1,XSIZE)][j]) / (2.0 * Gh);
    result.y = (func[i][checkbc(j+1,YSIZE)] - func[i][checkbc(j-1,YSIZE)]) / (2.0 * Gh);
    return (result);
}

// calculates laplacian of c

double laplac(double func[XSIZE][YSIZE], int i, int j) {
    double xderiv = func[checkbc(i+1,XSIZE)][j] - 2 * func[i][j] + func[checkbc(i-1,XSIZE)][j];
    double yderiv = func[i][checkbc(j+1,YSIZE)] - 2 * func[i][j] + func[i][checkbc(j-1,YSIZE)];
    return (xderiv+yderiv)/(Gh*Gh);
}

//calculates the divergence of a vector field at a point

double div(struct Vec v[XSIZE][YSIZE], int i, int j) {

    double diverg = 0.0;

    diverg = (v[checkbc(i+1, XSIZE)][j].x - v[checkbc(i-1,XSIZE)][j].x) / (2.0 * Gh);
    diverg += (v[i][checkbc(j+1,YSIZE)].y - v[i][checkbc(j-1,YSIZE)].y) / (2.0 * Gh);
    return (diverg);
}
