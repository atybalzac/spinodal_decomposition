#include <iostream>
#include "spinodal.hpp"


const double GM = 1;
const double Gh = 1;
const double Gdt = .2;
const double Geps2 = 1;

double Gc[XSIZE][YSIZE]; 


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


//calculates the first space derivative of c
double firstderiv(int i, int j){
    double xderiv = (Gc[checkbc(i+1, XSIZE)][j] - Gc[checkbc(i-1,XSIZE)][j]);
    double yderiv = (Gc[i][checkbc(j+1,YSIZE)] - Gc[i][checkbc(j-1,YSIZE)]);
    return (xderiv+yderiv)/(2*Gh);
}

//calculates second space deriv(laplacian) of c

double laplac(int i, int j){
    double xderiv = Gc[checkbc(i+1,XSIZE)][j] - 2 * Gc[i][j] + Gc[checkbc(i-1,XSIZE)][j];
    double yderiv = Gc[i][checkbc(j+1,YSIZE)] - 2 * Gc[i][j] + Gc[i][checkbc(j-1,YSIZE)];
    return (xderiv+yderiv)/(Gh*Gh);
}

//calculates third space deriv of c

double thirderiv(int i, int j){
    double xderiv = Gc[checkbc(i+2,XSIZE)][j] - 2*Gc[checkbc(i+1,XSIZE)][j] + 2*Gc[checkbc(i-1,XSIZE)][j] - Gc[checkbc(i-2,XSIZE)][j];
    double yderiv = Gc[i][checkbc(j+2,YSIZE)] - 2*Gc[i][checkbc(j+1,YSIZE)] + 2*Gc[i][checkbc(j-1,YSIZE)] - Gc[i][checkbc(j-2,YSIZE)];
    return (xderiv+yderiv)/(2*Gh*Gh*Gh);
}