#pragma once
#define XSIZE 100
#define YSIZE 100

extern const double Gh;
extern const double Gdt;

extern const double GM;
extern const double Geps2;

extern double Gc[XSIZE][YSIZE];
double firstderiv(int i, int j);
double laplac(int i, int j);
double thirderiv(int i, int j);