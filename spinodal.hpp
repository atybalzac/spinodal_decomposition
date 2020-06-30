#pragma once
#define XSIZE 100
#define YSIZE 100

struct Vec {
    double x;   // x component of a 2D vector
    double y;   // y component of a 2D vector
};

// Global variables

extern const double Gh;
extern const double Gdt;

extern const double GM;
extern const double Geps2;
extern const double GFluct;

extern double Gc[XSIZE][YSIZE];
extern struct Vec Gfield[XSIZE][YSIZE];

// Function declarations
int checkbc(int i, int lambda);
struct Vec grad(double func[XSIZE][YSIZE], int i, int j);
double laplac(double func[XSIZE][YSIZE], int i, int j);
double div(struct Vec v[XSIZE][YSIZE], int i, int j);
