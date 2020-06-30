//dd random number number generator to generate phase field with compositions ~.5 but between 0.48-0.52 kind of idea

#include <iostream>
#include "spinodal.hpp"
#include <stdio.h>
#include "EasyBMP.h"
#include <limits.h>

int main(int argc, char **argv){
    
    register int i, j;
    double endtime = 1000.0;
    double dcdt,conc,time,prefactor,bulkderiv;
    double temp[XSIZE][YSIZE],chempot[XSIZE][YSIZE];
    char oldtime[128],newtime[128];
    struct Vec term[XSIZE][YSIZE];

    //seed random number generator
    FILE *seed = fopen("/dev/urandom", "r");
    unsigned int random_data;
    size_t seedsize = fread(&random_data, sizeof(random_data), 1, seed);
    fclose(seed);

    srand(random_data);

    //initialize Gc
    double new_min = -GFluct;
    double new_max = GFluct;
    double new_range = new_max - new_min;
    double num = 0.0;

    for(i = 0; i < XSIZE; i++){
        for (j = 0; j < YSIZE; j++){
            Gc[i][j] = ((rand()*(new_range))/RAND_MAX)+new_min;
            Gfield[i][j].x = 0.0;
            Gfield[i][j].y = 0.0;
        }

    }

    //Generate a blank image that can be updated and output as an image as the program runs
    BMP* img = new BMP();
    img -> SetSize(XSIZE,YSIZE);
    //Declare a path length for naming the image files
    char imgpath[PATH_MAX];

    //Advance through time
    sprintf(oldtime,"0");
    for (time = 0; time <= endtime; time += Gdt){
        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                conc = Gc[i][j];
                bulkderiv  = conc * ((conc*conc) - 1.0);
                chempot[i][j] = bulkderiv - (Geps2 * laplac(Gc,i,j));
            }
        }

        // Now the inner scalar term is calculated everywhere, we next need
        // to create the inner vector field

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                Gfield[i][j] = grad(chempot,i,j);
                Gfield[i][j].x = Gfield[i][j].x;
                Gfield[i][j].y = Gfield[i][j].y;
            }
        }

        // Now all we need to do is take the divergence of the vector field Gfield,
        // already calculated at each point above, and then multiply by the mobility
       
        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){
                dcdt = GM * div(Gfield,i,j);
                temp[i][j] = Gc[i][j] + (dcdt * Gdt);
            }
        }

        // The new concentrations are stored in the temp array. Transfer them to the
        // concentration array.

        for (i = 0; i < XSIZE; i++){
            for (j = 0; j < YSIZE; j++){

                //set pixel colors, can replace with commented lines for grayscale instead of blue/pink
                img -> SetPixel(i, j, (RGBApixel){
                    .Blue = 27, //(ebmpBYTE)(255*Gphi[i][j]),
                    .Green = 163, //(ebmpBYTE)(255*Gphi[i][j]),
                    .Red = (ebmpBYTE)(((Gc[i][j]+2.0)*255.0)/(4.0)),
                    .Alpha = 0,
                });

                Gc[i][j] = temp[i][j];
            }
        }
       
        //write image to file iterating name based on timestep. Numbers are padded with leading zeroes so they sort correctly
        sprintf(newtime,"%d",(int)(time));
        if (argc > 1 && (strcmp(oldtime,newtime) || time < Gdt)) { 
            snprintf(imgpath, sizeof(imgpath), "%s_t=%04d.bmp", argv[1], atoi(newtime));
             img -> WriteToFile(imgpath);
            strcpy(oldtime,newtime);
        }
    }

    exit(0);

}
