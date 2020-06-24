//dd random number number generator to generate phase field with compositions ~.5 but between 0.48-0.52 kind of idea

#include <iostream>
#include "spinodal.hpp"
#include <stdio.h>
#include "EasyBMP.h"
#include <linux/limits.h>

int main(int argc, char **argv){
    register int i, j;

    double endtime = 100;
    double time;

    //seed random number generator
    FILE *seed = fopen("/dev/urandom", "r");
    unsigned int random_data;
    size_t seedsize = fread(&random_data, sizeof(random_data), 1, seed);
    fclose(seed);

    srand(random_data);

    //initialize Gc
    for(i = 0; i < XSIZE; i++){
        for (j = 0; j < YSIZE; j++){
            double new_min = 0.48;
            double new_max = 0.52;
            double new_range = new_max - new_min;
            double num = ((rand()*(new_range))/RAND_MAX)+new_min;
            Gc[i][j]=num;
        }

    }

    //Generate a blank image that can be updated and output as an image as the program runs
    BMP* img = new BMP();
    img -> SetSize(XSIZE,YSIZE);
    //Declare a path length for naming the image files
    char imgpath[PATH_MAX];

    //Advance through time
    for (time = 0; time <= endtime; time += Gdt){
        for (i = 0; i < XSIZE; i++){
            for(j = 0; j < YSIZE; j++){
                double newc = Gc[i][j] + GM * (firstderiv(i,j) -2 * Gc[i][j] * firstderiv(i,j)) * (firstderiv(i,j) * (2*Gc[i][j] * firstderiv(i,j)) - Geps2 * thirderiv(i,j));
                //set pixel colors, can replace with commented lines for grayscale instead of blue/pink
                img -> SetPixel(i, j, (RGBApixel){
                    .Blue = 27, //(ebmpBYTE)(255*Gphi[i][j]),
                    .Green = 163, //(ebmpBYTE)(255*Gphi[i][j]),
                    .Red = (ebmpBYTE)(((Gc[i][j]-0.48)*255.0)/(0.52-0.48)),
                    .Alpha = 0,
                });
                Gc[i][j] = newc;
            }
        }
       
        //write image to file iterating name based on timestep. Numbers are padded with leading zeroes so they sort correctly
        if (argc > 1){ 
            snprintf(imgpath, sizeof(imgpath), "%s_t=%03f.bmp", argv[1], time);
             img -> WriteToFile(imgpath);
           
        }
    }


}