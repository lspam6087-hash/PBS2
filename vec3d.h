#ifndef VEC3D_H_
#define VEC3D_H_

#include <math.h>

typedef struct Vec3D //define Vec3D
{
    double x;
    double y;
    double z;
} Vec3D;

//display as {x,y,z}
Vec3D v3 (double x , double y , double z );

//addition of vectors
Vec3D add(Vec3D a , Vec3D b);

//subtraction of vectors
Vec3D sub(Vec3D a , Vec3D b );

//scale vector by s
Vec3D scl (double s , Vec3D a );

//dot product of vectors
double dot (Vec3D a , Vec3D b );

//cross product of vectors
Vec3D cross (Vec3D a , Vec3D b );

//norm or absolute value of vector
double norm (Vec3D a );

//distance between two vectors
double dist (Vec3D a , Vec3D b );

#endif