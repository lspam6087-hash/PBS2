#ifndef VEC3D_H_
#define VEC3D_H_

#include <math.h>
#include "structs.h"

//display as {x,y,z}
struct Vec3D v3 (double x , double y , double z );

//addition of vectors
struct Vec3D add(struct Vec3D a , struct Vec3D b);

//subtraction of vectors
struct Vec3D sub(struct Vec3D a , struct Vec3D b );

//scale vector by s
struct Vec3D scl (double s , struct Vec3D a );

//dot product of vectors
double dot (struct Vec3D a , struct Vec3D b );

//cross product of vectors
struct Vec3D cross (struct Vec3D a , struct Vec3D b );

//norm or absolute value of vector
double norm (struct Vec3D a );

//distance between two vectors
double dist (struct Vec3D a , struct Vec3D b );

#endif