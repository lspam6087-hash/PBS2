#include "vec3d.h"
#include <math.h>

// typedef struct Vec3D //define Vec3D
// {
//     double x;
//     double y;
//     double z;
// } Vec3D;

//display as {x,y,z}
Vec3D v3(double x , double y , double z ) {
return (Vec3D){x , y , z } ;
}

//addition of vectors
Vec3D add(Vec3D a , Vec3D b) {
return v3 ( a . x + b . x , a . y + b . y , a . z + b . z ) ;
}

//subtraction of vectors
Vec3D sub(Vec3D a , Vec3D b ) {
return v3 ( a . x - b . x , a . y - b . y , a . z - b . z ) ;
}

//scale vector by s
Vec3D scl(double s , Vec3D a ) {
return v3 ( s * a . x , s * a . y , s * a . z ) ;
}

//dot product of vectors
double dot(Vec3D a , Vec3D b ) {
return a . x * b . x + a . y * b . y + a . z * b . z ;
}

//cross product of vectors
Vec3D cross(Vec3D a , Vec3D b ) {
return v3 ( a . y * b . z - a . z * b . y ,
            a . z * b . x - a . x * b . z ,
            a . x * b . y - a . y * b . x ) ;
}

//norm or absolute value of vector
double norm(Vec3D a ) {
return sqrt ( dot ( a , a ) ) ;
}

//distance between two vectors
double dist(Vec3D a , Vec3D b ) {
return sqrt ( (a.x - b.x ) * (a.x - b.x ) +
              (a.y - b.y ) * (a.y - b.y ) +
              (a.z - b.z ) * (a.z - b.z ) ) ;
}