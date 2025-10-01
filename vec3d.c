#include "structs.h"
#include "vec3d.h"

//display as {x,y,z}
struct Vec3D v3(double x , double y , double z ) {
return (struct Vec3D){x , y , z } ;
}

//addition of vectors
struct Vec3D add(struct Vec3D a , struct Vec3D b) {
return v3 ( a . x + b . x , a . y + b . y , a . z + b . z ) ;
}

//subtraction of vectors
struct Vec3D sub(struct Vec3D a , struct Vec3D b ) {
return v3 ( a . x - b . x , a . y - b . y , a . z - b . z ) ;
}

//scale vector by s
struct Vec3D scl(double s , struct Vec3D a ) {
return v3 ( s * a . x , s * a . y , s * a . z ) ;
}

//dot product of vectors
double dot(struct Vec3D a , struct Vec3D b ) {
return a . x * b . x + a . y * b . y + a . z * b . z ;
}

//cross product of vectors
struct Vec3D cross(struct Vec3D a , struct Vec3D b ) {
return v3 ( a . y * b . z - a . z * b . y ,
            a . z * b . x - a . x * b . z ,
            a . x * b . y - a . y * b . x ) ;
}

//norm or absolute value of vector
double norm(struct Vec3D a ) {
return sqrt ( dot ( a , a ) ) ;
}

//distance between two vectors
double dist(struct Vec3D a , struct Vec3D b ) {
return sqrt ( (a.x - b.x ) * (a.x - b.x ) +
              (a.y - b.y ) * (a.y - b.y ) +
              (a.z - b.z ) * (a.z - b.z ) ) ;
}