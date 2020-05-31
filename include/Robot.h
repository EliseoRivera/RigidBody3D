#ifndef ROBOT_H
#define ROBOT_H
#include "modelo3D.h"
#include<vector>
#include <cstdlib>
///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
class Robot
{
    public:
        Robot();
        ~Robot();

        modelo3D *placa;


void inicializar();
void renderizar();
void configurarTH();

void AplicarTHx(float theta, vector3d d);
void AplicarTHy(float theta, vector3d d);
void AplicarTHz(float theta, vector3d d);
Matrix THx,THy,THz,TH;

std::vector<Matrix> THList;
std::vector<vector3d> Origenes;
std::vector<modelo3D*> modelos;


/////////////////
float dh,tr;
void InitRungeKutta();
Matrix Y,Yv,Y0,Yv0; ///aq=G(q,vq,t);
Matrix F1,F2,F3,F4;
Matrix f1,f2,f3,f4; //furzas
vector3d p1,p2,p3,p4; //posición de punto de aplicación  cuerdas
Matrix a1,a2,a3,a4;///origen cuerda arriba
Matrix b1,b2,b3,b4;///origen cuerda abajo
float lc1,lc2,lc3,lc4; ///longitud cuerda sin resorte
float lt1,lt2,lt3,lt4; ///longitud cuerda con resorte
float le,k;
Matrix  u1,u2,u3,u4;
//dynamics
float phi,theta,psi,g,mass; //
Matrix R,T,Il;
Matrix Qvth(Matrix y);
Matrix QeR(Matrix y);
Matrix Qeth(Matrix y);
Matrix P(vector3d &p);
Matrix A(Matrix &y);
Matrix G(Matrix &y);
Matrix w(Matrix &y);
Matrix wl(Matrix &y);
Matrix dtheta(Matrix &y);
Matrix I(Matrix &y);
Matrix mth(Matrix &y);
Matrix Gl(Matrix &y);
Matrix dGl(Matrix &y);
Matrix M(Matrix &y);
/////////////////////////////////
///y1=q, y2=vq
void RungeKuttaIntegrar();
void RungeKuttaUpdate();
Matrix Frk(Matrix Y_,float t_);///12x1
Matrix Grk(Matrix Y_,float t_); ///6x1
/////////////////////////

private :
void DefinirTHx(float theta, vector3d d);
void DefinirTHy(float theta, vector3d d);
void DefinirTHz(float theta, vector3d d);
void  Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R) ;


};

#endif // ROBOT_H
