#ifndef SISTEMA_H
#define SISTEMA_H
#include "modelo3D.h"
class Sistema
{
    public:
        Sistema();
        virtual ~Sistema();
 modelo3D *bloque;
void inicializar();
void renderizar();
float dh,tr;
void InitRungeKutta();
Matrix Y,Yv,Y0,Yv0; ///aq=G(q,vq,t);
Matrix F1,F2,F3,F4;
Matrix f1,f2,f3,f4; //furzas
void RungeKuttaIntegrar();
void RungeKuttaUpdate();
Matrix Frk(Matrix Y_,float t_);///2N x 1
Matrix Grk(Matrix Y_,float t_); ///N x 1
float m; //masa
float k; ///constante de hooke
float b; ///amortiguamieno
};

#endif // SISTEMA_H
