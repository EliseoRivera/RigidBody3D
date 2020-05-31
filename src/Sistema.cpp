#include "Sistema.h"

Sistema::Sistema()
{
    //ctor
}

Sistema::~Sistema()
{
    //dtor
}
void Sistema::InitRungeKutta(){
m=1;
b=3;
k=2;
Y.entry(0,0)=PI/6.0;
Y.entry(1,0)=0;



}
void Sistema::inicializar(){
    ///y1=q, y2=vq;
dh=0.01; tr=0.0f; //paso de tiempo tiempo en runge kutta

///
Y.zero(2,1);
Yv.zero(2,1);

Y0.zero(2,1); ///condiciones iniciales
Yv0.zero(2,1);

F1.zero(2,1);
F2.zero(2,1);
F3.zero(2,1);
F4.zero(2,1);


bloque=new modelo3D();
bloque->leer("bloque.STL");

InitRungeKutta();
}
void Sistema::RungeKuttaIntegrar(){
F1=dh*Frk(Y,tr);
F2=dh*Frk(Y+0.5f*F1,tr+0.5f*dh);
F3=dh*Frk(Y+0.5f*F2,tr+0.5f*dh);
F4=dh*Frk(Y+F3,tr+dh);
Y=Y+ 0.166666666666666666666*(F1+2.0f*F2+2.0f*F3+F4);
};
void Sistema::RungeKuttaUpdate(){

};
Matrix Sistema::Frk(Matrix Y_,float t_){ ///Y 2N x 1
Yv.entry(0,0)=Y_.entry(1,0); ///velocidad

Matrix A_(Grk(Y_,t_)); ///SEGUNDA DERIVADA DESPEJADA (N x1)
//cout<<" A_"<<endl;
//A_.mostrar();
Yv.entry(1,0)=A_.entry(0,0); ///aceleracion


return Yv;
};
Matrix Sistema::Grk(Matrix Y_,float t_){ ///N x 1
///la aceleración depende de la posición y posiblemente de la velocidad
Matrix aq(1,1);
float a;

a=-(9.81/0.5)*sin(Y_.entry(0,0))-2*Y_.entry(1,0);//(M(Y_).inverse1())*(qe+qv);
aq.entry(0,0)=a;
return aq;

};


void Sistema::renderizar(){
RungeKuttaIntegrar();
glPointSize(50);
glBegin(GL_POINTS);
glVertex3f(0.5*sin(Y.entry(0,0)),0,-0.5*cos(Y.entry(0,0)));
glEnd();
 tr=tr+dh;
}
