#include "Robot.h"
#define PI 3.1415926535897932384626433832795
///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
Robot::Robot()
{

   THx.identity(4);
   THy.identity(4);
   THz.identity(4);
   int d=4;
   TH.identity(4);


}

Robot::~Robot()
{

delete placa;

    //dtor
}
Matrix Robot::A(Matrix &y){ //phi, theta, psi

      float a,b, c;  //phi, theta, psi
    a=y.entry(3,0);b=y.entry(4,0);c=y.entry(5,0);
Matrix m(3,3);

m.entry(0,0)=(cos(c)*cos(a)-cos(b)*sin(a)*sin(c));
m.entry(0,1)=-sin(c)*cos(a)-cos(b)*sin(a)*cos(c);
m.entry(0,2)=sin(b)*sin(a);

m.entry(1,0)=cos(c)*sin(a)+cos(b)*cos(a)*sin(c);
m.entry(1,1)=-sin(c)*sin(a)+cos(b)*cos(a)*cos(c);
m.entry(1,2)=-sin(b)*cos(a);

m.entry(2,0)=sin(b)*sin(c);
m.entry(2,1)=sin(b)*cos(c);
m.entry(2,2)=cos(b);
return m;
};
Matrix Robot::P(vector3d &p){

Matrix temp(3,1);

temp.entry(0,0)=p.x;
temp.entry(1,0)=p.y;
temp.entry(2,0)=p.z;

return A(Y)*temp;
}

Matrix Robot::G(Matrix &y){ //phi, theta, psi
    float a, b,  c;  //phi, theta, psi
    a=y.entry(3,0);b=y.entry(4,0);c=y.entry(5,0);

Matrix m(3,3);
m.entry(0,0)=0;
m.entry(0,1)=cos(a);
m.entry(0,2)=sin(b)*sin(a);

m.entry(1,0)=0;
m.entry(1,1)=sin(a);
m.entry(1,2)=-sin(b)*cos(a);

m.entry(2,0)=1;
m.entry(2,1)=0;
m.entry(2,2)=cos(b);
return m;
};
Matrix Robot::Gl(Matrix &y){ //phi, theta, psi
  float a, b,  c;  //phi, theta, psi
    a=y.entry(3,0);b=y.entry(4,0);c=y.entry(5,0);

Matrix m(3,3);
m.entry(0,0)=sin(b)*sin(c);
m.entry(0,1)=cos(c);
m.entry(0,2)=0;

m.entry(1,0)=sin(b)*cos(c);
m.entry(1,1)=-sin(c);
m.entry(1,2)=0;

m.entry(2,0)=cos(b);
m.entry(2,1)=0;
m.entry(2,2)=1;
return m;
};
Matrix Robot::dGl(Matrix &y){
  float a,b,  c;  //phi, theta, psi
  float p,q, r;   //derivadas de phi, theta, psi respecto del tiempo
 p=y.entry(9,0); q=y.entry(10,0); r=y.entry(11,0);
a=y.entry(3,0);b=y.entry(4,0);c=y.entry(5,0);
Matrix m(3,3);

m.entry(0,0)=q*cos(b)*sin(c)+r*sin(b)*cos(c);
m.entry(0,1)=-r*sin(c);
m.entry(0,2)=0;

m.entry(1,0)=q*cos(b)*cos(c)-r*sin(b)*sin(c);
m.entry(1,1)=-r*cos(c);
m.entry(1,2)=0;

m.entry(2,0)=-q*sin(b);
m.entry(2,1)=0;
m.entry(2,2)=0;
return m;
};
Matrix Robot::mth(Matrix &y){
Matrix m(3,3);
m=Gl(y).traspose()*Il*Gl(y);
//cout<<" mt  $%%%%%%%%%%%%%%%%%"<<endl;
//m.mostrar();
return m;
}
Matrix Robot::dtheta(Matrix &y){
Matrix m(3,1);
 float p, q,  r;   //derivadas de phi, theta, psi respecto del tiempo
 p=y.entry(9,0); q=y.entry(10,0); r=y.entry(11,0);
m.entry(0,0)=p;
m.entry(1,0)=q;
m.entry(2,0)=r;
return m;
}
Matrix Robot::Qvth(Matrix y){

Matrix qv(3,1);



qv=((-1)*((Gl(y)).traspose()))*(wl(y).crossProduct(Il*wl(y))+Il*dGl(y)*dtheta(y));

return qv;
}
Matrix Robot::QeR(Matrix y){
Matrix qe(3,1);

Matrix v(3,1);
v.entry(0,0)=pow(y.entry(6,0),1);
v.entry(1,0)=pow(y.entry(7,0),1);
v.entry(2,0)=pow(y.entry(8,0),1);
qe.entry(2,0)=-mass*9.81;
qe=qe+f1+f2-100*v;
return qe;
}
Matrix Robot::Qeth(Matrix y){
    Matrix qe(3,1);


  qe=G(y).traspose()*(P(p1).crossProduct(f1));
  qe=qe+G(y).traspose()*(P(p2).crossProduct(f2));
 // qe=qe+G(y).traspose()*(P(p3).crossProduct(f3));
  //qe=qe+G(y).traspose()*(P(p4).crossProduct(f4));

Matrix v(3,1);
v.entry(0,0)=pow(y.entry(9,0),1);
v.entry(1,0)=pow(y.entry(10,0),1);
v.entry(2,0)=pow(y.entry(11,0),1);

return qe-50*v;
}

Matrix Robot::w(Matrix &y){
Matrix w_(3,1);
Matrix m(3,1);
m.entry(0,0)=y.entry(9,0);
m.entry(1,0)=y.entry(10,0);
m.entry(2,0)=y.entry(11,0);

w_=G(y)*m;
return w_;
}
Matrix Robot::wl(Matrix &y){
Matrix w_(3,1);
Matrix m(3,1);
m.entry(0,0)=y.entry(9,0);
m.entry(1,0)=y.entry(10,0);
m.entry(2,0)=y.entry(11,0);

w_=Gl(y)*m;
return w_;
}
Matrix Robot::I(Matrix &y){
Matrix m(3,3);

m=A(y)*Il*A(y).traspose();
return m;
}

Matrix Robot::M(Matrix &y){
Matrix m(6,6);

m.entry(0,0)=mass;
m.entry(1,1)=mass;
m.entry(2,2)=mass;
Matrix mt(mth(y));

//mt.mostrar();


m.entry(3,3)=mt.entry(0,0);
m.entry(4,4)=mt.entry(1,1);
m.entry(5,5)=mt.entry(2,2);
return m;
}


void Robot::InitRungeKutta(){
    ///y1=q, y2=vq;
dh=0.01; tr=0.0f; //paso de tiempo tiempo en runge kutta
///
b1.zero(3,1);b2.zero(3,1);b3.zero(3,1);b4.zero(3,1);
a1.zero(3,1);a2.zero(3,1);a3.zero(3,1);a4.zero(3,1);
///
Y.zero(12,1);
Yv.zero(12,1);
Y0.zero(12,1);
Yv0.zero(12,1);
F1.zero(12,1);
F2.zero(12,1);
F3.zero(12,1);
F4.zero(12,1);
//Dynamics

Il.zero(3,3);

/////////////////////////////////////////////////
R.zero(4,1);
phi=0;  ///configuraciones iniciales
theta=0;
psi=0;
Il.entry(0,0)=2.5063;
Il.entry(1,1)=0.906;
Il.entry(2,2)=3.4;
 mass=3.0;
k=5;

configurarTH();
///
//puntos de aplicacion
////////////////////////
///configuRACIONES INICIALES
Y=Y0; Yv=Yv0;
Y.entry(3,0)=0;
Y.entry(4,0)=0;
Y.entry(5,0)=0;
p1={.3,.5,.025};p2={-.3,.5,.025};p3={-.3,-.5,.025};p4={.3,-.5,.025}; ///locales
f1.zero(3,1);f2.zero(3,1);f3.zero(3,1);f4.zero(3,1);
u1.zero(3,1);u2.zero(3,1);u3.zero(3,1);u4.zero(3,1);
a1.entry(0,0)=0.3;a1.entry(1,0)=0.5;a1.entry(2,0)=1;
a2.entry(0,0)=-0.3;a2.entry(1,0)=0.5;a2.entry(2,0)=1;
a3.entry(0,0)=-0.3;a3.entry(1,0)=-0.5;a3.entry(2,0)=1;
a4.entry(0,0)=0.3;a4.entry(1,0)=-0.5;a4.entry(2,0)=1;

a1=1000*a1;a2=1000*a2;a3=1000*a3;a4=1000*a4;

///longitudes cuerdas
R.entry(0,0)=1000*Y.entry(0,0);
R.entry(1,0)=1000*Y.entry(1,0);
R.entry(2,0)=1000*Y.entry(2,0);
Matrix R_(3,1);
R_.entry(0,0)=R.entry(0,0);
R_.entry(1,0)=R.entry(1,0);
R_.entry(2,0)=R.entry(2,0);
b1=(R_+1000*P(p1));
b2=(R_+1000*P(p2));
b3=(R_+1000*P(p3));
b4=(R_+1000*P(p4));
le=1000*0.05;
lt1=(b1-a1).column_2norm(0);
lt2=(b2-a2).column_2norm(0);
lt3=(b3-a3).column_2norm(0);
lt4=(b4-a4).column_2norm(0);
lc1=lt1-le; ///longitud inicial de cuerda
lc2=lt2-le;
lc3=lt3-le;
lc4=lt4-le;

cout<<" lc1  "<<lc1<<endl;
cout<<" lc2  "<<lc2<<endl;
cout<<" lc3 "<<lc3<<endl;
cout<<" lc4  "<<lc4<<endl;
}
Matrix Robot::Frk(Matrix  Y_,float t_){  ///12x1

Yv.entry(0,0)=Y_.entry(6,0);
Yv.entry(1,0)=Y_.entry(7,0);
Yv.entry(2,0)=Y_.entry(8,0);
Yv.entry(3,0)=Y_.entry(9,0);
Yv.entry(4,0)=Y_.entry(10,0);
Yv.entry(5,0)=Y_.entry(11,0);
Matrix A_(Grk(Y_,t_));
//cout<<" A_"<<endl;
//A_.mostrar();
Yv.entry(6,0)=A_.entry(0,0);
Yv.entry(7,0)=A_.entry(1,0);
Yv.entry(8,0)=A_.entry(2,0);
Yv.entry(9,0)=A_.entry(3,0);
Yv.entry(10,0)=A_.entry(4,0);
Yv.entry(11,0)=A_.entry(5,0);

return Yv;
}
Matrix Robot::Grk(Matrix Y_,float t_){ ///6x1
///la aceleración depende de la posición y posiblemente de la velocidad

Matrix aq(6,1);
Matrix qe(6,1);
Matrix qv(6,1);
lt1=(b1-a1).column_2norm(0);
lt2=(b2-a2).column_2norm(0);
lt3=(b3-a3).column_2norm(0);
lt4=(b4-a4).column_2norm(0);
cout<<" lt1 =" <<lt1-(lc1+le)<<endl;
u1=(a1-b1).unit(0);
u2=(a2-b2).unit(0);
u3=(a3-b3).unit(0);
u4=(a4-b4).unit(0);
f1=k*(lt1-(lc1+le))*u1;
f2=k*(lt2-(lc2+le))*u2;
f3=k*(lt3-(lc3+le))*u3;
f4=k*(lt4-(lc4+le))*u4;



qe.entry(0,0)=QeR(Y_).entry(0,0);
qe.entry(1,0)=QeR(Y_).entry(1,0);
qe.entry(2,0)=QeR(Y_).entry(2,0);
qe.entry(3,0)=Qeth(Y_).entry(0,0);
qe.entry(4,0)=Qeth(Y_).entry(1,0);
qe.entry(5,0)=Qeth(Y_).entry(2,0);

qv.entry(0,0)=0;
qv.entry(1,0)=0;
qv.entry(2,0)=0;
qv.entry(3,0)=Qvth(Y_).entry(0,0);
qv.entry(4,0)=Qvth(Y_).entry(1,0);
qv.entry(5,0)=Qvth(Y_).entry(2,0);
//cout<<"qv"<<endl;
//M(Y_).inverse1().mostrar();

aq=(M(Y_).inverse1())*(qe+qv);
return aq;
}
void Robot::RungeKuttaIntegrar(){

F1=dh*Frk(Y,tr);
F2=dh*Frk(Y+0.5f*F1,tr+0.5f*dh);
F3=dh*Frk(Y+0.5f*F2,tr+0.5f*dh);
F4=dh*Frk(Y+F3,tr+dh);
Y=Y+ 0.166666666666666666666*(F1+2.0f*F2+2.0f*F3+F4);

phi=Y.entry(3,0);
theta=Y.entry(4,0);
psi=Y.entry(5,0);
AplicarTHz(phi*180/PI,{0,0,0.0});
THList[0]=(THz);
AplicarTHx(theta*180/PI,{0,0,0.0}); //b2
THList[1]=(THx);
AplicarTHz(psi*180/PI,{0,0,0.0}); //b3
THList[2]=(THz);
R.entry(0,0)=1000*Y.entry(0,0);
R.entry(1,0)=1000*Y.entry(1,0);
R.entry(2,0)=1000*Y.entry(2,0);
///graficar el movimiento del eslabón3D
Matrix R_(3,1);
R_.entry(0,0)=R.entry(0,0);
R_.entry(1,0)=R.entry(1,0);
R_.entry(2,0)=R.entry(2,0);
b1=(R_+1000*P(p1));
b2=(R_+1000*P(p2));
b3=(R_+1000*P(p3));
b4=(R_+1000*P(p4));
}
void Robot::inicializar(){

placa=new modelo3D();
placa->leer("placa.STL");
modelos.push_back(placa);

InitRungeKutta();


}
void Robot::configurarTH(){

AplicarTHz(phi*180/PI,{0,0,0.0}); //base
THList.push_back(THz);
AplicarTHx(theta*180/PI,{0,0,0.0}); //b2
THList.push_back(THx);
AplicarTHz(psi*180/PI,{0,0,0.0}); //b3
THList.push_back(THz);

}


void Robot::renderizar(){
RungeKuttaIntegrar();
glBegin(GL_LINES);
glVertex3f(a1.entry(0,0),a1.entry(1,0),a1.entry(2,0));
glVertex3f(b1.entry(0,0),b1.entry(1,0),b1.entry(2,0));

glVertex3f(a2.entry(0,0),a2.entry(1,0),a2.entry(2,0));
glVertex3f(b2.entry(0,0),b2.entry(1,0),b2.entry(2,0));
/*
glVertex3f(a3.entry(0,0),a3.entry(1,0),a3.entry(2,0));
glVertex3f(b3.entry(0,0),b3.entry(1,0),b3.entry(2,0));

glVertex3f(a4.entry(0,0),a4.entry(1,0),a4.entry(2,0));
glVertex3f(b4.entry(0,0),b4.entry(1,0),b4.entry(2,0));
*/
glEnd();

tr=tr+dh;
TH.resetIdentity();
modelo3D *model;

for (int m=0;m<modelos.size();m++){
    model=modelos[m];
    TH= THList[m]*THList[m+1]*THList[m+2];


vector3d ux,uy,uz,O;
ux={1,0,0};
uy={0,1,0};
uz={0,0,1};

Matrix ux4(ux,1),uy4(uy,1),uz4(uz,1),O4(O,1);


ux4=TH*ux4-TH*O4;
uy4=TH*uy4-TH*O4;
uz4=TH*uz4-TH*O4;
O4=TH*O4+R;


ux={ux4.aij[0][0],ux4.aij[1][0],ux4.aij[2][0]};
uy={uy4.aij[0][0],uy4.aij[1][0],uy4.aij[2][0]};
uz={uz4.aij[0][0],uz4.aij[1][0],uz4.aij[2][0]};
O={O4.aij[0][0],O4.aij[1][0],O4.aij[2][0]};

//if (m<2){
         Drawarrow3D(O,O+100*ux,{1,0.1,0.2},1,10);
         Drawarrow3D(O,O+100*uy,{.1,1,0.2},1,10);
         Drawarrow3D(O,O+100*uz,{0.1,0.2,1},1,10);
       //  }
         glColor4f(fabs(cos(m*PI/modelos.size())),fabs(sin(20*(m-5)*PI/modelos.size())),0.2,0.5);

glEnable(GL_BLEND);
 glBegin(GL_TRIANGLES);

  glFrontFace(GL_FRONT_AND_BACK);
    for (int i=0;i<model->ntriangles;i++){

vector3d v1=model->triangulos[i].vertices[0];   //posiciones locales
vector3d v2=model->triangulos[i].vertices[1];
vector3d v3=model->triangulos[i].vertices[2];
Matrix v14(v1,1),v24(v2,1),v34(v3,1);

v14=TH*v14+R;
v24=TH*v24+R;
v34=TH*v34+R;
v1={v14.entry(0,0),v14.entry(1,0),v14.entry(2,0)};
v2={v24.entry(0,0),v24.entry(1,0),v24.entry(2,0)};
v3={v34.entry(0,0),v34.entry(1,0),v34.entry(2,0)};



Matrix N(4,1),d14(4,1),d24(4,1);
d14=v24-v14;
d24=v34-v14;
vector3d d1,d2,n;
d1={d14.entry(0,0),d14.entry(1,0),d14.entry(2,0)};
d2={d24.entry(0,0),d24.entry(1,0),d24.entry(2,0)};
n=d1*d2;  ///devuelve el producto vectorial
n.normalize();



        glNormal3f(n.x,n.y,n.z);
        glVertex3f(v1.x,v1.y,v1.z);
        glVertex3f(v2.x,v2.y,v2.z);
        glVertex3f(v3.x,v3.y,v3.z);
    }
glEnd();
// }
 glDisable(GL_BLEND);


///DIBUJAR EJES

//}
}

}

void Robot::DefinirTHx(float dtheta, vector3d d){

THx.aij[0][0]=1;
THx.aij[0][1]=0;
THx.aij[0][2]=0;
THx.aij[0][3]=d.x;

THx.aij[1][0]=0;
THx.aij[1][1]=cos(dtheta);
THx.aij[1][2]=-sin(dtheta);
THx.aij[1][3]=d.y;

THx.aij[2][0]=0;
THx.aij[2][1]=sin(dtheta);
THx.aij[2][2]=cos(dtheta);
THx.aij[2][3]=d.z;

THx.aij[3][0]=0;
THx.aij[3][1]=0;
THx.aij[3][2]=0;
THx.aij[3][3]=1;

}
void Robot::DefinirTHy(float dtheta, vector3d d){


THy.aij[0][0]=cos(dtheta);
THy.aij[0][1]=0;
THy.aij[0][2]=sin(dtheta);

THy.aij[1][0]=0;
THy.aij[1][1]=1;
THy.aij[1][2]=0;

THy.aij[2][0]=-sin(dtheta);
THy.aij[2][1]=0;
THy.aij[2][2]=cos(dtheta);

THy.aij[3][0]=0;
THy.aij[3][1]=0;
THy.aij[3][2]=0;
THy.aij[3][3]=1;

THy.aij[0][3]=d.x;
THy.aij[1][3]=d.y;
THy.aij[2][3]=d.z;
}
void Robot::DefinirTHz(float dtheta, vector3d d){

THz.aij[0][0]=cos(dtheta);
THz.aij[0][1]=-sin(dtheta);
THz.aij[0][2]=0;
THz.aij[0][3]=d.x;

THz.aij[1][0]=sin(dtheta);
THz.aij[1][1]=cos(dtheta);
THz.aij[1][2]=0;
THz.aij[1][3]=d.y;

THz.aij[2][0]=0;
THz.aij[2][1]=0;
THz.aij[2][2]=1;
THz.aij[2][3]=d.z;

THz.aij[3][0]=0;
THz.aij[3][1]=0;
THz.aij[3][2]=0;
THz.aij[3][3]=1;

}

void Robot::AplicarTHx(float theta, vector3d d){
theta=theta*PI/180.0;

DefinirTHx(theta,d);

}
void Robot::AplicarTHy(float theta, vector3d d){
theta=theta*PI/180.0;
DefinirTHy(theta,d);

}
void Robot::AplicarTHz(float theta, vector3d d){
theta=theta*PI/180.0;
DefinirTHz(theta,d);

}

void Robot::Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R)
{

double color1,color2,color3,a,b,c,d,e;



color1=color.x;//abs(color1/255);
color2=color.y;//abs(color2/255);
color3=color.z;//abs(color3/255);

glColor3f( color1,color2, color3);

vector3d n=B-A,np,vertex[10],normallight;
n.normalize();
if(n.z!=0)np={1,1,(-1/n.z)*(n.x+n.y)};
else if(n.y!=0)np={1,(-1/n.y)*(n.x+n.z),1};
else np={(-1/n.x)*(n.y+n.z),1,1};

np.normalize();
vertex[0]=R*np;
vertex[2]=R*(n*np).normalize();
vertex[1]=R*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=R*(n*vertex[2]).normalize();
vertex[3]=R*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=R*(n*vertex[4]).normalize();
vertex[5]=R*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=R*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
int nx=8;
double d_thetha,fraccion=0.1,radioflecha=R+.7*R;
d_thetha=2.0f*PI/nx;


  ///tubos
 glBegin( GL_TRIANGLE_STRIP );

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+A.x,vertex[i].y+A.y,vertex[i].z+A.z);

glVertex3f(vertex[i].x+B.x-fraccion*(B.x-A.x),vertex[i].y+B.y-fraccion*(B.y-A.y),vertex[i].z+B.z-fraccion*(B.z-A.z));

    // top face

                }

glEnd();



//flecha
vertex[0]=radioflecha*np;
vertex[2]=radioflecha*(n*np).normalize();
vertex[1]=radioflecha*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=radioflecha*(n*vertex[2]).normalize();
vertex[3]=radioflecha*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=radioflecha*(n*vertex[4]).normalize();
vertex[5]=radioflecha*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=radioflecha*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
vector3d Ap(B-fraccion*(B-A));



 glBegin( GL_TRIANGLE_STRIP );  //flecha

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+Ap.x,vertex[i].y+Ap.y,vertex[i].z+Ap.z);


glNormal3f(n.x, n.y, n.z);
glVertex3f(Ap.x+fraccion*(B-A).x,Ap.y+fraccion*(B-A).y,Ap.z+fraccion*(B-A).z);

    // top face

                }

glEnd();


}
