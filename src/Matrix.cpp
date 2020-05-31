#include "Matrix.h"
using namespace std;
Matrix Matrix::unit(int col){
Matrix u(n,1);
float norm = column_2norm(col);
	for (int i = 0; i<n; i++) {
		 u.entry(i, 0)=entry(i, col)/norm;
	}
	return u;
}
Matrix Matrix::crossProduct(Matrix A){

float a,b,c,d,e,f;
a=entry(0,0);
b=entry(1,0);
c=entry(2,0);
d=A.entry(0,0);
e=A.entry(1,0);
f=A.entry(2,0);
Matrix p(3,1);
p.entry(0,0)=b*f-c*e;
p.entry(1,0)=c*d-a*f;
p.entry(2,0)=a*e-b*d;
return p;
}
Matrix Matrix::traspose(){
    Matrix T(m,n);
    for( int i=0;i<n;i++)
        for( int j=0;j<m;j++)
{
    T.entry(j,i)=entry(i,j);
}
return T;
}
Matrix::Matrix(vector3d v){
int nn=3; int mm=1;

if (nn>0&&mm>0){
  aij=new double*[nn];
  for(int i = 0; i < nn; i++)  aij[i] = new double[mm];
  n=nn;
  m=mm;

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
    //    coutincorrect <<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
           // cin>>ij(i,j);
aij[i][j]=0;
    };
}
ij(0,0)=v.x;
ij(1,0)=v.y;
ij(2,0)=v.z;
}
Matrix::Matrix(vector3d v,float w){
int nn=4; int mm=1;

if (nn>0&&mm>0){
  aij=new double*[nn];
  for(int i = 0; i < nn; i++)  aij[i] = new double[mm];
  n=nn;
  m=mm;

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
    //    coutincorrect <<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
           // cin>>ij(i,j);
aij[i][j]=0;
    };
}

ij(0,0)=v.x;
ij(1,0)=v.y;
ij(2,0)=v.z;
ij(3,0)=1;
}
double Matrix::cofactor (int i, int j){
double c;
c=pow (-1,i+j)*Mij(i,j).determinante();  //cofactor
return c;
};
Matrix Matrix::inversa(){
double d=(*this).determinante();

if (d==0) { Matrix B(1,1); cout<<" no existe la inversa"<<endl; return B;}
 if(n>2) {
    Matrix B(n,n);
     for( int i = 0; i < n; i++)
    for(int j = 0; j < n; j++){
         B.aij[j][i]=(1.0/d)*cofactor(i,j);


    }
    return B;

}
else if(n==2) {
    Matrix B(n,n);
     B.aij[0][0]=(1.0/d)*aij[1][1];

      B.aij[0][1]=-(1.0/d)*aij[0][1];
       B.aij[1][0]=-(1.0/d)*aij[1][0];
      B.aij[1][1]=(1.0/d)*aij[0][0];
    return B;

}
else if(n==1) {
    Matrix B(1,1);
     B.aij[0][0]=1/aij[0][0];


    return B;

}

};
float Matrix::column_2norm(int col){
	float norm = 0.0L;
	for (int i = 0; i<n; i++) {
		norm = norm + entry(i, col)*entry(i, col);
	}
	return sqrt(norm);
};
////////////
void Matrix::Crout_decomposition(Matrix &L, Matrix &U) {
	//A must be square

	for (int i = 0; i<n; i++) U.entry(i, i) = 1.00;
	for (int i = 0; i<n; i++) { L.entry(i, 0) = entry(i, 0); }

	for (int j = 1; j<n; j++) { U.entry(0, j) = entry(0, j) / L.entry(0, 0); }

	for (int j = 1; j<(n - 1); j++) {

		for (int i = j; i<n; i++) {
		double S = 0;
			for (int k = 0; k <= (j - 1); k++)S = S + L.entry(i, k)*U.entry(k, j);
			L.entry(i, j) = entry(i, j) - S;
		};

		for (int k = (j + 1); k<n; k++) {
			double S1 = 0;
			for (int i = 0; i <= (j - 1); i++)S1 = S1 + L.entry(j, i)*U.entry(i, k);
			U.entry(j, k) = (entry(j, k) - S1) / (L.entry(j, j));
		}

	}



	double S2 = 0;
	for (int k = 0; k <= (n - 2); k++)S2 = S2 + L.entry(n - 1, k)*U.entry(k, n - 1);

	L.entry(n - 1, n - 1) = entry(n - 1, n - 1) - S2;




};
 Matrix Matrix::inverse1() {
//correcto


	//A must be square, uses LU
	Matrix L(n, n), U(n, n), Y(n, n), Z(n, n);  ///Y is  L inverse, Z is U inverse
	Crout_decomposition(L, U);
	for (int j = 0; j<n; j++) {
		Y.entry(j, j) = 1.0L / L.entry(j, j);

		for (int i = j + 1; i<n; i++) {
			double S = 0;
			for (int k = j; k <= (i - 1); k++) { S = S + L.entry(i, k)*Y.entry(k, j); }
			Y.entry(i, j) = -S / L.entry(i, i);
		}


	}

	for (int j = 0; j<n; j++) {
		Z.entry(j, j) = 1.0L / U.entry(j, j);

		for (int i = (j - 1); i>-1; i--) {
			double S = 0;
			for (int k = (i + 1); k <= j; k++) { S = S + U.entry(i, k)*Z.entry(k, j); }
			Z.entry(i, j) = -S / U.entry(i, i);
		}

	}

	return (Z*Y);

	}
/////////////

    int Matrix::size () const{return n*m;};
    Matrix Matrix::Mij(int a,int b){ //menor de una matriz
        if (n==m&&size()>0){ // eliminar fila a y columna b

    int s=n-1;
    Matrix B(s,s); //matriz de tamaño reducido
    for (int i=0; i<a;i++)
    for (int j=0; j<b;j++)

    {
      B.aij[i][j]=aij[i][j]; //primera parte
    }

      for (int i=0; i<a;i++)
    for (int j=b+1; j<m;j++)

    {
      B.aij[i][j-1]=aij[i][j];  //se recorren las columnas
    }

        for (int i=a+1; i<n;i++)
    for (int j=0; j<b;j++)

    {
      B.aij[i-1][j]=aij[i][j];  //se recorren las columnas
    }


        for (int i=a+1; i<n;i++)
    for (int j=b+1; j<m;j++)

    {
      B.aij[i-1][j-1]=aij[i][j];   //se recorren las columnas
    }
    return B;
}
};
vector3d Matrix::operator *(const vector3d &P){//operacion multiplicacion por punto
Matrix A(3,1);
A.aij[0][0]=P.x;
A.aij[1][0]=P.y;
A.aij[2][0]=P.z;

Matrix T((*this)*A);

vector3d Pr;
Pr.x=T.aij[0][0];
Pr.y=T.aij[1][0];
Pr.z=T.aij[2][0];
return Pr;
}
Matrix::Matrix(int nn,int mm){ //constructor a partir de n y m
if (nn>0&&mm>0){
  aij=new double*[nn];
  for(int i = 0; i < nn; i++)  aij[i] = new double[mm];
  n=nn;
  m=mm;

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
    //    coutincorrect <<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
           // cin>>ij(i,j);
aij[i][j]=0;
    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};
Matrix::Matrix(){



n=m=0;
  //numero de filas

};
void Matrix::zero(int nn, int mm){
     n=nn;
        m=mm;
     //constructor de mat


if (n>0&&m>0){



  aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];


  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){

        ij(i,j)=0.0;;

    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};
 void Matrix::resetIdentity(){

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
if(  i==j)ij(i,j)=1.0;
else ij(i,j)=0.0;

    };
 }

void Matrix::identity(int nn){
     n=nn;
        m=nn;
     //constructor de mat


if (n>0&&m>0){



  aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];


  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
if(  i==j)ij(i,j)=1.0;
else ij(i,j)=0.0;

    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};

void Matrix::definir(){
    cout <<"Creando Matriz "<<endl;
     cout <<"introduce el numero de filas "<<endl;
     cin>> n;
         cout <<"introduce  el numero de columnas "<<endl;
              cin>> m;
     //constructor de mat


if (n>0&&m>0){



  aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];


  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
        cout<<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
        cin>>ij(i,j);

    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};
void Matrix::mostrar(){
  if (this->n>0&&this->m>0){


  cout<<"\n"<<endl;
for( int i = 0; i < n; i++){
    for( int j = 0; j < m; j++){

      cout<<ij(i,j)<<"  "<<flush;
    }
    cout<<"\n"<<endl;
    }
}
}
double& Matrix::ij(int i,int j){ //devuelve la direccion
if (i>=0&&j>=0&&i<n&&j<m) return aij[i][j]; //c++ numera las ijs de una matrix a partir de cero, nosotros a partir de 1

};
Matrix::Matrix(const Matrix &B){
if(B.size()!=0){
        n=B.n;
        m=B.m;


        aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];

    for( int i = 0; i < n; i++)
    for( int j = 0; j < m; j++){

       aij[i][j]=B.aij[i][j];

    }

}
};
Matrix Matrix::operator +( const Matrix &A){
if (n>0&&m>0){
 if (n==A.n&&m==A.m){

Matrix B(n,m);

for( int i = 0; i < n; i++)
    for( int j = 0; j < m; j++){

        B.aij[i][j]= aij[i][j]+A.aij[i][j];

    }

      return B;
    }
    else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return B;};

    }


};
Matrix& Matrix::operator =( const Matrix &A){
if (A.n>0&&A.m>0&&n==A.n&&m==A.m){

Matrix B(n,m);
for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){

       aij[i][j]=A.aij[i][j];
    }
      return *this;

     }   else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return *this;}

};
Matrix Matrix::operator -( Matrix A){
if (this->n>0&&this->m>0){
 if (this->n==A.n&&this->m==A.m){
Matrix B(n,m);
for( int i = 0; i < this->n; i++)
    for(int j = 0; j < this->m; j++){

        B.ij(i,j)= ij(i,j)-A.ij(i,j);
    }
      return B;
    }    else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return B;}
}

};
Matrix Matrix::operator *( Matrix B){ //multiplicaciion de matrices
if (size()>0&&B.size()>0){
 if (m==B.n){
Matrix A(n,B.m);

for( int i = 0; i < n; i++)
    for(int j = 0; j < B.m; j++){
            float s=0;
    for(int k=0;k<m;k++){s=s+ij(i,k)*B.ij(k,j); }
        A.ij(i,j)= s;
    }
      return A;
    }    else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return B;} }


};
Matrix operator*(double tt, Matrix A){
if (A.n>0&&A.m>0){

Matrix B(A.n,A.m);
for(int i = 0; i < A.n; i++)
    for( int j = 0; j < A.m; j++){

        B.ij(i,j)= tt*A.ij(i,j);
    }


    return B;
     }


};
Matrix::~Matrix(){


			if (n>0){
				for (int i = 0; i<n; i++){
					delete[] aij[i];  //se elimnan las columnas
				}
				delete[] aij; //eliminar las filas
				aij = NULL;
			}
	}

double Matrix::determinante(){  //se define en forma recursiva
double d=0;
if (n!=m) {cout<<"la matriz no es cuadrada"<<endl; return 0; }
 if (n>2){
for (int j=0; j<m;j++){

    d=d+aij[0][j]*pow (-1,j)*Mij(0,j).determinante();

}
return d;
} else if (n==2) {
    d=aij[0][0]*aij[1][1]-aij[0][1]*aij[1][0];
    return d;

} else if(n==1){

return aij[0][0];
};

}

double& Matrix::entry(int i, int j){
return aij[i][j];
}
