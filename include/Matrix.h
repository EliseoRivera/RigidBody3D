#ifndef MATRIX_H
#define MATRIX_H
#include"vector3d.h"
#include<iostream>
class Matrix
{
public:
   double **aij;
   int n;
   int m;
   Matrix(vector3d v);
   Matrix (vector3d v, float w);
   double cofactor (int i, int j);
   Matrix inversa();
   Matrix crossProduct(Matrix A);
   Matrix inverse1();
   void Crout_decomposition(Matrix &L, Matrix &U);
     float column_2norm(int col);
     Matrix unit(int col);
   int size () const;
   Matrix Mij(int a,int b);
   vector3d operator *(const vector3d &P);
   Matrix(int nn,int mm);
   Matrix();
   void zero(int nn, int mm);
   void identity(int nn);
   void resetIdentity();
   void definir();
   void mostrar();
   double& ij(int i,int j);
   Matrix(const Matrix &B);
   Matrix operator +( const Matrix &A);
   Matrix& operator =( const Matrix &A);
   Matrix operator -( Matrix A);
   Matrix operator *( Matrix B);
   double& entry(int i, int j);
   Matrix traspose();
   friend  Matrix operator*(double tt, Matrix A);
  ~Matrix();
   double determinante();
};

#endif // MATRIX_H
