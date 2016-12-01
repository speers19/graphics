#include <D2d_matrix.h>



/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( 1 )          (1)

instead of (x',y',1) = (x,y,1) * M  

*/



int D2d_print_mat (double a[3][3])
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
    for (c = 0 ; c < 3 ; c++ ) {
     printf(" %12.4lf ",a[r][c]) ;
   }
   printf("\n") ;
 }

 return 1 ;
} 




int D2d_copy_mat (double a[3][3], double b[3][3])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
    for (c = 0 ; c < 3 ; c++ ) {
     a[r][c] = b[r][c] ;
   }
 }

 return 1 ;
} 



int D2d_mat_mult (double res[3][3], double a[3][3], double b[3][3]) 
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// D2d_mat_mult(p,  p,q) or D2d_mat_mult(p,  q,p) or  D2d_mat_mult(p, p,p)
{
  double tmp[3][3];
  int r, c ;
  for (r = 0 ; r < 3; r++) {
    for (c = 0 ; c < 3 ; c++) {
      tmp[r][c] = (a[r][0]*b[0][c] + a[r][1]*b[1][c] + a[r][2]*b[2][c]);
    }
  }

  D2d_copy_mat(res, tmp);

  return 1;
}


int D2d_make_identity (double a[3][3])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
    for (c = 0 ; c < 3 ; c++ ) {
     if (r == c) a[r][c] = 1.0 ;
     else    a[r][c] = 0.0 ;
   }
 }

 return 1 ;
} 




/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


int D2d_translate (double a[3][3], double b[3][3], double dx, double dy)
// a = translation*a  
// b = b*translation_inverse  
{
  double t[3][3] ;

  D2d_make_identity(t) ;

  t[0][2] =  dx ;  t[1][2] = dy ;  
  D2d_mat_mult(a,  t,a) ;

  t[0][2] = -dx ;  t[1][2] = -dy ;
  D2d_mat_mult(b,  b,t) ;

  return 1 ;
}

int D2d_scale (double a[3][3], double b[3][3], double sx, double sy)
// a = scale*a
// b = b*scale_inverse
{
  double t[3][3] ; 

  D2d_make_identity(t) ;

  t[0][0] = sx ; t[1][1] = sy;
  D2d_mat_mult(a, t,a) ;

  t[0][0] = (1/sx) ; t[1][1] = (1/sy);
  D2d_mat_mult(b, b,t) ;

  return 1;
}

int D2d_rotate (double a[3][3], double b[3][3], double radians)
// a = rotate*a
// b = b*rotate_inverse
{
  double t[3][3] ; 

  D2d_make_identity(t) ;

  t[0][0] = cos(radians) ; t[0][1] = -sin(radians) ;
  t[1][0] = sin(radians) ; t[1][1] = cos(radians) ;
  D2d_mat_mult(a, t,a) ;

  t[0][0] = cos(-radians) ; t[0][1] = -sin(-radians) ;
  t[1][0] = sin(-radians) ; t[1][1] = cos(-radians) ;
  D2d_mat_mult(b, b,t) ;

  return 1;
}

int D2d_negate_x (double a[3][3], double b[3][3]) 
// a = reflect*a
// b = b*reflect_inverse
{  
  double t[3][3] ;

  D2d_make_identity(t) ;

  t[0][0] = -1; t[0][1] = -1; t[0][2] = -1;
  D2d_mat_mult(a, t,a);

  t[0][0] = 1; t[0][1] = 1; t[0][2] = 1;
  D2d_mat_mult(b, b,t);

  return 1;
}

int D2d_negate_y (double a[3][3], double b[3][3]) 
// a = reflect*a
// b = b*reflect_inverse
{  
  double t[3][3] ;

  D2d_make_identity(t) ;

  t[1][0] = -1; t[1][1] = -1; t[1][2] = -1;
  D2d_mat_mult(a, t,a);

  t[1][0] = 1; t[1][1] = 1; t[1][2] = 1;
  D2d_mat_mult(b, b,t);

  return 1;      
}

int D2d_mat_mult_points (double *X, double *Y,
 double m[3][3],
 double *x, double *y, int numpoints)
{

  int i;
  double tmp;

  for (i = 0; i < numpoints; i++) {

    tmp = (m[0][0] * x[i]) + (m[0][1] * y[i]) + m[0][2];
    Y[i] = (m[1][0] * x[i]) + (m[1][1] * y[i]) + m[1][2];
    X[i] = tmp;

  }


return 1;
}


