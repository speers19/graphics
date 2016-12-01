#include <D3d_matrix.h>



//NOTE: Does not contain negate_x, negate_y, negate_z

/////////////////////////////////////////////////////////
//////////VECTOR OPERATIONS//////////////////////////////
/////////////////////////////////////////////////////////


int D3d_x_product (double res[3], double a[3], double b[3]) 
//res = a x b, cross product
//make it SAFE
{
  double tmpX = ( (a[1] * b[2]) - (a[2] * b[1]) ) ;
  double tmpY = - ( (a[0] * b[2]) - (a[2] * b[0]) ) ;
  double tmpZ = ( (a[0] * b[1]) - (a[1] * b[0]) ) ;

  res[0] = tmpX ;
  res[1] = tmpY ;
  res[2] = tmpZ ;

  return 1;
}

double D3d_dot_product (double a[3], double b[3]) 
//returns a double 'result' from dot product
{
  double result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] ;
  return result;
}

double D3d_vector_distance (double a[3]) 
// returns a double 'result' from distance formula
{
  double result = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] ) ;
  return result;
}

/////////////////////////////////////////////////////////
//////////BASIC MATRIX OPERATIONS////////////////////////
/////////////////////////////////////////////////////////

int D3d_print_mat (double a[4][4])
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
    for (c = 0 ; c < 4 ; c++ ) {
     printf(" %12.4lf ",a[r][c]) ;
   }
   printf("\n") ;
 }

 return 1 ;
}




int D3d_copy_mat (double a[4][4], double b[4][4])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
    for (c = 0 ; c < 4 ; c++ ) {
     a[r][c] = b[r][c] ;
   }
 }

 return 1 ;
} 



int D3d_mat_mult (double res[4][4], double a[4][4], double b[4][4])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// D3d_mat_mult(p,  p,q) or D3d_mat_mult(p,  q,p) or  D3d_mat_mult(p, p,p)
{
  double tmp[4][4];
  int r, c ;
  for (r = 0 ; r < 4; r++) {
    for (c = 0 ; c < 4 ; c++) {
      tmp[r][c] = (a[r][0]*b[0][c] + a[r][1]*b[1][c] + a[r][2]*b[2][c] + a[r][3]*b[3][c]);
    }
  }

  D3d_copy_mat(res, tmp);

  return 1;
}


int D3d_make_identity (double a[4][4])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
    for (c = 0 ; c < 4 ; c++ ) {
     if (r == c) a[r][c] = 1.0 ;
     else    a[r][c] = 0.0 ;
   }
 }

 return 1 ;
} 




/////////////////////////////////////////////////////////
//////////TRANSLATE SCALE AND ROTATE/////////////////////
/////////////////////////////////////////////////////////


int D3d_translate (double a[4][4], double b[4][4], double dx, double dy, double dz)
// a = translation*a  
// b = b*translation_inverse  
{
  double t[4][4] ;

  D3d_make_identity(t) ;

  t[0][3] =  dx ;  t[1][3] = dy ; t[2][3] = dz;  
  D3d_mat_mult(a,  t,a) ;

  t[0][3] =  -dx ;  t[1][3] = -dy ; t[2][3] = -dz;  
  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}

int D3d_scale (double a[4][4], double b[4][4], double sx, double sy, double sz)
// a = scale*a  
// b = b*scale_inverse  
{
  double t[4][4] ; 

  D3d_make_identity(t) ;

  t[0][0] = sx ; t[1][1] = sy; t[2][2] = sz;
  D3d_mat_mult(a, t,a) ;

  t[0][0] = (1/sx) ; t[1][1] = (1/sy); t[2][2] = (1/sz);
  D3d_mat_mult(b, b,t) ;

  return 1;
}

int D3d_rotate_x (double a[4][4], double b[4][4], double radians)
// a = rotate*a  
// b = b*rotate_inverse  
{
  double t[4][4] ; 

  D3d_make_identity(t) ;

  t[1][1] = cos(radians) ; t[1][2] = -sin(radians) ;
  t[2][1] = sin(radians) ; t[2][2] = cos(radians) ;
  D3d_mat_mult(a, t,a) ;

  t[1][1] = cos(-radians) ; t[1][2] = -sin(-radians) ;
  t[2][1] = sin(-radians) ; t[2][2] = cos(-radians) ;
  D3d_mat_mult(b, b,t) ;

  return 1;
}

int D3d_rotate_y (double a[4][4], double b[4][4], double radians)
// a = rotate*a  
// b = b*rotate_inverse  
{
  double t[4][4] ; 

  D3d_make_identity(t) ;

  t[0][0] = cos(radians) ; t[0][2] = sin(radians) ;
  t[2][0] = -sin(radians) ; t[2][2] = cos(radians) ;
  D3d_mat_mult(a, t,a) ;

  t[0][0] = cos(-radians) ; t[0][2] = sin(-radians) ;
  t[2][0] = -sin(-radians) ; t[2][2] = cos(-radians) ;
  D3d_mat_mult(b, b,t) ;

  return 1;
}

int D3d_rotate_z (double a[4][4], double b[4][4], double radians)
// a = rotate*a  
// b = b*rotate_inverse  
{
  double t[4][4] ; 

  D3d_make_identity(t) ;

  t[0][0] = cos(radians) ; t[0][1] = -sin(radians) ;
  t[1][0] = sin(radians) ; t[1][1] = cos(radians) ;
  D3d_mat_mult(a, t,a) ;

  t[0][0] = cos(-radians) ; t[0][1] = -sin(-radians) ;
  t[1][0] = sin(-radians) ; t[1][1] = cos(-radians) ;
  D3d_mat_mult(b, b,t) ;

  return 1;
}

int D3d_mat_mult_points (double *X, double *Y, double *Z,
                         double m[4][4],
                         double *x, double *y, double *z, int numpoints) {

  int i;
  double tmp;

  for (i = 0; i < numpoints; i++) {

    double tmpX = (m[0][0] * x[i]) + (m[0][1] * y[i]) + (m[0][2] * z[i]) + m[0][3];
    double tmpY = (m[1][0] * x[i]) + (m[1][1] * y[i]) + (m[1][2] * z[i]) + m[1][3];
    Z[i] = (m[2][0] * x[i]) + (m[2][1] * y[i]) + (m[2][2] * z[i]) + m[2][3];
    X[i] = tmpX;
    Y[i] = tmpY;

  }


return 1;
}