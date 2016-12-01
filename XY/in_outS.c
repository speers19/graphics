#include <FPT.h>

int swidth, sheight ;

int click_and_save (double *x, double *y)
{
  int n ;
  double P[2] ;

  G_rgb(0,1,0.5) ;
  G_fill_rectangle(0,0,swidth,20) ;

  G_rgb(1,0,0) ;
  G_wait_click(P);

  n = 0 ;
  while (P[1] > 20) {
    x[n] = P[0] ;
    y[n] = P[1] ;
    G_circle(x[n],y[n],2) ;
    if (n > 0) { G_line(x[n-1],y[n-1], x[n],y[n]) ;}
    n++ ;
    G_wait_click(P) ;
  }

  return n ;
}



int in_out (double *x, double *y, int n, double P[2])
// return 1 if point P is inside the convex polygon
// else return 0
{

  int i;
  double m;
  double b;
  int answer = 0; 

  double cx = 0;
  double cy = 0; 

  for (i = 0; i < n; i++) {
    cx = cx + x[i];
    cy = cy + y[i];
  }

  cx = cx/n;
  cy = cy/n;

  for (i = 0; i < n; i++) {

    if (i < n) {
      m = ( y[i] - y[i+1] ) / ( x[i] - x[i+1] ) ;
      b = y[i] - ( m * x[i] ) ;
    }

    if (i == (n-1)) {
      m = ( y[i] - y[0] ) / ( x[i] - x[0] ) ;
      b = y[i] - ( m * x[i] ) ;
    } 

    printf("%lf, %lf \n", m, b);
    
    if ( ( m*P[0] + b - P[1] ) > 0 && ( m*cx + b - cy ) > 0 ) {
      answer += 1;
    }

    if ( ( m*P[0] + b - P[1] ) < 0 && ( m*cx + b - cy ) < 0 ) {
      answer += 1; 
    }

     printf("%d    %d \n,", answer, i);

  }

  if ((answer / n) == 1) {
    return 1;
  }
  else {
    return 0;
  }


}

int main()
{
  double xp[1000],yp[1000] ;
  int n,q ;
  double P[2] ;


  swidth = 700 ; sheight = 700 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  G_rgb(1,0,0) ;
  n = click_and_save(xp,yp) ;
  G_rgb(0,1,0) ;
  G_fill_polygon(xp,yp,n) ;

  while (q != 'q') {

  G_wait_click(P) ;
  
  int color = in_out(xp, yp, n, P);

  if (color == 1) {
    G_rgb(0,0,1);
    G_circle(P[0], P[1], 2);
  }
  else {
    G_rgb(1,1,0);
    G_circle(P[0], P[1], 2);
  }

  q = G_wait_key() ;
  }
}