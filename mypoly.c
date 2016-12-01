#include <FPT.h>

double calculate_x(double y0, double x0, double y1, double x1, double y) {
  double m = (y1 - y0) / ((x1 - x0) + .01);
  double b = y0 - (m * x0);

  double x = (y - b) / m;
  return x;
}

void sort(double *x, int n) {

  int s, k, i;
  double t; //temp 

  for (k = 0; k < n; k++) {

    s = k; //location of smallest value

    for (i = k+1; i < n; i++) {
      if (x[i] < x[s]) { 
        s = i;
      }
    }

    t = x[k];
    x[k] = x[s];
    x[s] = t;

  }
}

void my_fill_polygon(double *xp, double *yp, int n) {

  double y; //row of pixels
  int i; //counter for finding y intersections
  int c; //index of ycross
  int q; //counter for coloring lines

  for (y = 0.1; y < 599; y++) { //big for loop for each row of pixels
   double ycross[1000];
   c = 0;
   
   for (i = 1; i < n; i++) { 
     if ( ((yp[i] > y) && (yp[i-1] < y)) || ((yp[i] < y) && (yp[i-1] > y)) ) {
       ycross[c] = calculate_x(yp[i], xp[i], yp[i-1], xp[i-1], y);
       c = c+1;
     }
   }

   if ( ((yp[0] > y) && (yp[n-1] < y)) || ((yp[0] < y) && (yp[n-1] > y)) ) {
     ycross[c] = calculate_x(yp[0], xp[0], yp[n-1], xp[n-1], y);
     c = c+1;
   }


   sort(ycross, c);

   
   for (q = 0; q < c; q = q + 2) {
    G_line(ycross[q], y, ycross[q+1], y);
  }


}


} 

void grid() {

  int i;

  for (i = 1; i < 60; i++) {
    G_line(i*10, 0, i*10, 599);
    G_line(0, i*10, 599, i*10);
  }
  

}

int click_and_save(double *x, double *y) {

  double p[2];
  int n = 0;

  G_wait_click(p);

  while (p[0] > 50 && p[1] > 20) {

    if (fmod(p[0], 10) >= 5) {
      x[n] = round(p[0]) + (10 - fmod(p[0], 10)); 
    }
    else {
      x[n] = round(p[0]) - fmod(p[0], 10);
    }


    if (fmod(p[1], 10) >= 5) {
      y[n] = round(p[1]) + (10 - fmod(p[1], 10)); 
    }
    else {
      y[n] = round(p[1]) - fmod(p[1], 10);
    }

    G_fill_circle(x[n], y[n], 4);
    n++;
    G_wait_click(p);
  }

  return n;

}



int main() {

  double a[100], b[100], c[100], d[100];
  int m, n;
  double p[2];

  G_init_graphics(600, 600);
  grid();  

  G_rgb(0, 0, 0);
  G_fill_rectangle(0, 0, 50, 20);
  G_rgb(1, 1, 1);
  G_draw_string("Fill", 5, 0);

  G_rgb (1, 0, 0);
  m = click_and_save(a, b);
  my_fill_polygon(a, b, m);

  G_wait_click(p);
  G_close();
}