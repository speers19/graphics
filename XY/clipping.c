#include <D2d_matrix.h>
#include <FPT.h>

//compile program with "acom -I . D2d_matrixS.c clipping.c"
//click points of convex clipping window, hit CLIP, use arrow keys, space, z and x to manipulate image 

int numobjects;
int numpoints[100];

double xp[100][3000];
double yp[100][3000];

double cxp[100][3000];
double cyp[100][3000];
int numclippoints[100];

int numpolys[100];
int psize[100][2000];

int con[100][2000][20];

double red[100][200];
double green[100][200];
double blue[100][200];

int firsttime = 1;

double findxhi(int onum) {

  int hi = xp[onum][0];
  int i;

  for (i = 0; i < numpoints[onum]; i++) {
    if (xp[onum][i] > hi) {
      hi = xp[onum][i];
    }
  }

  return hi;
}

double findyhi(int onum) {

  int hi = yp[onum][0];
  int i;

  for (i = 0; i < numpoints[onum]; i++) {
    if (yp[onum][i] > hi) {
      hi = yp[onum][i];
    }
  }

  return hi;
}

double findxlo(int onum) {

  int lo = xp[onum][0];
  int i;

  for (i = 0; i < numpoints[onum]; i++) {
    if (xp[onum][i] < lo) {
      lo = xp[onum][i];
    }
  }

  return lo;
}

double findylo(int onum) {

  int lo = yp[onum][0];
  int i;

  for (i = 0; i < numpoints[onum]; i++) {
    if (yp[onum][i] < lo) {
      lo = yp[onum][i];
    }
  }

  return lo;
}

void readobject(FILE *f, int onum) { 

  int i;
  int k;

  fscanf(f, "%d",&numpoints[onum]);
  
  ///// LABEL POINTS
  for (i = 0; i < numpoints[onum]; i++) {
    fscanf(f, "%lf %lf", &xp[onum][i], &yp[onum][i]);
  }
  /////


  fscanf(f, "%d",&numpolys[onum]);

  ///// COLLECT POLYGON SIDES AND LINES
  for (i = 0; i < numpolys[onum]; i++) {
    fscanf(f, "%d",&psize[onum][i]);

    for (k = 0; k < psize[onum][i]; k++) {
      fscanf(f, "%d", &con[onum][i][k]);
    }
  }
  /////


  ///// COLLECT COLORS
  for (i = 0; i < numpolys[onum]; i++) {
    fscanf(f, "%lf %lf %lf", &red[onum][i], &green[onum][i], &blue[onum][i]);
  }
  //////
}

int in_out (int onum, double p[2], int c) {
  // return 1 if point P is "IN"
  // return 0 if point P is "OUT"


  int i;
  double m;
  double b;
  int n = numclippoints[onum];
  int cn = (c+1) % n;

  double mx = 0;
  double my = 0; 

  for (i = 0; i < n; i++) { //CENTER OF MASS
    mx = mx + cxp[onum][i];
    my = my + cyp[onum][i];
  }

  mx = mx/n;
  my = my/n;

    m = ( cyp[onum][c] - cyp[onum][cn] ) / ( cxp[onum][c] - cxp[onum][cn]  + .001 ) ;
    b = cyp[onum][c] - ( m * cxp[onum][c] ) ;

    if ( ( m * p[0] + b ) < p[1] && ( m * mx + b ) < my) {
      return 1;
    }
    else if ( ( m * p[0] + b ) > p[1] && ( m * mx + b ) > my) {
      return 1;
    }
    else {
      return 0; 
    }

  }


double clipIntersectX(int onum, double p0[2], double p1[2], int c) {

  int n = numclippoints[onum];
  int cn = ( c + 1 ) % n;

  double lm = ( p0[1] - p1[1] ) / ( p0[0] - p1[0] + .01) ;
  double cm = ( cyp[onum][c] - cyp[onum][cn] ) / ( cxp[onum][c] - cxp[onum][cn] + .01) ;

  double lb = p0[1] - ( lm * p0[0] ) ;
  double cb = cyp[onum][c] - ( cm * cxp[onum][c] ) ; 

  return ( ( cb - lb ) / ( lm - cm ) ) ;
}

double clipIntersectY(int onum, double p0[2], double p1[2], int c) { 

  int n = numclippoints[onum];
  int cn = ( c + 1 ) % n;

  double lm = ( p0[1] - p1[1] ) / ( p0[0] - p1[0] + .01 ) ;
  double cm = ( cyp[onum][c] - cyp[onum][cn] ) / ( cxp[onum][c] - cxp[onum][cn] + .01) ;

  double lb = p0[1] - ( lm * p0[0] ) ;
  double cb = cyp[onum][c] - ( cm * cxp[onum][c] ) ;

  return ( cm * clipIntersectX(onum, p0, p1, c) + cb  ) ;
}

int clip(double x[3000], double y[3000], int n, int onum) {

  double newx[3000];
  double newy[3000];

  double temp1[2]; //holds x and y for first coordinate in in_out test
  double temp2[2]; //holds x and y for second coordinate in in_out test

  int k; //x and y counter
  int i; //newx and newy counter, return as number of sides in clipped polygon
  int c; //clipping line counter

  for (c = 0; c < numclippoints[onum]; c++) {

    for (k = 0; k < n; k++) {

      int kn = (k+1) % n;

      temp1[0] = x[k];
      temp1[1] = y[k];
      
      temp2[0] = x[kn];
      temp2[1] = y[kn];
      
      if (in_out(onum, temp1, c) == 1 && in_out(onum, temp2, c) == 1) { //in to in: destination
        newx[i] = x[kn];

        newy[i] = y[kn];

        i++;
      }
      
      else if (in_out(onum, temp1, c) == 1 && in_out(onum, temp2, c) == 0) { // in to out, intersection
        newx[i] = clipIntersectX(onum, temp1, temp2, c);

        newy[i] = clipIntersectY(onum, temp1, temp2, c);

        i++;
      }

      else if (in_out(onum, temp1, c) == 0 && in_out(onum, temp2, c) == 1) {//intersection, then destination
        newx[i] = clipIntersectX(onum, temp1, temp2, c);
        newx[i + 1] = x[kn];

        newy[i] = clipIntersectY(onum, temp1, temp2, c);
        newy[i + 1] = y[kn];

        i = i + 2;
      }

    }

    for (k = 0; k < i; k++) { //turn thisx and thisy into newx and newy
      x[k] = newx[k];
      y[k] = newy[k];
    }
    
    n = i;
    i = 0;

  }

  return n;
}

void drawobject(int onum) {

  int i;
  int k;
  int key;

  double p[2];
  double thisx[3000];
  double thisy[3000];

  int psides;


  G_rgb(0, 0, 0);
  G_fill_rectangle(0, 0, 1000, 1000);


  for (i = 0; i < numpolys[onum]; i++) { //BIG LOOP FOR EACH POLYGON

    G_rgb(red[onum][i], green[onum][i], blue[onum][i]); //SET COLOR
    
    for (k = 0; k < psize[onum][i]; k++) {
      thisx[k] = xp[onum][con[onum][i][k]];
      thisy[k] = yp[onum][con[onum][i][k]]; 
    }

    psides = psize[onum][i];

    if (firsttime > 1) { //CLIP POINTS IF CLIPPING WINDOW IS SET 
      psides = clip(thisx, thisy, psize[onum][i], onum);
    }

    G_fill_polygon(thisx, thisy, psides);

  }
} 

void position_object(int onum) {

  double xhi, xlo, yhi, ylo;  
  double xm;
  double ym;

  double width;
  double height;
  double xs, ys;

  double m[3][3];
  double minv[3][3];

  D2d_make_identity(m) ;  D2d_make_identity(minv) ;

  xhi = findxhi(onum);
  xlo = findxlo(onum);

  yhi = findyhi(onum);
  ylo = findylo(onum);

  xm = (xhi + xlo) / 2;
  ym = (yhi + ylo) / 2;

  width = xhi - xlo;
  height = yhi - ylo;

  xs = 600/width;
  ys = 600/height;

  D2d_translate(m, minv, -xm, -ym);
  D2d_scale(m, minv, xs, ys);
  D2d_translate(m, minv, 300, 300); 

  D2d_mat_mult_points(xp[onum],yp[onum],  m, xp[onum],yp[onum],numpoints[onum]) ;
}

int click_clipping_window(int onum) {

  double p[2];
  int n = 0;

  G_wait_click(p);

  while (p[0] > 50 && p[1] > 20) {

    cxp[onum][n] = p[0];
    cyp[onum][n] = p[1];
    
    G_point(cxp[onum][n], cyp[onum][n]);

    if (n > 0) {
      G_line(cxp[onum][n], cyp[onum][n], cxp[onum][n-1], cyp[onum][n-1]);
    }

    n++;

    G_wait_click(p);

  }

  G_line(cxp[onum][n-1], cyp[onum][n-1], cxp[onum][0], cyp[onum][0]);

  return n;
}

void draw_clipping_window(int onum) {

  G_rgb(1, 1, 1);

  int n;

  for (n = 1; n < numclippoints[onum]; n++) {
    G_line(cxp[onum][n], cyp[onum][n], cxp[onum][n-1], cyp[onum][n-1]);
  } 

  G_line(cxp[onum][n-1], cyp[onum][n-1], cxp[onum][0], cyp[onum][0]);
}


int main(int argc, char **argv) {

  int onum;
  int n;
  char c;
  int access;
  double p[2];
  FILE *f;
  int key = 48;

  numobjects = argc - 1;

  G_init_graphics(600,600);

  for (onum = 1; onum <= numobjects; onum++) {

    f = fopen(argv[onum], "r");

    if (f == NULL) {
      printf("Can't open %s\n", argv[onum]);
      exit(0); 
    }

    readobject(f, onum);
    position_object(onum);

  }

  onum = 1;

  while (key != 'q') {

    drawobject(onum);
    draw_clipping_window(onum);

    if (firsttime == 1) {
      G_rgb(0, 0, 0);
      G_fill_rectangle(0, 0, 50, 20);
      G_rgb(1, 1, 1);
      G_draw_string("CLIP", 5, 0);
      numclippoints[onum] = click_clipping_window(onum);
    }

    double m[3][3];
    double minv[3][3];

    D2d_make_identity(m) ;  D2d_make_identity(minv) ;


    key = G_wait_key();

    if (key == 65361) { //LEFT
      D2d_translate(m, minv, -10, 0);
    }
    else if (key == 65362) { //UP
      D2d_translate(m, minv, 0, 10);
    }
    else if (key == 65363) { //RIGHT
      D2d_translate(m, minv, 10, 0);
    }
    else if (key == 65364) { //DOWN
      D2d_translate(m, minv, 0, -10);
    }
    else if (key == 32) {//ROTATE with space
      D2d_translate(m, minv, -300, -300);
      D2d_rotate(m, minv, .05);
      D2d_translate(m, minv, 300, 300);
    }
    else if (key == 'z') {//ZOOM IN
      D2d_translate(m, minv, -300, -300);
      D2d_scale(m, minv, 1.1, 1.1);
      D2d_translate(m, minv, 300, 300);
    }
    else if (key == 'x') {//ZOOM OUT
      D2d_translate(m, minv, -300, -300);
      D2d_scale(m, minv, .9, .9);
      D2d_translate(m, minv, 300, 300);
    }

    if (key > 47 && key < 58) {
      onum = key - 48;
    }

    firsttime++;

    D2d_mat_mult_points(xp[onum],yp[onum],  m, xp[onum],yp[onum],numpoints[onum]) ;
  }

}
