#include <D2d_matrix.h>
#include <FPT.h>

//compile program with "acom -I . D2d_matrixS.c xyreader.c"
//use arrow keys, space, z and x to manipulate image 

int numobjects;
int numpoints[100];

double xp[100][3000];
double yp[100][3000];

int numpolys[100];
int psize[100][2000];

int con[100][2000][20];

double red[100][200];
double green[100][200];
double blue[100][200];

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

void drawobject(int onum) {

  int i;
  int k;
  int key;

  double p[2];
  double thisx[3000];
  double thisy[3000];


  G_rgb(0, 0, 0);
  G_fill_rectangle(0, 0, 1000, 1000);


  for (i = 0; i < numpolys[onum]; i++) { //BIG LOOP FOR EACH POLYGON

    G_rgb(red[onum][i], green[onum][i], blue[onum][i]); //SET COLOR

    for (k = 0; k < psize[onum][i]; k++) {
      thisx[k] = xp[onum][con[onum][i][k]];
      thisy[k] = yp[onum][con[onum][i][k]]; 
    }

    G_fill_polygon(thisx, thisy, psize[onum][i]);

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
  drawobject(onum); //remove this if something fucks up?

  while (key != 'q') {

    drawobject(onum);

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


  D2d_mat_mult_points(xp[onum],yp[onum],  m, xp[onum],yp[onum],numpoints[onum]) ;

  }

}