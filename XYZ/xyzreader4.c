#include <D3d_matrix.h>
#include <FPT.h>

//This version has clipping
//compile program with "acom -I . D3d_matrixS.c xyzreader4.c"
//use arrow keys, x, y, z, i and o to manipulate object
//use h and j to move hither plane forward and back

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
POLYARRAY ;

int numobjects;
int numpoints[10];

double xp[10][10000];
double yp[10][10000];
double zp[10][10000];

double xm[10];
double ym[10];
double zm[10]; //averages of the objects' points

double xc[10];
double yc[10];
double zc[10]; //centers of the objects' points

int numpolys[10];
int psize[10][10000];

int con[10][10000][100];

POLYARRAY allpolys[10000];

double l[3]; //light coordinates
double inherent[3]; //inherent color - need this to work for different objects

//all these light values are from 0 to 1
double ambient;
double maxdiffuse;
double diffuse;
double specular;

double specularpower;

double hither;
double yon;

  int in_out (double p[3], int c, double plane[6][4]) {
  // return 1 if point P is "IN"
  // return 0 if point P is "OUT"

    double mx = 0;
    double my = 0;
    double mz = (yon + hither) / 2;

    double pointinsert  = plane[c][0] * p[0] + plane[c][1] * p[1] + plane[c][2] * p[2] + plane[c][3];
    double centerinsert = plane[c][0] * mx   + plane[c][1] * my   + plane[c][2] * mz   + plane[c][3];
  
    if (pointinsert > 0 && centerinsert > 0) {
      return 1;
    }
    
    else if (pointinsert < 0 && centerinsert < 0) {
      return 1;
    }
    else if (pointinsert == 0 && centerinsert == 0) {
      return 1;
    }
    else {
      // printf("test - c is %d\n", c);
      return 0;
    }

}

double clipIntersect(double p0[3], double p1[3], int c, int which, double plane[6][4]) {

//finding t value of system of equations between parametric and plane equations

  int i = 0;
  double r = 0;
  double l = 0;
  double t;

for (i = 0; i < 3; i++) { //add up non-t values
  r += - (plane[c][i] * p0[i]);
}

  r += - plane[c][3]; //d value in matrix equation

for (i = 0; i < 3; i++) { //add up t coefficients
  l += plane[c][i] * (p1[i] - p0[i]);
}

 t = r / l;

//returns calculation based on 0, 1 or 2 for x, y or z
return p0[which] + t * (p1[which] - p0[which]);

}


int clip(double x[3000], double y[3000], double z[3000], int n, int onum) {

 double newx[3000];
 double newy[3000];
 double newz[3000];

  double temp1[3]; //holds x y and z for first coordinate in in_out test
  double temp2[3]; //holds x y and z for second coordinate in in_out test

  int k; //x y and z counter
  int i; //newx newy and newz counter, return as number of sides in clipped polygon
  int c; //clipping plane counter

  
    double plane[6][4] = {
    { 0, 1, -tan(M_PI/9), 0 } , // top
    { 0, -1, -tan(M_PI/9), 0 } , // bottom
    { 1, 0, tan(M_PI/9), 0 } , // right
    { -1, 0, tan(M_PI/9), 0 } , // left
    { 0, 0, 1, -yon } , // yon
    { 0, 0, 1, -hither } //hither
    } ;    

  for (c = 0; c < 6; c++) { //each plane

    i = 0 ;
    for (k = 0; k < n; k++) { //each point in the polygon
      int kn = (k+1) % n;

      temp1[0] = x[k];
      temp1[1] = y[k];
      temp1[2] = z[k];
      
      temp2[0] = x[kn];
      temp2[1] = y[kn];
      temp2[2] = z[kn];

      if (in_out(temp1, c, plane) == 1 && in_out(temp2, c, plane) == 1) { //in to in: destination
       newx[i] = x[kn];
       newy[i] = y[kn];
       newz[i] = z[kn];

       i++;

     }

      else if (in_out(temp1, c, plane) == 1 && in_out(temp2, c, plane) == 0) { // in to out, intersection
	      newx[i] = clipIntersect(temp1, temp2, c, 0, plane);
	      newy[i] = clipIntersect(temp1, temp2, c, 1, plane);
	      newz[i] = clipIntersect(temp1, temp2, c, 2, plane);

	      

        i++;

      }

      else if (in_out(temp1, c, plane) == 0 && in_out(temp2, c, plane) == 1) { // out to in, intersection, then destination
	     newx[i] = clipIntersect(temp1, temp2, c, 0, plane);
	     newx[i + 1] = x[kn];

	     newy[i] = clipIntersect(temp1, temp2, c, 1, plane);
	     newy[i + 1] = y[kn];

	     newz[i] = clipIntersect(temp1, temp2, c, 2, plane);
	     newz[i + 1] = z[kn];

	     i = i + 2;

     }

    } // end for k

   for (k = 0; k < i; k++) {
    x[k] = newx[k];
    y[k] = newy[k];
    z[k] = newz[k];
  }

  n = i;

  } // end for c

  //printf("n is %d\n", n);

return n;
}

double light_calculate(int onum, int polyindex) {

  /////////DIFFUSE//////////////

  double Nv[3]; //normal vector
  //unit vectors:
  double Lu[3];
  double Nu[3];
  double Eu[3];
  double Ru[3]; 

  double Av[3] = 
  { xp[onum][con[onum][polyindex][1]] - xp[onum][con[onum][polyindex][0]], 
    yp[onum][con[onum][polyindex][1]] - yp[onum][con[onum][polyindex][0]], 
    zp[onum][con[onum][polyindex][1]] - zp[onum][con[onum][polyindex][0]] } ;

    double Bv[3] = 
    { xp[onum][con[onum][polyindex][2]] - xp[onum][con[onum][polyindex][0]], 
      yp[onum][con[onum][polyindex][2]] - yp[onum][con[onum][polyindex][0]], 
      zp[onum][con[onum][polyindex][2]] - zp[onum][con[onum][polyindex][0]] } ;

      double Lv[3] = 
      { l[0] - xp[onum][con[onum][polyindex][0]],
        l[1] - yp[onum][con[onum][polyindex][0]],
        l[2] - zp[onum][con[onum][polyindex][0]] } ;

        D3d_x_product(Nv, Av, Bv);

        double Lvlength = D3d_vector_distance(Lv);
        double Nvlength = D3d_vector_distance(Nv);

        Lu[0] = Lv[0]/Lvlength;
        Lu[1] = Lv[1]/Lvlength;
        Lu[2] = Lv[2]/Lvlength;

        Nu[0] = Nv[0]/Nvlength;
        Nu[1] = Nv[1]/Nvlength;
        Nu[2] = Nv[2]/Nvlength;

        double anglea = (D3d_dot_product(Nu, Lu));

        if (anglea < 0) {
          Nu[0] = -Nu[0];
          Nu[1] = -Nu[1];
          Nu[2] = -Nu[2];
          anglea = (D3d_dot_product(Nu, Lu));
        }

        diffuse = maxdiffuse * anglea;

        if (diffuse < 0) {
          diffuse = 0;
        } 


////////SPECULAR////////////////////

double Ev[3] = //eye
{ 0 - xp[onum][con[onum][polyindex][0]],
  0 - yp[onum][con[onum][polyindex][0]],
  0 - zp[onum][con[onum][polyindex][0]] } ;


double Rv[3] = //reflection from Lu
{ 2 * anglea * Nu[0] - Lu[0],
  2 * anglea * Nu[1] - Lu[1],
  2 * anglea * Nu[2] - Lu[2] } ;

  double Evlength = D3d_vector_distance(Ev);
  double Rvlength = D3d_vector_distance(Rv);

  Eu[0] = Ev[0]/Evlength;
  Eu[1] = Ev[1]/Evlength;
  Eu[2] = Ev[2]/Evlength;

  Ru[0] = Rv[0]/Rvlength;
  Ru[1] = Rv[1]/Rvlength;
  Ru[2] = Rv[2]/Rvlength;

  double angleb = D3d_dot_product(Eu, Ru);

  if (angleb < 0) {
    specular = 0;
  }
  else {
    specular = (1 - maxdiffuse - ambient) * pow(angleb, specularpower);
  }
  
////////INTENSITY////////////

  if (D3d_dot_product(Eu, Nu) < 0) {
    return ambient;
  }
  else {
    return ambient + diffuse + specular;
  }

}

void print_array(int n) {

  int i ;
  for (i = 0 ; i < n ; i++) {
    printf("%d %d %lf\n",allpolys[i].objnum, allpolys[i].polynum, allpolys[i].dist) ;
  }
  printf("\n") ;
}

int compare (const void *p, const void *q) {

  POLYARRAY *a, *b ;

  a = (POLYARRAY*)p ;
  b = (POLYARRAY*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}


int init_array() { //NOTE: ALLPOLYS OBJNUM COUNTS FROM 1, allpolysN returns

  int allpolysN = 0;
  int i;
  int j;

  for (i = 1; i <= numobjects; i++) { //COUNT TOTAL NUMBER OF POLYGONS IN "BASKET", INITIALIZE 
    for (j = 0; j < numpolys[i]; j++) {
      allpolys[allpolysN].objnum = i;
      allpolys[allpolysN].polynum = j;
      allpolys[allpolysN].dist = zp[i][con[i][j][0]];
      allpolysN = allpolysN + 1;
    }  
  }

  qsort (allpolys, allpolysN, sizeof(POLYARRAY), compare);

  return allpolysN;

}

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

double findzhi(int onum) {

  int hi = zp[onum][0];
  int i;

  for (i = 0; i < numpoints[onum]; i++) {
    if (zp[onum][i] > hi) {
      hi = zp[onum][i];
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

double findzlo(int onum) {

  int lo = zp[onum][0];
  int i;

  for (i = 0; i < numpoints[onum]; i++) {
    if (zp[onum][i] < lo) {
      lo = zp[onum][i];
    }
  }

  return lo;
}

void update_center(int onum) {

  double xhi = findxhi(onum);
  double xlo = findxlo(onum);

  double yhi = findyhi(onum);
  double ylo = findylo(onum);

  double zhi = findzhi(onum);
  double zlo = findzlo(onum);

  xc[onum] = (xhi + xlo) / 2;
  yc[onum] = (yhi + ylo) / 2;
  zc[onum] = (zhi + zlo) / 2;

}

void readobject(FILE *f, int onum) { 

  int i;
  int k;

  fscanf(f, "%d",&numpoints[onum]);
  
  for (i = 0; i < numpoints[onum]; i++) {  //LABEL POINTS
    fscanf(f, "%lf %lf %lf", &xp[onum][i], &yp[onum][i], &zp[onum][i]);
  }



  fscanf(f, "%d",&numpolys[onum]);

  for (i = 0; i < numpolys[onum]; i++) { //COLLECT POLYGON SIDES AND LINES
    fscanf(f, "%d", &psize[onum][i]);

    for (k = 0; k < psize[onum][i]; k++) {
      fscanf(f, "%d", &con[onum][i][k]);
    }
  }

}

void drawobjects() {

  int i;
  int j;
  int k;
  int l;
  int onum;

  double intensity;
  double bigK = ambient + maxdiffuse;
  double intensityratio;

  double thisx[10000];
  double thisy[10000];

  int allpolysN = init_array();

  G_rgb(0, 0, 0);
  G_fill_rectangle(0, 0, 1000, 1000);

  double cxp[3000], cyp[3000], czp[3000];

  int psides;

  int q;

  for (j = allpolysN; j >= 0; j--) { //BIG LOOP FOR EACH POLYGON

    onum = allpolys[j].objnum;
    i = allpolys[j].polynum;

    for (l = 0; l < psize[onum][i]; l++) { //copy array for polygon
      cxp[l] = xp[onum][con[onum][i][l]];
      cyp[l] = yp[onum][con[onum][i][l]];
      czp[l] = zp[onum][con[onum][i][l]];
    }

    psides = clip(cxp, cyp, czp, psize[onum][i], onum);

    for (k = 0; k < psides; k++) { //project 2D onto 3D
      thisx[k] = (300 / tan(M_PI/9)) * ( cxp[k] / czp[k] ) + 300;
      thisy[k] = (300 / tan(M_PI/9)) * ( cyp[k] / czp[k] ) + 300; 
    }

    intensity = light_calculate(onum, i);
    intensityratio = intensity / bigK;

    G_rgb(intensityratio * inherent[0], intensityratio * inherent[1], intensityratio * inherent[2]);
    G_fill_polygon(thisx, thisy, psides);

    //printf("polygon number %d has %d sides\n", j, psides);
    //G_rgb(0,0,0);
    //G_polygon(thisx, thisy, psides);

  }
}


void update_average(int onum) {

  int i;
  double tempx = 0;
  double tempy = 0; 
  double tempz = 0;

  for (i = 0; i < numpoints[onum]; i++) {
    tempx = tempx + xp[onum][i];
    tempy = tempy + yp[onum][i];
    tempz = tempz + zp[onum][i];
  }

  xm[onum] =  tempx / numpoints[onum];
  ym[onum] =  tempy / numpoints[onum];
  zm[onum] =  tempz / numpoints[onum];

}

void initial_position_object(int onum) {

  double m[4][4];
  double minv[4][4];

  update_center(onum);

  D3d_make_identity(m) ;  D3d_make_identity(minv) ;
  D3d_translate(m, minv, -xc[onum], -yc[onum], -zc[onum]);
  D3d_translate(m, minv, 0, 0, 20);
  D3d_mat_mult_points(xp[onum],yp[onum],zp[onum],  m, xp[onum],yp[onum],zp[onum],numpoints[onum]) ;

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
    initial_position_object(onum);
    
  }
  
  onum = 1;

  ///Quick Start:///
  
  ambient = .2; 
  maxdiffuse = .5; 
  specularpower = 100; 
  l[0] = 100 ; 
  l[1] = 100 ; 
  l[2] = -10 ; 

  inherent[0] = 0;
  inherent[1] = 1;
  inherent[2] = 0; 

  hither = .1;
  yon = 50;

  ///User Inputs:///

  // printf("Please enter an ambient value, a maximum diffuse value, and a specular power value\n");
  // scanf("%lf, %lf, %lf", &ambient, &maxdiffuse, &specularpower);

  // printf("Please enter the x, y and z coordinates of the light\n");
  // scanf("%lf, %lf, %lf", &l[0], &l[1], &l[2]);

  // printf("Please enter the rgb information for inherent color\n");
  // scanf("%lf, %lf, %lf", &inherent[0], &inherent[1], &inherent[2]);

  // printf("Please enter the initial hither z-value\n");
  // scanf("%lf," &hither);

  //printf("Please enter the yon z-value\n");
  //scanf("%lf," &yon);

  while (key != 'q') {

    drawobjects();

    double m[4][4];
    double minv[4][4];

    D3d_make_identity(m) ;  D3d_make_identity(minv) ;


    key = G_wait_key();

    if (key == 65361) { //LEFT
      D3d_translate(m, minv, -.3, 0, 0);
    }
    else if (key == 65362) { //UP
      D3d_translate(m, minv, 0, .3, 0);
    }
    else if (key == 65363) { //RIGHT
      D3d_translate(m, minv, .3, 0, 0);
    }
    else if (key == 65364) { //DOWN
      D3d_translate(m, minv, 0, -.3, 0);
    }

    else if (key == 'x') {//ROTATE X AXIS
      update_average(onum);      
      D3d_translate(m, minv, -xm[onum], -ym[onum], -zm[onum]);
      D3d_rotate_x(m, minv, 1*(M_PI/180));
      D3d_translate(m, minv, xm[onum], ym[onum], zm[onum]);
    }

    else if (key == 'y') {//ROTATE Y AXIS 
      update_average(onum);
      D3d_translate(m, minv, -xm[onum], -ym[onum], -zm[onum]);
      D3d_rotate_y(m, minv, 1*(M_PI/180));
      D3d_translate(m, minv, xm[onum], ym[onum], zm[onum]);
    }

    else if (key == 'z') {//ROTATE Z AXIS
      update_average(onum);      
      D3d_translate(m, minv, -xm[onum], -ym[onum], -zm[onum]);
      D3d_rotate_z(m, minv, 1*(M_PI/180));
      D3d_translate(m, minv, xm[onum], ym[onum], zm[onum]);
    }

    else if (key == 'h') {//HITHER FORWARD
      hither = hither + .7;
    }

    else if (key == 'j') {//HITHER BACKWARD
      hither = hither - .7;
    }

    else if (key == 'i') {//ZOOM IN
      D3d_translate(m, minv, 0, 0, -.7);
    }
    else if (key == 'o') {//ZOOM OUT
      D3d_translate(m, minv, 0, 0, .7);
    } 

    if (key > 47 && key < 58) {
      onum = key - 48;
    }

    D3d_mat_mult_points(xp[onum],yp[onum],zp[onum],  m, xp[onum],yp[onum],zp[onum],numpoints[onum]) ;

  }

}
