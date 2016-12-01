#include <D3d_matrix.h>
#include <FPT.h>

//This version renders XYZ objects with simple wireframes
//Some XYZ files render with viewpoint inside
//compile program with "acom -I . D3d_matrixS.c xyzreader.c"
//use arrow keys, x, y, z, i and o to manipulate object

int numobjects;
int numpoints[100];

double xp[100][10000];
double yp[100][10000];
double zp[100][10000];

double xm[100];
double ym[100];
double zm[100]; //averages of the objects' points

int numpolys[10000];
int psize[100][10000];

int con[100][10000][100];

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

void drawobject(int onum) {

  int i;
  int k;

  double thisx[10000];
  double thisy[10000];


  G_rgb(0, 0, 0);
  G_fill_rectangle(0, 0, 1000, 1000);

  G_rgb(0, 1, 0);


  for (i = 0; i < numpolys[onum]; i++) { //BIG LOOP FOR EACH POLYGON

    for (k = 0; k < psize[onum][i]; k++) {
      thisx[k] = (300 / tan(M_PI/9)) * ( xp[onum][con[onum][i][k]] / zp[onum][con[onum][i][k]] ) + 300;
      thisy[k] = (300 / tan(M_PI/9)) * ( yp[onum][con[onum][i][k]] / zp[onum][con[onum][i][k]] ) + 300; 
    }

    G_polygon(thisx, thisy, psize[onum][i]);
				  
  }
}

    

void initial_position_object(int onum) {

  double m[4][4];
  double minv[4][4];

  D3d_make_identity(m) ;  D3d_make_identity(minv) ;
  D3d_translate(m, minv, 0, 0, 10);
  D3d_mat_mult_points(xp[onum],yp[onum],zp[onum],  m, xp[onum],yp[onum],zp[onum],numpoints[onum]) ;

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

  while (key != 'q') {

    drawobject(onum);

    double m[4][4];
    double minv[4][4];

    D3d_make_identity(m) ;  D3d_make_identity(minv) ;


    key = G_wait_key();

    if (key == 65361) { //LEFT
      D3d_translate(m, minv, -.1, 0, 0);
    }
    else if (key == 65362) { //UP
      D3d_translate(m, minv, 0, .1, 0);
    }
    else if (key == 65363) { //RIGHT
      D3d_translate(m, minv, .1, 0, 0);
    }
    else if (key == 65364) { //DOWN
      D3d_translate(m, minv, 0, -.1, 0);
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

    else if (key == 'i') {//ZOOM IN
      D3d_translate(m, minv, 0, 0, -.5);
    }
    else if (key == 'o') {//ZOOM OUT
      D3d_translate(m, minv, 0, 0, .5);
     } 

    if (key > 47 && key < 58) {
      onum = key - 48;
    }

    D3d_mat_mult_points(xp[onum],yp[onum],zp[onum],  m, xp[onum],yp[onum],zp[onum],numpoints[onum]) ;

  }

}
