#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;

#include "GL/gl.h"
#include "GL/glu.h"
#include "SDL/SDL.h"
#include "gl2ps.h"

//Screen attributes
const int SCREEN_WIDTH = 480;
const int SCREEN_HEIGHT = 480;
const int SCREEN_BPP = 32;

//The frame rate
const int FRAMES_PER_SECOND = 20;

//Event handler
SDL_Event event;

double user_theta  = 0.1;
double user_phi = 0.;
double user_dist = 15.;
double usr_x = 0., usr_y = 0.01*user_dist, usr_z = 0.99*user_dist;
double usr_d = user_dist;
bool h_res = false; //high resolution

double ma_x, ma_y, ma_z; //section limits

int num = 0;

float col_list[5][3] = {{1.,0.3,0.3},{0.3,1.,0.3},{0.3,0.3,1.},{0.8,0.3,0.8},{0.8,0.8,0.3}};


class Particle_g {
private:
  //id
  int id;
  //radius
  double r;
  //position
  double x, y, z;
  //latitudes and longitudes
  int lats, longs;
  //vertices
  vector<double> vert_x;
  vector<double> vert_y;
  vector<vector<double> > vert_z;
  //alpha
  float A;

public:
  Particle_g();
  Particle_g(int, int);
  ~Particle_g() {}

  void set_id(int ID) { id = ID; }
  int get_id() { return id; }
  void set_r(double rr) { r = rr; }
  double get_r() { return r; }
  void set_pos(const double p[3]) { x = p[0]; y = p[1]; z = p[2]; }
  void get_pos(double p[3]) { p[0] = x; p[1] = y; p[2] = z; }
  void set_lalo(const int l[2]) { lats = l[0]; longs = l[1]; }
  void mul_lalo() { lats *= 2; longs *= 3; }
  void div_lalo() { lats /= 2; longs /= 3; }
  void vert_calc();
  void set_A(float a) { A = a; }
  float get_A() { return A; }
  void draw();
};

//The timer
class Timer {
    private:
    //The clock time when the timer started
    int startTicks;
    
    //The ticks stored when the timer was paused
    int pausedTicks;
    
    //The timer status
    bool paused;
    bool started;
    
    public:
    //Initializes variables
    Timer();
    
    //The various clock actions
    void start();
    void stop();
    void pause();
    void unpause();
    
    //Gets the timer's time
    int get_ticks();
    
    //Checks the status of the timer
    bool is_started();
    bool is_paused();    
};

///********************************************************************

Particle_g::Particle_g() : id(0), r(0.5), x(0.), y(0.), z(0.), 
		       lats(10), longs(12) {
  A = 1.;
  vert_calc();
}

Particle_g::Particle_g(int lts, int lngs) : id(0), r(0.5), x(0.), y(0.), z(0.), 
			  lats(lts), longs(lngs) {
  A = 1.;
  vert_calc();
}

void Particle_g::vert_calc() {
  vert_x.resize(longs+1);
  vert_y.resize(longs+1);
  vert_z.resize(lats+2);
  for(int i = 0; i < vert_z.size(); ++i)
    vert_z[i].resize(2);

  double rsq = r*r;
  double p4 = 4. * M_PI;

  double Pet[21];
  for(int i = 0; i <= lats+1; ++i) {
    double lat = M_PI * (-0.5 + static_cast<double>(i - 1) / lats);
    vert_z[i][0] = r * sin(lat);
    vert_z[i][1] = r * cos(lat);
  }
  for(int j = 0; j <= longs; ++j) {
    double lng = 2 * M_PI * static_cast<double>(j - 1) / longs;
    vert_x[j] = cos(lng);
    vert_y[j] = sin(lng);
  }
}

void Particle_g::draw() {
  double rsq = r * r;
  glPushMatrix();
  glTranslated(x, y, z);
  for(int i = 0; i <= lats; ++i) {
    glBegin(GL_QUAD_STRIP);
    for(int j = 0; j <= longs; ++j) {
      double x0 = vert_x[j] * vert_z[i][1];
      double x1 = vert_x[j] * vert_z[i+1][1];
      double y0 = vert_y[j] * vert_z[i][1];
      double y1 = vert_y[j] * vert_z[i+1][1];
      
      glColor4f(col_list[id][0], col_list[id][1], col_list[id][2], A);
      glNormal3f(x0 / rsq, y0 / rsq, vert_z[i][0] / rsq);
      glVertex3f(x0, y0, vert_z[i][0]);
      glNormal3f(x1 / rsq, y1 / rsq, vert_z[i+1][0] / rsq);
      glVertex3f(x1, y1, vert_z[i+1][0]);
    }
    glEnd();
  }
  glPopMatrix();
}


Timer::Timer() {
  //Initialize the variables
  startTicks = 0;
  pausedTicks = 0;
  paused = false;
  started = false;    
}

void Timer::start() {
  //Start the timer
  started = true;
  
  //Unpause the timer
  paused = false;
  
  //Get the current clock time
  startTicks = SDL_GetTicks();    
}

void Timer::stop() {
  //Stop the timer
  started = false;
  
  //Unpause the timer
  paused = false;    
}

void Timer::pause() {
  //If the timer is running and isn't already paused
  if( ( started == true ) && ( paused == false ) ){
    //Pause the timer
    paused = true;
    
    //Calculate the paused ticks
    pausedTicks = SDL_GetTicks() - startTicks;
  }
}

void Timer::unpause() {
  //If the timer is paused
  if( paused == true ){
    //Unpause the timer
    paused = false;
    
    //Reset the starting ticks
    startTicks = SDL_GetTicks() - pausedTicks;
    
    //Reset the paused ticks
    pausedTicks = 0;
  }
}

int Timer::get_ticks() {
  //If the timer is running
  if( started == true ){
    //If the timer is paused
    if( paused == true ){
      //Return the number of ticks when the timer was paused
      return pausedTicks;
    }else{
      //Return the current time minus the start time
      return SDL_GetTicks() - startTicks;
    }    
  }
  
  //If the timer isn't running
  return 0;    
}

bool Timer::is_started() {
  return started;    
}

bool Timer::is_paused() {
  return paused;    
}





void computeLocation() {
  glMatrixMode(GL_PROJECTION); // Set projection parameters.
  //glPushMatrix();
  glLoadIdentity();
  gluPerspective(60., 1., 0.1, 100.);
  glTranslatef(0., 0., - usr_d);
  gluLookAt(usr_x, usr_y, usr_z,  0, 0, 0,  0, 0, 1);
  //glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}

// Initializes information for drawing within OpenGL.
bool init_GL() {
  GLfloat light0_direction[] = { 20.0, 20.0, 20.0, 0.0 };
  GLfloat light0_intensity[] = { 0.05, 0.05, 0.05, 1.0 };
  //GLfloat light1_direction[] = { -20.0, -20.0, 20.0, 0.0 };
  //GLfloat light1_intensity[] = { 0.005, 0.005, 0.005, 1.0 };
  GLfloat ambient_intensity[] = { 0.6, 0.6, 0.6, 1.0 };
  GLfloat mat_specular[] = { 0.1, 0.1, 0.1, 1.0 };
  GLfloat mat_shine[] = { 0.1 };
  
  glClearColor(1.0, 1.0, 1.0, 0.0);   // Set window color to white.
  computeLocation();
  
  glEnable(GL_DEPTH_TEST);            // Draw only closest surfaces
  
  glEnable(GL_LIGHTING);              // Set up ambient light.
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_intensity);
  
  glEnable(GL_LIGHT0);                // Set up sunlight.
  glLightfv(GL_LIGHT0, GL_POSITION, light0_direction);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_intensity);
  //glEnable(GL_LIGHT1);                // Set up sunlight.
  //glLightfv(GL_LIGHT1, GL_POSITION, light1_direction);
  //glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_intensity);
  
  
  glEnable(GL_COLOR_MATERIAL);        // Configure glColor().
  //glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shine);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //If there was any errors
  if( glGetError() != GL_NO_ERROR ){
    return false;    
  }
  //If everything initialized
  return true;
}

//gets input
bool input(int& Ncomp, vector<int>& Npart, vector<double>& rad){
  Npart.clear();
  rad.clear();
  string line;
  
  ifstream In("input.txt");
  if(! In){
    cout << "Couldn't open " << "input.txt" << endl;
    In.close();
    return false;
  }
  In >> line >> line >> line >> line;
  In >> line >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line;
  In >> line >> line >> line >> Ncomp;
  In >> line >> line >> line;
  for(int i = 0; i < Ncomp; ++i){
    int a;
    In >> a;
    Npart.push_back(a);
  }
  In >> line;
  for(int i = 0; i < Ncomp; ++i)
    In >> line;
  In >> line;
  for(int i = 0; i < Ncomp; ++i){
    double a;
    In >> a;
    rad.push_back(a);
  }
  In.close();
  return true;
}

//charges configuration
bool chrg_conf(vector<vector<Particle_g> >& part, 
	       const vector<int>& Npart, int Ncomp, 
	       double L[3]){
  char name[30];
  sprintf(name, "conf.dat");
  ifstream InPos(name);
  if(! InPos){
    InPos.close();
    return false;
  }
  double pi[3], dummy;
  InPos >> dummy;
  for(int i = 0; i < 10; ++i)
    InPos >> dummy;
  for(int n = 0; n < Ncomp; ++n)
    InPos >> dummy;
  InPos >> L[0] >> L[1] >> L[2];
  for(int n = 0; n < Ncomp; ++n){
    int NN = Npart[n];
    for(int i = 0; i < NN; ++i){
      InPos >> pi[0] >> pi[1] >> pi[2]
	    >> dummy >> dummy >> dummy;
      for(int a = 0; a < 3; ++a)
	pi[a] -= L[a] * floor(pi[a] / L[a] + 0.5);
      part[n][i].set_pos(pi);
      part[n][i].vert_calc();
    }
  }
  InPos.close();
  return true;
}

void draw_box(const double L[3]){
  double L2[3] = {L[0]/2, L[1]/2, L[2]/2};
  glPushMatrix();
  glBegin(GL_LINE_LOOP);
  glVertex3f(-L2[0], -L2[1], -L2[2]);
  glVertex3f(L2[0], -L2[1], -L2[2]);
  glVertex3f(L2[0], L2[1], -L2[2]);
  glVertex3f(-L2[0], L2[1], -L2[2]);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f(-L2[0], -L2[1], L2[2]);
  glVertex3f(L2[0], -L2[1], L2[2]);
  glVertex3f(L2[0], L2[1], L2[2]);
  glVertex3f(-L2[0], L2[1], L2[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(-L2[0], -L2[1], -L2[2]);
    glVertex3f(-L2[0], -L2[1], L2[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(L2[0], -L2[1], -L2[2]);
    glVertex3f(L2[0], -L2[1], L2[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(L2[0], L2[1], -L2[2]);
    glVertex3f(L2[0], L2[1], L2[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(-L2[0], L2[1], -L2[2]);
    glVertex3f(-L2[0], L2[1], L2[2]);
  glEnd();
  glPopMatrix();
}

// Draws the current image.
void show(vector<vector<Particle_g> >& part, 
	  const vector<int>& Npart, int Ncomp, const double L[3], 
	  bool toggle[]){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear window.
  glColor3f(1., 1., 1.);
  glShadeModel(/*GL_FLAT*/GL_SMOOTH);
  for(int n = 0; n < Ncomp; ++n){
    if(toggle[n]){
      int NN = Npart[n];
      double p[3];
      for(int i = 0; i < NN; ++i){
	part[n][i].get_pos(p);
	if(p[0] <= ma_x && p[1] <= ma_y && p[2] <= ma_z){
	  part[n][i].draw();
	}
      }
    }
  }
  glColor3f(0.6, 0.3, 0.3);
  draw_box(L);
}

void calc_usr(){
  double cost = cos(user_theta);
  double sint = sin(user_theta);
  double cosp = cos(user_phi);
  double sinp = sin(user_phi);
  usr_x = user_dist * sint * cosp;
  usr_y = user_dist * sint * sinp;
  usr_z = user_dist * cost;
  usr_d = sqrt(usr_x * usr_x + usr_y * usr_y + usr_z * usr_z); // distance to origin
  //Adjust the view
  computeLocation();
}

//writes an output ps
void writefile(int format, int sort, int options, int nbcol, char *filename, 
	       const char *extension, 
	       vector<vector<Particle_g> >& part, 
	       const vector<int>& Npart, int Ncomp, double L[3], 
	       bool toggle[]){
  FILE *fp;
  char file[256];
  int state = GL2PS_OVERFLOW, buffsize = 0;
  GLint viewport[4];

  strcpy(file, filename);
  strcat(file, ".");
  strcat(file, extension);

  viewport[0] = 0;
  viewport[1] = 0;
  viewport[2] = SCREEN_WIDTH;
  viewport[3] = SCREEN_HEIGHT;

  fp = fopen(file, "wb");

  if(!fp){
    printf("Unable to open file %s for writing\n", file);
    exit(1);
  }

  printf("Saving image to file %s... ", file);
  fflush(stdout);

  while(state == GL2PS_OVERFLOW){
    buffsize += 1024*1024;
    gl2psBeginPage(file, "gl2psTest", viewport, format, sort, options,
                   GL_RGBA, 0, NULL, nbcol, nbcol, nbcol,
                   buffsize, fp, file);
    show(part, Npart, Ncomp, L, toggle);
    state = gl2psEndPage();
  }

  fclose(fp);

  printf("Done!\n");
  fflush(stdout);
}

void handle_input(vector<vector<Particle_g> >& part, 
		  const vector<int>& Npart, int Ncomp, double L[3], 
		  bool toggle[]){
  char ext[32];
  char name[32];
  float a, val;
  int NN;
  //If a key was pressed
  if( event.type == SDL_KEYDOWN ){
    //Adjust the velocity
    switch( event.key.keysym.sym ){
    case SDLK_UP: 
      user_theta += 0.1;
      if(user_theta > M_PI)
	user_theta = M_PI;
      calc_usr();
      break;
    case SDLK_DOWN:
      user_theta -= 0.1;
      if(user_theta < 0.1)
	user_theta = 0.1;
      calc_usr();
      break;
    case SDLK_LEFT: 
      user_phi  += 0.1;
      calc_usr();
      break;
    case SDLK_RIGHT:
      user_phi  -= 0.1;
      calc_usr();
      break;
    case SDLK_PAGEUP:
      user_dist += 0.5;
      calc_usr();
      break;
    case SDLK_PAGEDOWN: 
      user_dist -= 0.5;
      calc_usr();
      break;
    case SDLK_SPACE:
      if(! chrg_conf(part, Npart, Ncomp, L))
    	cout << "problem opening ptPosF.dat" << endl;
      break;
    case SDLK_h:
      if(! h_res){
	for(int n = 0; n < Ncomp; ++n){
	  int NN = Npart[n];
	  for(int i = 0; i < NN; ++i){
	    part[n][i].mul_lalo();
	    part[n][i].vert_calc();
	  }
	}
	h_res = true;
      }else{
	for(int n = 0; n < Ncomp; ++n){
	  int NN = Npart[n];
	  for(int i = 0; i < NN; ++i){
	    part[n][i].div_lalo();
	    part[n][i].vert_calc();
	  }
	}
	h_res = false;
      }
      break;
    case SDLK_KP0:
      if(Ncomp > 0){
	if(toggle[0])
	  toggle[0] = false;
	else
	  toggle[0] = true;
      }
      break;
    case SDLK_KP1:
      if(Ncomp > 1){
	if(toggle[1])
	  toggle[1] = false;
	else
	  toggle[1] = true;
      }
      break;
    case SDLK_KP2:
      if(Ncomp > 2){
	if(toggle[2])
	  toggle[2] = false;
	else
	  toggle[2] = true;
      }
      break;
    case SDLK_KP3:
      if(Ncomp > 3){
	if(toggle[3])
	  toggle[3] = false;
	else
	  toggle[3] = true;
      }
      break;

    case SDLK_KP4:
      if(Ncomp > 4){
	if(toggle[4])
	  toggle[4] = false;
	else
	  toggle[4] = true;
      }
      break;
    case SDLK_KP5:
      if(Ncomp > 5){
	if(toggle[5])
	  toggle[5] = false;
	else
	  toggle[5] = true;
      }
      break;
    case SDLK_KP6:
      if(Ncomp > 6){
	if(toggle[6])
	  toggle[6] = false;
	else
	  toggle[6] = true;
      }
      break;

    case SDLK_p:
      sprintf(name, "configuration%d", num);
      strcpy(ext, gl2psGetFileExtension(GL2PS_EPS));
      writefile(GL2PS_EPS, GL2PS_SIMPLE_SORT, GL2PS_NONE, 0, name, ext, 
		part, Npart, Ncomp, L, toggle);
      num++;
      break;
    case SDLK_i:
      if(ma_x > -L[0]/2)
	ma_x -= L[0]/20;
      break;
    case SDLK_j:
      if(ma_x < L[0]/2)
	ma_x += L[0]/20;
      break; 
    case SDLK_m:
      if(ma_y > -L[1]/2)
	ma_y -= L[1]/20;
      break;
    case SDLK_n:
      if(ma_y < L[1]/2)
	ma_y += L[1]/20;
      break;
    case SDLK_k:
      if(ma_z > -L[2]/2)
	ma_z -= L[2]/20;
      break;
    case SDLK_l:
      if(ma_z < L[2]/2)
	ma_z += L[2]/20;
      break;
    case SDLK_F1:
      if(Ncomp > 0){
	a = part[0][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[0];
	for(int i = 0; i < NN ; ++i){
	  part[0][i].set_A(val);
	}
      }
      break;
    case SDLK_F2:
      if(Ncomp > 1){
	a = part[1][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[1];
	for(int i = 0; i < NN ; ++i){
	  part[1][i].set_A(val);
	}
      }
      break;
    case SDLK_F3:
      if(Ncomp > 2){
	a = part[2][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[2];
	for(int i = 0; i < NN ; ++i){
	  part[2][i].set_A(val);
	}
      }
      break;
    case SDLK_F4:
      if(Ncomp > 3){
	a = part[3][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[3];
	for(int i = 0; i < NN ; ++i){
	  part[3][i].set_A(val);
	}
      }
      break;

    case SDLK_F5:
      if(Ncomp > 4){
	a = part[4][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[4];
	for(int i = 0; i < NN ; ++i){
	  part[4][i].set_A(val);
	}
      }
      break;
    case SDLK_F6:
      if(Ncomp > 5){
	a = part[5][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[5];
	for(int i = 0; i < NN ; ++i){
	  part[5][i].set_A(val);
	}
      }
      break;
    case SDLK_F7:
      if(Ncomp > 6){
	a = part[6][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[6];
	for(int i = 0; i < NN ; ++i){
	  part[6][i].set_A(val);
	}
      }
      break;
      
    }

    //Show the scene
    show(part, Npart, Ncomp, L, toggle);
  }
}

bool init(){
  //Initialize SDL
  if( SDL_Init( SDL_INIT_EVERYTHING ) < 0 ){
    return false;    
  }
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);//disable(0)-enable(1) double buffering
  //Create Window
  if( SDL_SetVideoMode( SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, 
			SDL_OPENGL ) == NULL ){
    return false;
  }
  //Initialize OpenGL
  if( init_GL() == false ){
    return false;    
  }
  //Set caption
  SDL_WM_SetCaption( "Visualization", NULL );
  return true;
}

void clean_up(){
  //Quit SDL
  SDL_Quit();
}


int main( int argc, char *argv[] ){
  //Quit flag
  bool quit = false;
  //Initialize
  if( init() == false ){
    return 1;    
  }

  int Ncomp;
  vector<int> Npart;
  vector<double> rad;
  input(Ncomp, Npart, rad);
  // //*****
  // cout << Ncomp << endl;
  // for(int i = 0; i < Ncomp; ++i)
  //   cout << Npart[i] << " ";
  // cout << endl;
  // for(int i = 0; i < Ncomp; ++i)
  //   cout << rad[i] << " ";
  // cout << endl;
  // //*****
  
  //The frame rate regulator
  Timer fps;

  //the particles
  vector<vector<Particle_g> > part(Ncomp);
  for(int n = 0; n < Ncomp; ++n){
    Particle_g pt;
    for(int i = 0; i < Npart[n]; ++i)
      part[n].push_back(pt);
  }
  double r;
  int l[2];
  bool toggle[Ncomp]; //for toggling components

  for(int n = 0; n < Ncomp; ++n){
    if(n==0)
      toggle[n] = true;
    else
      toggle[n] = false;
    int NN = Npart[n];
    r = rad[n];
    l[0] = 10; l[1] = 10;
    for(int i = 0; i < NN; ++i){
      part[n][i].set_id(n);
      part[n][i].set_r(r);
      part[n][i].set_lalo(l);
    }
  }

  //charge initiel configuration
  double L[3] = {0.1, 0.1, 0.1};
  if(! chrg_conf(part, Npart, Ncomp, L))
    cout << "problem opening conf.dat" << endl;
  //section limits
  ma_x = L[0]/2; ma_y = L[1]/2; ma_z = L[2]/2;   
  //Show the scene
  show(part, Npart, Ncomp, L, toggle);
  //Update screen
  SDL_GL_SwapBuffers();

  //Wait for user exit
  while( quit == false ){
    //Start the frame timer
    fps.start();
    //While there are events to handle
    while( SDL_PollEvent( &event ) ){
      //Handle key presses
      handle_input(part, Npart, Ncomp, L, toggle);

      if( event.type == SDL_QUIT ){
	quit = true;
      }
    }
    
    // //to animate
    // if(! chrg_conf(part, Npart, Ncomp, L))
    //   cout << "problem opening conf.dat" << endl;
    // show(part, Npart, Ncomp, L, toggle);
    // //

    //Update screen
    SDL_GL_SwapBuffers();

    //Cap the frame rate
    if( fps.get_ticks() < 1000 / FRAMES_PER_SECOND ){
      SDL_Delay( ( 1000 / FRAMES_PER_SECOND ) - fps.get_ticks() );
    }
  }

  //Clean up
  clean_up();

  return 0;
}
