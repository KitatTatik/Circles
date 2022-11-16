#include <SDL2/SDL.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>

#define bool char
#define true 1
#define false 0
#define maxn 100
#define WIN_CAPTION " - CURRENT BALLS"


//  STRUCTURES

typedef struct {
    double       ball_x;
    double       ball_y;
    double       speed_x;
    double       speed_y;
    double       rad;
	double       mass;
	double       cr;
	double       cg;
	double       cb;
    } shar;

shar ball[maxn];

typedef struct {
    SDL_Window     *window;
    SDL_Renderer   *rend;
    SDL_Event       event;
    double          deltat;
	double          t0;
	double          gy;
	double          k_visc1;
	double          k_visc2;
    shar            ball[maxn];
    } Application;

//   FUNCTIONS - DRAWING PREPARATION

int init_window(Application* app) // window to show
{
    if(SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        printf("Failed to initialize the SDL2 library\n");
        printf("SDL2 Error: %s\n", SDL_GetError());
        return 0;
    }
    
    app->window = SDL_CreateWindow("Shariki",
                                   SDL_WINDOWPOS_CENTERED,
                                   SDL_WINDOWPOS_CENTERED,
                                   680, 480,
                                   0);
    if(!app->window)
    {
        printf("Failed to create window\n");
        printf("SDL2 Error: %s\n", SDL_GetError());
        return 0;
    }
    SDL_Delay(50);
    return 1;
}

int get_rend(Application* app) //render to rend into the window
{
    app->rend = SDL_CreateRenderer(app->window, -1, SDL_RENDERER_ACCELERATED);  

    if(!app->rend)
    {
        printf("Failed to get the renderer to the window\n");
        printf("SDL2 Error: %s\n", SDL_GetError());
        return 0;
    }
    SDL_SetRenderDrawColor(app->rend, 97, 205, 207, 0.5);

    return 1;
}


int init_ball(Application* app,  int k) // balls parameters
{
    double rk=0.0;
    srand(time(NULL));

    app->ball[k].speed_x = 1.0*(1.+2*k);
    app->ball[k].speed_y = 0.0;
    
    rk = 8.0+3*k; // temporary - balls size
    if (k>8)rk=8.0+2*k;
    app->ball[k].rad=rk;
    app->ball[k].cr=rand()%(255);
    app->ball[k].cg=rand()%(255);
    app->ball[k].cb=rand()%(255);
    app->ball[k].mass = 4.189*rk*rk*rk*0.25; // temporary - balls density
    return(0);
}

void add_ball(Application* app,  int k, int x,int y)
{
    init_ball(app,k);
    app->ball[k].ball_x = x;
    app->ball[k].ball_y = y;
}



//      FUNCTIONS - evolution description


void border(Application* app, int k)   //  reflections from the borders
{
    double a,b,c,d;
    double r;
    double x,y,vx,vy;
    double mod,spot; // spot - ne pro speed, a pro overlap
    double a0 = 678.0; // borders' actual screen postions
    double c0 = 2.0;
    double b0 = 478.0;
    double d0 = 2.0;

    x=app->ball[k].ball_x;
    y=app->ball[k].ball_y;
    vx=app->ball[k].speed_x;
    vy=app->ball[k].speed_y;
    r=app->ball[k].rad;

    a=x+r;
    b=y+r;
    c=x-r;
    d=y-r;    
    
    if(a > a0)
    {
        mod=fabs(vx);
        spot=(a-a0);
        if (vx>0.)vx=-vx;
        x=x+vx*spot/mod;
    }
    if(c < c0)
    {
        mod=fabs(vx);
        spot=(c0-c);
        if (vx<0.)vx=-vx;
        x=x+vx*spot/mod;
    }
    if(b >b0) 
    {
        mod=fabs(vy);
        spot=(b-b0);
        if (vy>0.)vy=-vy;
        y=y+vy*spot/mod;
    }
    if(d < d0) 
    {
        mod=fabs(vy);
        spot=(d0-d);
        if (vy<0.)vy=-vy;
        y=y+vy*spot/mod;
    }
    
    app->ball[k].speed_x=vx; 
    app->ball[k].ball_x=x;
    app->ball[k].speed_y=vy;
    app->ball[k].ball_y=y;
}

void rotation (double *x,double *y, double sin, double cos,int sign) // matrix rotation need below in collision
{
    double x1,y1;

    x1=(*x*cos - sign*(*y)*sin);
    y1=(*y*cos + sign*(*x)*sin);
    *x=x1;
    *y=y1;
}

void collision(Application* app,  int i, int j) //handles balls' collision
{
    double xi, xj, yi, yj, spxi, spyi, spxj, spyj,ri,rj,mi,mj;
    double dx,dy,dist,angle,sine,cosine;
    double sptot,mod,spot;
    double dt=app->deltat;
    double xi1,xj1,yi1,yj1,dx1,dy1;
   
    xi=app->ball[i].ball_x;
    xj=app->ball[j].ball_x;
    yi=app->ball[i].ball_y;
    yj=app->ball[j].ball_y;
    
    spxi=app->ball[i].speed_x;
    spxj=app->ball[j].speed_x;
    spyi=app->ball[i].speed_y;
    spyj=app->ball[j].speed_y;
    
    ri=app->ball[i].rad;
    rj=app->ball[j].rad;
    
    mi=app->ball[i].mass;
    mj=app->ball[j].mass;

    xi1=xi+dt*spxi;
    yi1=yi+dt*spyi;
    xj1=xj+dt*spxj;
    yj1=yj+dt*spyj;

    dx=xj-xi;
    dy=yj-yi;

    dx1=xj1-xi1;
    dy1=yj1-yi1;
    dist=sqrt(dx1*dx1+dy1*dy1); 
   
    if (dist < (ri+rj))
    {
        angle=atan2(dy,dx);
        sine=sin(angle);
        cosine=cos(angle);
  
        xi=0.0;
        yi=0.0;
        rotation(&dx,&dy,sine,cosine,-1);
        xj=dx;
        yj=dy;
        rotation(&spxi,&spyi,sine,cosine,-1);
        rotation(&spxj,&spyj,sine,cosine,-1);
        sptot=spxi-spxj;
        spxi=((mi-mj)*spxi+2.0*mj*spxj)/(mi+mj);
        spxj= spxi+sptot;

        mod=fabs(spxi)+fabs(spxj);  //avoided crossing
        spot=(ri+rj)-fabs(dx);
        xi=xi+spxi*spot/mod;
        xj=xj+spxj*spot/mod;

        rotation(&xi,&yi,sine,cosine,1);
        rotation(&xj,&yj,sine,cosine,1);	    
        rotation(&spxi,&spyi,sine,cosine,1);
        rotation(&spxj,&spyj,sine,cosine,1);

        app->ball[i].ball_x=app->ball[i].ball_x+xi;
        app->ball[i].ball_y=app->ball[i].ball_y+yi;
        app->ball[j].ball_x=app->ball[i].ball_x+xj;
        app->ball[j].ball_y=app->ball[i].ball_y+yj;

        app->ball[i].speed_x = spxi;
        app->ball[j].speed_x = spxj;
        app->ball[i].speed_y = spyi;
        app->ball[j].speed_y = spyj;
    }
}

void shift(Application* app, int i)  // provides dynamics
{
    
    double dt=app->deltat;
    double r=app->ball[i].rad;
    double x=app->ball[i].ball_x;
    double y=app->ball[i].ball_y;
    double vx=app->ball[i].speed_x;
    double vy=app->ball[i].speed_y;
    
    double g=app->gy;        // gravitational acceleration
    double k1=app->k_visc1;
    double k2=app->k_visc2;
    double kroot;            // visc.coefs specifics: k1 - linear R, k2 - linear S=pi*R*R 

    k1=k1*r;
    k2=k2*r*r;
    kroot=sqrt(vx*vx+vy*vy);
    x= x + vx*dt;
    y = y + vy*dt;
    vy = vy*(1.0-k1-k2*kroot) + g*dt;         // gravity+viscosity
    vx = vx*(1.0-k1-k2*kroot);                // viscosity
					      
    app->ball[i].ball_x=x;
    app->ball[i].ball_y=y;
    app->ball[i].speed_x=vx;
    app->ball[i].speed_y=vy;
}

//       FUNCTIONS - REAL DRAWING


void draw(Application* app, int i)  // rendering the i-th ball (render, i-th texture, i-th texture spot, i-th current coords
    {
       filledCircleRGBA(app->rend, app->ball[i].ball_x, app->ball[i].ball_y, app->ball[i].rad, app->ball[i].cr, app->ball[i].cg, app->ball[i].cb,254);
    }


//      FUNCTIONS - FINALIZING

void finish(Application *app)  // TOTAL DESTROYING
{
   SDL_DestroyRenderer(app->rend);
   app->rend = NULL;
   SDL_DestroyWindow(app->window);
   app->window = NULL;

   SDL_Quit();
   IMG_Quit();
}

//       MAIN (IS MAIN)

int main(int argc, char *argv[])
{
    Application *app = (Application*) malloc(sizeof(Application));
    
    int m=0; 
    int l=0;	  
    int j=0;
    int i;
    int xadd,yadd;

    char wname[50];
    const char *winname;  
   
    IMG_Init(IMG_INIT_PNG);

    init_window(app);
    get_rend(app);

    app->deltat=0.04;         // time step (60c^-1)
    app->gy =50.0;                 // gravitation acceleration
    app->k_visc1 = 0.00;         // linear visc.
    app->k_visc2 = 0.00;	  // squaric visc.
    
    SDL_RenderClear(app->rend);
    for (l=0;l<m;l++)
        {
            init_ball(app,l);
        }
    bool keep_window_open = true;
    while(keep_window_open)
    {
        while(SDL_PollEvent(&(app->event)) > 0) // app->event <==> (*app).event
        {
            switch(app->event.type)
            {
                case SDL_QUIT:
                    keep_window_open = false;
                    break;
                 
	        case SDL_MOUSEBUTTONDOWN:
		     SDL_GetMouseState(&xadd, &yadd);  
		     if(m<maxn)
		      {
		          add_ball(app,m, xadd, yadd);    // new ball arising
		          m++;
		      }
		      else
		      break;

                 case SDL_KEYDOWN:                         // keybord event
                         switch(app->event.key.keysym.sym)
                         {
                             case SDLK_UP:                      // arrow up g grows
                                 app->gy += 0.1;
		                         printf("G increased is %f\n", app->gy);	      
				             break;
                             case SDLK_DOWN:                    // and vice versa
                                 app->gy -=0.1;
		                         printf("G decreased is %f\n", app->gy);	      
				             break;
                             default:
                             break;
			              }

            }       
        }
	SDL_SetRenderDrawColor(app->rend,97,205,207,1);
    SDL_RenderClear(app->rend);
//            printf(" DEBUG ERASING %d\t%d\t%d\n",kd,m,destroy);	    
    for (l=0;l<m;l++)
        {
            draw(app, l);
        } 
        for (l=0;l<m;l++)
            {   
                border(app,  l);
                for(j=0;j<l;j++)
                {  
                    collision(app, l,j);
                }
	        shift(app,l);
            }
    SDL_RenderPresent(app->rend);
    SDL_RenderClear(app->rend);

    }
    finish(app);
}
