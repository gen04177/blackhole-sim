#include <SDL2/SDL.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define MAX_GALAXIES 20
#define DELTAT 0.02
#define QCONS 0.001
#define GALAXYRANGESIZE 0.1
#define GALAXYMINSIZE 0.15
#define COLORBASE 16

double camera_distance = 1.0;
double camera_rot_speed = 0.03;
double camera_zoom_speed = 0.02;

Uint32 blackhole_last_switch = 0;
int blackhole_active = 1;
int blackhole_manual_override = 0;

const Uint32 BH_ACTIVE_TIME = 60000;
const Uint32 BH_INACTIVE_TIME = 300000;



typedef struct
{
  double pos[3];
  double vel[3];
} Star;

typedef struct
{
  int mass;
  int nstars;
  Star *stars;
  SDL_Rect *oldpoints;
  SDL_Rect *newpoints;
  double pos[3];
  double vel[3];
  int galcol;
  Uint8 r, g, b;
} Galaxy;

typedef struct
{
  double pos[3];
  double mass;
} BlackHole;

BlackHole blackhole = { {0, 0, 0}, 1e6 };

Star *free_stars = NULL;
int free_star_count = 0;

typedef struct
{
  double mat[3][3];
  double scale;
  int midx, midy;
  double size;
  double rot_y, rot_x;
  Galaxy galaxies[MAX_GALAXIES];
  int ngalaxies;
  int f_hititerations;
  int step;
  int pscale;
} Universe;

Universe universe;

double
floatRand ()
{
  return (double) rand () / (double) RAND_MAX;
}


void
startover (int width, int height)
{

  for (int i = 0; i < universe.ngalaxies; i++)
    {
      free (universe.galaxies[i].stars);
      free (universe.galaxies[i].oldpoints);
      free (universe.galaxies[i].newpoints);
    }

  universe.step = 0;
  universe.rot_x = 0;
  universe.rot_y = 0;

  universe.scale = fmin (width, height) / 500.0;
  universe.midx = width / 2;
  universe.midy = height / 2;
  universe.pscale = 1;

  universe.ngalaxies = MAX_GALAXIES;
  universe.f_hititerations = 7500;

  for (int i = 0; i < universe.ngalaxies; i++)
    {
      Galaxy *g = &universe.galaxies[i];

      g->nstars = 500 + rand () % 10001;
      g->stars = malloc (sizeof (Star) * g->nstars);
      g->oldpoints = malloc (sizeof (SDL_Rect) * g->nstars);
      g->newpoints = malloc (sizeof (SDL_Rect) * g->nstars);

      g->r = rand () % 256;
      g->g = rand () % 256;
      g->b = rand () % 256;

      g->mass = g->nstars * (1 + floatRand ());

      g->vel[0] = floatRand () * 2.0 - 1.0;
      g->vel[1] = floatRand () * 2.0 - 1.0;
      g->vel[2] = floatRand () * 2.0 - 1.0;

      g->pos[0] =
	(-g->vel[0] * DELTAT * universe.f_hititerations +
	 (floatRand () - 0.5)) * 2.0;
      g->pos[1] =
	(-g->vel[1] * DELTAT * universe.f_hititerations +
	 (floatRand () - 0.5)) * 2.0;
      g->pos[2] =
	(-g->vel[2] * DELTAT * universe.f_hititerations +
	 (floatRand () - 0.5)) * 2.0;

      universe.size = (GALAXYRANGESIZE * floatRand () + GALAXYMINSIZE) * 5.0;

      double w1 = 2.0 * M_PI * floatRand ();
      double w2 = 2.0 * M_PI * floatRand ();
      double sinw1 = sin (w1), cosw1 = cos (w1);
      double sinw2 = sin (w2), cosw2 = cos (w2);

      universe.mat[0][0] = cosw2;
      universe.mat[0][1] = -sinw1 * sinw2;
      universe.mat[0][2] = cosw1 * sinw2;
      universe.mat[1][0] = 0.0;
      universe.mat[1][1] = cosw1;
      universe.mat[1][2] = sinw1;
      universe.mat[2][0] = -sinw2;
      universe.mat[2][1] = -sinw1 * cosw2;
      universe.mat[2][2] = cosw1 * cosw2;

      for (int j = 0; j < g->nstars; j++)
	{
	  Star *st = &g->stars[j];

	  double theta = 2.0 * M_PI * floatRand ();
	  double r = universe.size * sqrt (floatRand ());
	  double z = ((floatRand () - 0.5) * universe.size / 10.0);

	  st->pos[0] =
	    universe.mat[0][0] * r * cos (theta) +
	    universe.mat[1][0] * r * sin (theta) + universe.mat[2][0] * z +
	    g->pos[0];
	  st->pos[1] =
	    universe.mat[0][1] * r * cos (theta) +
	    universe.mat[1][1] * r * sin (theta) + universe.mat[2][1] * z +
	    g->pos[1];
	  st->pos[2] =
	    universe.mat[0][2] * r * cos (theta) +
	    universe.mat[1][2] * r * sin (theta) + universe.mat[2][2] * z +
	    g->pos[2];

	  double v = sqrt (g->mass * QCONS / (r + 0.01));
	  st->vel[0] =
	    (-universe.mat[0][0] * v * sin (theta) +
	     universe.mat[1][0] * v * cos (theta) + g->vel[0]) * DELTAT;
	  st->vel[1] =
	    (-universe.mat[0][1] * v * sin (theta) +
	     universe.mat[1][1] * v * cos (theta) + g->vel[1]) * DELTAT;
	  st->vel[2] =
	    (-universe.mat[0][2] * v * sin (theta) +
	     universe.mat[1][2] * v * cos (theta) + g->vel[2]) * DELTAT;

	  g->oldpoints[j] = (SDL_Rect)
	  {
	  0, 0, universe.pscale, universe.pscale};
	  g->newpoints[j] = (SDL_Rect)
	  {
	  0, 0, universe.pscale, universe.pscale};
	}
    }
}



void
draw_galaxy (SDL_Renderer * renderer)
{
  double cox = cos (universe.rot_y);
  double six = sin (universe.rot_y);
  double cor = cos (universe.rot_x);
  double sir = sin (universe.rot_x);

  for (int i = 0; i < universe.ngalaxies; i++)
    {
      Galaxy *g = &universe.galaxies[i];

      g->pos[0] += g->vel[0];
      g->pos[1] += g->vel[1];
      g->pos[2] += g->vel[2];

      for (int j = 0; j < g->nstars; j++)
	{
	  Star *st = &g->stars[j];
	  SDL_Rect *newp = &g->newpoints[j];

	  st->pos[0] += st->vel[0];
	  st->pos[1] += st->vel[1];
	  st->pos[2] += st->vel[2];

	  newp->x =
	    (int) (((cox * st->pos[0]) -
		    (six * st->pos[2])) * universe.scale * universe.pscale) +
	    universe.midx;
	  newp->y =
	    (int) (((cor * st->pos[1]) -
		    (sir * ((six * st->pos[0]) + (cox * st->pos[2])))) *
		   universe.scale * universe.pscale) + universe.midy;
	  newp->w = newp->h = universe.pscale;

	  double dx = blackhole.pos[0] - st->pos[0];
	  double dy = blackhole.pos[1] - st->pos[1];
	  double dz = blackhole.pos[2] - st->pos[2];
	  double dist = sqrt (dx * dx + dy * dy + dz * dz);

	  Uint8 intensity = 255;
	  if (dist < 0.5)
	    {
	      intensity = (Uint8) (255 * (dist / 0.5));
	    }

	  SDL_SetRenderDrawColor (renderer,
				  g->r * intensity / 255,
				  g->g * intensity / 255,
				  g->b * intensity / 255, 255);

	  SDL_RenderFillRect (renderer, newp);
	}
    }
}


void
handle_collisions ()
{
  for (int i = 0; i < universe.ngalaxies; i++)
    {
      for (int j = i + 1; j < universe.ngalaxies; j++)
	{
	  Galaxy *g1 = &universe.galaxies[i];
	  Galaxy *g2 = &universe.galaxies[j];

	  if (g1->mass == 0 || g2->mass == 0)
	    continue;

	  double dx = g2->pos[0] - g1->pos[0];
	  double dy = g2->pos[1] - g1->pos[1];
	  double dz = g2->pos[2] - g1->pos[2];
	  double dist2 = dx * dx + dy * dy + dz * dz;
	  double min_dist = 0.5;

	  if (dist2 < min_dist * min_dist)
	    {
	      double m1 = g1->mass, m2 = g2->mass;

	      for (int k = 0; k < 3; k++)
		{
		  double v1_new =
		    (g1->vel[k] * (m1 - m2) + 2 * m2 * g2->vel[k]) / (m1 +
								      m2);
		  double v2_new =
		    (g2->vel[k] * (m2 - m1) + 2 * m1 * g1->vel[k]) / (m1 +
								      m2);
		  g1->vel[k] = v1_new;
		  g2->vel[k] = v2_new;
		}

	      Galaxy *big = g1->mass >= g2->mass ? g1 : g2;
	      Galaxy *small = g1->mass >= g2->mass ? g2 : g1;

	      for (int s = 0; s < small->nstars; s++)
		{
		  Star *st = &small->stars[s];

		  st->vel[0] += (floatRand () - 0.5) * 0.05;
		  st->vel[1] += (floatRand () - 0.5) * 0.05;
		  st->vel[2] += (floatRand () - 0.5) * 0.05;

		  st->pos[0] += (big->pos[0] - small->pos[0]);
		  st->pos[1] += (big->pos[1] - small->pos[1]);
		  st->pos[2] += (big->pos[2] - small->pos[2]);
		}

	      big->mass += small->mass;
	      small->mass = 0;
	    }
	}
    }
}


void
update_blackhole_state ()
{
  Uint32 now = SDL_GetTicks ();

  if (!blackhole_manual_override)
    {
      if (blackhole_active && now - blackhole_last_switch >= BH_ACTIVE_TIME)
	{
	  blackhole_active = 0;	
	  blackhole_last_switch = now;
	}
      else if (!blackhole_active
	       && now - blackhole_last_switch >= BH_INACTIVE_TIME)
	{
	  blackhole_active = 1;	
	  blackhole_last_switch = now;
	  startover (universe.midx * 2, universe.midy * 2);
	}
    }
}


void
handle_blackhole ()
{
  if (!blackhole_active)
    return;

  for (int i = 0; i < universe.ngalaxies; i++)
    {
      Galaxy *g = &universe.galaxies[i];
      if (g->mass == 0)
	continue;

      for (int j = 0; j < g->nstars; j++)
	{
	  Star *st = &g->stars[j];

	  double dx = blackhole.pos[0] - st->pos[0];
	  double dy = blackhole.pos[1] - st->pos[1];
	  double dz = blackhole.pos[2] - st->pos[2];
	  double dist2 = dx * dx + dy * dy + dz * dz;
	  double dist = sqrt (dist2);

	  double acc = blackhole.mass * 0.01 / (dist2 + 0.0001);

	  st->vel[0] += acc * dx / dist * DELTAT;
	  st->vel[1] += acc * dy / dist * DELTAT;
	  st->vel[2] += acc * dz / dist * DELTAT;

	  if (dist < 0.05)
	    {
	      st->pos[0] = st->pos[1] = st->pos[2] = blackhole.pos[0];
	      st->vel[0] = st->vel[1] = st->vel[2] = 0;
	    }
	}
    }
}


void
reset_blackhole ()
{
  blackhole_active = 1;
  blackhole_last_switch = SDL_GetTicks ();
  blackhole_manual_override = 0;
}

void 
printHelp(void)
{

    printf("\n\nblackhole-sim by gen04177 - v0.1 - Running with SDL %d.%d.%d (compiled with %d.%d.%d)\n", SDL_MAJOR_VERSION, SDL_MINOR_VERSION, SDL_PATCHLEVEL, SDL_COMPILEDVERSION >> 24, (SDL_COMPILEDVERSION >> 16) & 0xFF, SDL_COMPILEDVERSION & 0xFFFF);

}

int
main (int argc, char *argv[])
{
        printHelp();
	
  srand (time (NULL));
  if (SDL_Init (SDL_INIT_VIDEO | SDL_INIT_GAMECONTROLLER) < 0)
    {
      fprintf (stderr, "Could not initialize SDL: %s\n", SDL_GetError ());
      return 1;
    }

  SDL_GameController *controller = NULL;

  if (SDL_NumJoysticks () > 0)
    {
      if (SDL_IsGameController (0))
	{
	  controller = SDL_GameControllerOpen (0);
	  if (!controller)
	    {
	      printf ("Failed to open controller: %s\n", SDL_GetError ());
	    }
	}
    }

  int width = 1920, height = 1080;
  SDL_Window *window =
    SDL_CreateWindow ("fake blackhole-sim", SDL_WINDOWPOS_CENTERED,
		      SDL_WINDOWPOS_CENTERED, width, height,
		      SDL_WINDOW_SHOWN);
  SDL_Renderer *renderer =
    SDL_CreateRenderer (window, -1, SDL_RENDERER_ACCELERATED);

  startover (width, height);

  int running = 1;
  SDL_Event event;
  while (running)
    {
      while (SDL_PollEvent (&event))
	{
	  if (event.type == SDL_QUIT)
	    running = 0;

	  if (event.type == SDL_CONTROLLERBUTTONDOWN)
	    {
	      if (event.cbutton.button == SDL_CONTROLLER_BUTTON_RIGHTSHOULDER)
		{
		  startover (width, height);
		}
	    }

	  if (event.type == SDL_CONTROLLERBUTTONDOWN)
	    {
	      if (event.cbutton.button == SDL_CONTROLLER_BUTTON_LEFTSHOULDER)
		{
		  reset_blackhole ();
		}
	    }


	}

      if (controller)
	{

	  Sint16 lx =
	    SDL_GameControllerGetAxis (controller, SDL_CONTROLLER_AXIS_LEFTX);
	  Sint16 ly =
	    SDL_GameControllerGetAxis (controller, SDL_CONTROLLER_AXIS_LEFTY);

	  double nx = lx / 32768.0;
	  double ny = ly / 32768.0;

	  if (fabs (nx) < 0.1)
	    nx = 0;
	  if (fabs (ny) < 0.1)
	    ny = 0;

	  universe.rot_y += nx * camera_rot_speed;

	  camera_distance += ny * camera_zoom_speed;

	  if (camera_distance < 0.2)
	    camera_distance = 0.2;

	  if (camera_distance > 100.0)
	    camera_distance = 100.0;

	  universe.scale = camera_distance * (fmin (width, height) / 500.0);
	}


      update_blackhole_state ();
      handle_blackhole ();
      handle_collisions ();

      SDL_SetRenderDrawColor (renderer, 0, 0, 0, 255);
      SDL_RenderClear (renderer);

      draw_galaxy (renderer);

      SDL_RenderPresent (renderer);

      SDL_Delay (16);
    }

  SDL_DestroyRenderer (renderer);
  SDL_DestroyWindow (window);
  if (controller)
    {
      SDL_GameControllerClose (controller);
    }
  SDL_Quit ();
  return 0;
}
