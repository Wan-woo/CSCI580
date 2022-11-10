/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */
    if (u <= 0.0 || u >= 1.0)    return GZ_FAILURE;
    if (v <= 0.0 || v >= 1.0)    return GZ_FAILURE;
    float x = u * (xs - 1);
    float y = v * (ys - 1);
    float s = x - floor(x);
    float t = y - floor(y);
    GzColor A, B, C, D;
    memcpy((void*)A, (void*)image[(int)floor(x) + (int)floor(y) * xs], sizeof(GzColor));
    memcpy((void*)B, (void*)image[(int)ceil(x) + (int)floor(y) * xs], sizeof(GzColor));
    memcpy((void*)C, (void*)image[(int)ceil(x) + (int)ceil(y) * xs], sizeof(GzColor));
    memcpy((void*)D, (void*)image[(int)floor(x) + (int)ceil(y) * xs], sizeof(GzColor));
    color[0] = s * t * C[0] + (1 - s) * t * D[0] + s * (1 - t) * B[0] + (1 - s) * (1 - t) * A[0];
    color[1] = s * t * C[1] + (1 - s) * t * D[1] + s * (1 - t) * B[1] + (1 - s) * (1 - t) * A[1];
    color[2] = s * t * C[2] + (1 - s) * t * D[2] + s * (1 - t) * B[2] + (1 - s) * (1 - t) * A[2];
    return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
    float x = u * 10;
    float y = v * 10;
    if ((int)floor(x) % 2 == 0 && (int)floor(y) % 2 == 0 || (int)floor(x) % 2 == 1 && (int)floor(y) % 2 == 1)
    {
        color[0] = 0;
        color[1] = 0;
        color[2] = 1;
    }
    else
    {
        color[0] = 1;
        color[1] = 0;
        color[2] = 0;
    }
	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

