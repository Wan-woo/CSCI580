#include	"gz.h"
#ifndef GZRENDER_
#define GZRENDER_


/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */
#define	MAX_TRIANGLES	1000		/* how many triangles allowed */

class GzRender{			/* define a renderer */
  

public:
	unsigned short	xres;
	unsigned short	yres;
	GzPixel		*pixelbuffer;		/* frame buffer array */
	char* framebuffer;

	GzCamera		m_camera;
	short		    matlevel;	        /* top of stack - current xform */
	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
	GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
	GzColor		flatcolor;          /* color state for flat shaded triangles */
	int			interp_mode;
	int			numlights;
	GzLight		lights[MAX_LIGHTS];
	GzLight		ambientlight;
	GzColor		Ka, Kd, Ks;
	float		    spec;		/* specular power */
	GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
	
	/* Raytracing */
	GzTri*		trianglebuffer;
	int			triIndex;
	GzRay		ray;

  	// Constructors
	GzRender(int xRes, int yRes);
	~GzRender();

	// Function declaration

	// HW1: Display methods
	int GzDefault();
	int GzBeginRender();
	int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	int GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth	*z);

	int GzFlushDisplay2File(FILE* outfile);
	int GzFlushDisplay2FrameBuffer();

	// HW2: Render methods
	int GzPutAttribute(int numAttributes, GzToken *nameList, GzPointer *valueList);
	int GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList);

	// HW3
	int GzPutCamera(GzCamera camera);
	int GzPushMatrix(GzMatrix	matrix);
	int GzPopMatrix();

	// Extra methods: NOT part of API - just for general assistance */
	inline int ARRAY(int x, int y){return (x+y*xres);}	/* simplify fbuf indexing */
	inline short	ctoi(float color) {return(short)((int)(color * ((1 << 12) - 1)));}		/* convert float color to GzIntensity short */

	// Object Translation
	int GzRotXMat(float degree, GzMatrix mat);
	int GzRotYMat(float degree, GzMatrix mat);
	int GzRotZMat(float degree, GzMatrix mat);
	int GzTrxMat(GzCoord translate, GzMatrix mat);
	int GzScaleMat(GzCoord scale, GzMatrix mat);

	// User add function compute color at a normal
	void GzComputeColor(GzCoord norm, GzColor color);
	// User add function compute (R,G,B) plane coeff or (N1,N2,N3) coeff
	void GzComputePlane(GzCoord* vertices, GzCoord* coeff);
	void GzComputePlane(GzCoord* vertices, GzCoord* coeff, float (*plane)[4]);
	float planes[3][4];


	//Ray tracing
	int GzRaytracing();
	int PointAtTValue(float t, GzCoord coord);
	int RayIntersection(GzTri triangle);
	int Rasterize(GzTri triangle);
	int AssignTriangleToPixel();
	bool IsPixelInTriangle(int i, int j, GzPixel pixel, GzTri triangle, GzCoord hitPoint);
	int ComputeRaycastColor();
	bool IsPixelAShadow(GzPixel pixel, GzLight light);
	bool GzFindFrontestIntersection(GzTri*& triangle, GzCoord intersection);
	float GzRender::dotProduct(GzCoord s1, GzCoord s2)
	{
		float res = 0;
		for (int i = 0; i < 3; i++)
		{
			res += s1[i] * s2[i];
		}
		return res;
	}
	void GzRender::crossProduct(GzCoord s1, GzCoord s2, GzCoord res)
	{
		res[0] = s1[1] * s2[2] - s1[2] * s2[1];
		res[1] = s1[2] * s2[0] - s1[0] * s2[2];
		res[2] = s1[0] * s2[1] - s1[1] * s2[0];
	}
	void GzRender::normalize(GzCoord s, GzCoord res)
	{
		float value = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
		for (int i = 0; i < 3; i++)
		{
			res[i] = s[i] / value;
		}
	}
	void GzRender::minus(GzCoord s1, GzCoord s2, GzCoord res)
	{
		for (int i = 0; i < 3; i++)
		{
			res[i] = s1[i] - s2[i];
		}
	}

};
#endif