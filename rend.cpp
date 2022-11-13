/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"raytracing.h"

#define PI (float) 3.14159265358979323846
/* HW1 methods: copy here the methods from HW1 */
#ifndef XYZ2NUM
#define X	0
#define Y	1
#define Z	2
#endif
GzMatrix IDENTIFY = { { 1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
GzTextureIndex uv[2];
GzCoord E = { 0, 0, -1 };

GzRender::DDA::DDA(float* start, float* end)
{
	this->start = start;
	this->end = end;
	this->current[X] = start[X];
	this->current[Y] = start[Y];
	this->current[Z] = start[Z];
}
void GzRender::DDA::initEdge()
{
	this->slopeX = (end[X] - start[X]) / (end[Y] - start[Y]);
	this->slopeZ = (end[Z] - start[Z]) / (end[Y] - start[Y]);

	float deltaY = ceil(start[Y]) - start[Y];
	current[X] = start[X] + deltaY * slopeX;
	current[Y] = start[Y] + deltaY;
	current[Z] = start[Z] + deltaY * slopeZ;

}
void GzRender::DDA::iterateSpan(GzRender* render, float abcd[5][4])
{
	std::array<float, 3> point;

	while (current[X] < end[X])//only include left
	{

		point = { current[X], current[Y], current[Z] };
		GzColor color;
		//result.push_back(coord);
		float Vz = point[2] / (INT_MAX - point[2]);

		if (render->interp_mode == GZ_FLAT)
		{
			for (int i = 0; i < 3; i++)
			{
				color[i] = render->flatcolor[i];
			}
		}
		else if (render->interp_mode == GZ_COLOR)
		{
			GzTextureIndex midUV;
			for (int i = 3; i < 5; i++)
			{
				midUV[i - 3] = (Vz + 1) * ((abcd[i][0] * point[0] + abcd[i][1] * point[1] + abcd[i][3]) / (-1 * abcd[i][2]));
			}
			GzColor K = { 1,1,1 };
			if (render->tex_fun != 0)
			{
				render->tex_fun(midUV[0], midUV[1], K);
			}

			for (int i = 0; i < 3; i++)
			{
				color[i] = K[i] * ((abcd[i][0] * point[0] + abcd[i][1] * point[1] + abcd[i][3]) / (-1 * abcd[i][2]));
			}
		}
		else if (render->interp_mode == GZ_NORMALS)
		{
			GzCoord midNorm;
			
			for (int i = 0; i < 3; i++)
			{
				midNorm[i] = (abcd[i][0] * point[0] + abcd[i][1] * point[1] + abcd[i][3]) / (-1 * abcd[i][2]);
			}
			render->normalize(midNorm, midNorm);
			
			GzTextureIndex midUV;
			for (int i = 3; i < 5; i++)
			{
				
				midUV[i - 3] = (Vz + 1) * ((abcd[i][0] * point[0] + abcd[i][1] * point[1] + abcd[i][3]) / (-1 * abcd[i][2]));
			}
			GzColor K;
			
			if (render->tex_fun != 0)
			{
				render->tex_fun(midUV[0], midUV[1], K);
				for (int i = 0; i < 3; i++)
				{
					render->Kd[i] = K[i];
					render->Ka[i] = K[i];
				}
			}
			render->computeColor(midNorm, color);
		}
		GzDepth z;
		render->GzGet(point[0], point[1], NULL, NULL, NULL, NULL, &z);
		if ((int)point[2] < z)
		{
			render->GzPut(point[0], point[1], render->ctoi(color[0]), render->ctoi(color[1]), render->ctoi(color[2]), 0, point[2]);
		}
		current[X] += 1;
		current[Z] += slopeZ;

	}


}
void GzRender::DDA::initSpan()
{
	this->slopeZ = (end[Z] - start[Z]) / (end[X] - start[X]);
	float deltaX = ceil(start[X]) - start[X];
	current[X] = start[X] + deltaX;
	current[Y] = start[Y];
	current[Z] = start[Z] + deltaX * slopeZ;

}
void GzRender::DDA::nextY()
{
	current[X] += slopeX;
	current[Y] += 1;
	current[Z] += slopeZ;
}

GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */
	if (xRes > MAXXRES || yRes > MAXYRES)	exit(GZ_FAILURE);

	xres = xRes;
	yres = yRes;
	framebuffer = new char[3 * xres * yres];
	pixelbuffer = new GzPixel[xres * yres];
	trianglebuffer = new GzTri[MAX_TRIANGLES];

	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/

	GzMatrix sp = {
		{xres / 2, 0, 0, xres / 2},
		{0, yres / -2, 0, yres / 2},
		{0, 0, INT_MAX, 0},
		{0, 0, 0, 1} };
	std::copy(&sp[0][0], &sp[0][0] + 16, &Xsp[0][0]);

	m_camera.FOV = DEFAULT_FOV;
	float lookat[] = { 0, 0, 0 };
	std::copy(lookat, lookat + 3, m_camera.lookat);
	float position[] = { DEFAULT_IM_X, DEFAULT_IM_Y, DEFAULT_IM_Z };
	std::copy(position, position + 3, m_camera.position);
	float worldUp[] = { 0, 1, 0 };
	std::copy(worldUp, worldUp + 3, m_camera.worldup);

	matlevel = -1;
	triIndex = 0;


}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	if (framebuffer != NULL)
	{
		delete[] framebuffer;
	}
	if (pixelbuffer != NULL)
	{
		delete[] pixelbuffer;
	}



}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */

	for (int i = 0; i < xres * yres; i++)
	{
		pixelbuffer[i] = {
			1000,
			1000,
			1000,
			4095,
			INT_MAX
		};
	}

	for (int i = 0; i < MAX_TRIANGLES; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				trianglebuffer[triIndex].vertices[j][k] = 0;
				trianglebuffer[triIndex].imageVerts[j][k] = 0;
			}
		}
	}



	return GZ_SUCCESS;
}


int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	float oneOverD = tan(m_camera.FOV / 2 * PI / 180);
	GzMatrix pi = { { 1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, oneOverD, 0}, {0, 0, oneOverD, 1} };
	std::copy(&pi[0][0], &pi[0][0] + 16, &m_camera.Xpi[0][0]);

	GzCoord cameraZ;
	GzCoord cameraY;
	GzCoord cameraX;
	minus(m_camera.lookat, m_camera.position, cameraZ);
	normalize(cameraZ, cameraZ);

	GzCoord tmpUp;
	scalarProduct(cameraZ, dotProduct(m_camera.worldup, cameraZ), tmpUp);
	minus(m_camera.worldup, tmpUp, cameraY);
	normalize(cameraY, cameraY);

	crossProduct(cameraY, cameraZ, cameraX);

	GzMatrix iw = {
		{cameraX[0], cameraX[1], cameraX[2], (-1) * dotProduct(cameraX, m_camera.position)},
		{cameraY[0], cameraY[1], cameraY[2], (-1) * dotProduct(cameraY, m_camera.position)},
		{cameraZ[0], cameraZ[1], cameraZ[2], (-1) * dotProduct(cameraZ, m_camera.position)},
		{0, 0, 0, 1} };
	std::copy(&iw[0][0], &iw[0][0] + 16, &m_camera.Xiw[0][0]);

	GzPushMatrix(Xsp, IDENTIFY);
	GzPushMatrix(m_camera.Xpi, IDENTIFY);
	GzPushMatrix(m_camera.Xiw, m_camera.Xiw);
	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/	m_camera = camera;

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix, GzMatrix norm)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	GzMatrix pureNorm;
	float scale = sqrt(norm[0][0] * norm[0][0] + norm[0][1] * norm[0][1] + norm[0][2] * norm[0][2]);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 3 && j == 3)
			{
				pureNorm[i][j] = 1;
			}
			else if (i == 3 || j == 3)
			{
				pureNorm[i][j] = 0;
			}
			else
			{
				pureNorm[i][j] = norm[i][j] / scale;
			}
		}
	}
	matlevel++;
	if (matlevel >= MATLEVELS)
	{
		exit(GZ_FAILURE);
	}

	if (matlevel == 0)
	{
		std::copy(&matrix[0][0], &matrix[0][0] + 16, &Ximage[matlevel][0][0]);
		std::copy(&pureNorm[0][0], &pureNorm[0][0] + 16, &Xnorm[matlevel][0][0]);
	}
	else
	{
		GzMatrix res;
		matrixProduct(Ximage[matlevel - 1], matrix, res);
		std::copy(&res[0][0], &res[0][0] + 16, &Ximage[matlevel][0][0]);

		matrixProduct(Xnorm[matlevel - 1], pureNorm, res);
		std::copy(&res[0][0], &res[0][0] + 16, &Xnorm[matlevel][0][0]);
	}


	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel < 0)
	{
		exit(GZ_FAILURE);
	}
	matlevel--;
	return GZ_SUCCESS;
}
int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */
	if (i < xres && j < yres && i >= 0 && j >= 0)
	{
		pixelbuffer[ARRAY(i, j)] = {
			clamp12BitRGB(r),
			clamp12BitRGB(g),
			clamp12BitRGB(b),
			a,
			z
		};
	}

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i < xres && j < yres && i >= 0 && j >= 0)
	{
		GzPixel pixel = pixelbuffer[ARRAY(i, j)];
		if (r != NULL)	*r = pixel.red;
		if (g != NULL)	*g = pixel.green;
		if (b != NULL)	*b = pixel.blue;
		if (a != NULL)	*a = pixel.alpha;
		if (z != NULL)	*z = pixel.z;
	}
	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	char header[30];
	sprintf(header, "P6 %d %d 255\r", xres, yres);
	fputs(header, outfile);
	char* rgb = new char[3 * xres * yres];
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			int pos = ARRAY(i, j);
			GzPixel pixel = pixelbuffer[pos];
			fputc(pixel.red >> 4, outfile);
			fputc(pixel.green >> 4, outfile);
			fputc(pixel.blue >> 4, outfile);


		}
	}

	fflush(outfile);
	delete[] rgb;
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			int pos = ARRAY(i, j);
			GzPixel pixel = pixelbuffer[pos];
			framebuffer[pos * 3 + 0] = pixel.blue >> 4;
			framebuffer[pos * 3 + 1] = pixel.green >> 4;
			framebuffer[pos * 3 + 2] = pixel.red >> 4;

		}
	}
	return GZ_SUCCESS;
}



/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/
	for (int i = 0; i < numAttributes; i++)
	{
		switch (nameList[i])
		{
		case GZ_RGB_COLOR:
		{
			flatcolor[0] = clampFloat(((float*)valueList[i])[0]);
			flatcolor[1] = clampFloat(((float*)valueList[i])[1]);
			flatcolor[2] = clampFloat(((float*)valueList[i])[2]);
			break;
		}
		case GZ_AMBIENT_LIGHT:
		{
			ambientlight = *(GzLight*)valueList[i];
			break;
		}
		case GZ_DIRECTIONAL_LIGHT:
		{
			numlights = numAttributes;
			lights[i] = *(GzLight*)valueList[i];
			break;
		}
		case GZ_AMBIENT_COEFFICIENT:
		{
			Ka[0] = ((float*)valueList[i])[0];
			Ka[1] = ((float*)valueList[i])[1];
			Ka[2] = ((float*)valueList[i])[2];
			break;
		}
		case GZ_DIFFUSE_COEFFICIENT:
		{
			Kd[0] = ((float*)valueList[i])[0];
			Kd[1] = ((float*)valueList[i])[1];
			Kd[2] = ((float*)valueList[i])[2];
			break;
		}
		case GZ_SPECULAR_COEFFICIENT:
		{
			Ks[0] = ((float*)valueList[i])[0];
			Ks[1] = ((float*)valueList[i])[1];
			Ks[2] = ((float*)valueList[i])[2];
			break;
		}
		case GZ_DISTRIBUTION_COEFFICIENT:
		{
			spec = *(float*)valueList[i];
			break;
		}
		case GZ_INTERPOLATE:
		{
			interp_mode = *(int*)valueList[i];
			break;
		}
		case GZ_TEXTURE_MAP:
		{
			tex_fun = *(GzTexture)valueList[i];
			break;
		}
		default:
			break;
		}



	}
	return GZ_SUCCESS;
}
void GzRender::computeColor(GzCoord norm, GzCoord res)
{
	GzCoord specularColor = { 0.0, 0.0, 0.0 };
	GzCoord diffuseColor = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < numlights; i++)
	{
		GzLight light = lights[i];

		GzCoord tmp;
		GzCoord R;
		float NE = dotProduct(norm, E);
		float NL = dotProduct(norm, light.direction);
		if (NE < 0 && NL < 0)
		{
			scalarProduct(norm, -1, norm);
			NL = dotProduct(norm, light.direction);
		}
		else if (NE > 0 && NL > 0)
		{
			//nothing
		}
		else
		{
			continue;
		}


		scalarProduct(norm, 2 * NL, tmp);
		minus(tmp, light.direction, R);


		float RES = clampFloat(pow(dotProduct(R, E), spec));
		scalarProduct(light.color, RES, tmp);
		add(specularColor, tmp, specularColor);

		scalarProduct(light.color, NL, tmp);
		add(diffuseColor, tmp, diffuseColor);

	}
	for (int i = 0; i < 3; i++)
	{
		specularColor[i] *= Ks[i];
		diffuseColor[i] *= Kd[i];
		res[i] = specularColor[i] + diffuseColor[i] + Ka[i] * ambientlight.color[i];
	}

}

int GzRender::GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/
	GzCoord triNorm[3];
	for (int i = 0; i < numParts; i++)
	{
		switch (nameList[i])
		{
		case GZ_POSITION:
		{
			GzCoord triCoord[3];
			GzCoord rawCoord[3];
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					rawCoord[j][k] = ((GzCoord*)valueList[i])[j][k];
				}
			}

			//transform from model space to screen space
			for (int coord = 0; coord < 3; coord++)
			{
				boolean inRangeZ = transform(Ximage[matlevel], rawCoord[coord], triCoord[coord]);
				if (!inRangeZ)	return GZ_SUCCESS;

				float Vz = triCoord[coord][2] / (INT_MAX - triCoord[coord][2]);
				uv[coord][0] /= (Vz + 1);
				uv[coord][1] /= (Vz + 1);
			}

			//sort
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3 - j - 1; k++)
				{
					if ((triCoord[k][Y] > triCoord[k + 1][Y])
						|| ((triCoord[k][Y] == triCoord[k + 1][Y]) && (triCoord[k][X] > triCoord[k + 1][X])))
					{
						for (int m = 0; m < 3; m++)
						{
							float tmp = triCoord[k][m];
							triCoord[k][m] = triCoord[k + 1][m];
							triCoord[k + 1][m] = tmp;
						}
					}
				}
			}

			//trianglebuffer
			if (triIndex < MAX_TRIANGLES)
			{
				GzCoord imageCoord[3];
				for (int coord = 0; coord < 3; coord++)
				{
					transform(Xnorm[matlevel], rawCoord[coord], imageCoord[coord]);
				}
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						trianglebuffer[triIndex].vertices[j][k] = triCoord[j][k];
						trianglebuffer[triIndex].imageVerts[j][k] = imageCoord[j][k];
					}
				}
			}
			break;
		}

		case GZ_NORMAL:
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					triNorm[j][k] = ((GzCoord*)valueList[i])[j][k];
				}
			}

			//transform normals from model space to image space
			for (int coord = 0; coord < 3; coord++)
			{
				transform(Xnorm[matlevel], triNorm[coord], triNorm[coord]);
			}

			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					trianglebuffer[triIndex].normals[j][k] = triNorm[j][k];
				}
			}

			break;
		}
		case GZ_TEXTURE_INDEX:
		{

			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 2; k++)
				{

					uv[j][k] = ((GzTextureIndex*)valueList[i])[j][k];

				}
			}

			break;
		}
		case GZ_NULL_TOKEN:
		{
			break;
		}



		default:
			break;
		}
	}

	//get ready for next triangle
	triIndex++;



	return GZ_SUCCESS;
}


int GzRender::GzRaytracing()
{
	ray.origin[X] = m_camera.position[0];
	ray.origin[Y] = m_camera.position[1];
	ray.origin[Z] = m_camera.position[2];

	for (int i = 0; i < triIndex; i++)
	{
		GzTri tri = trianglebuffer[i];

		//compute triangle center and then ray direction to it
		/*ray.direction[0] = (tri.vertices[0][X]+ tri.vertices[1][X] + tri.vertices[2][X]) / 3.0 - ray.origin[X];
		ray.direction[1] = (tri.vertices[0][Y] + tri.vertices[1][Y] + tri.vertices[2][Y]) / 3.0 - ray.origin[Y];
		ray.direction[2] = (tri.vertices[0][Z] + tri.vertices[1][Z] + tri.vertices[2][Z]) / 3.0 - ray.origin[Z];*/
		GzCoord centroid;
		centroid[X] = ((tri.imageVerts[0][X] + tri.imageVerts[1][X] + tri.imageVerts[2][X]) / 3.0);
		centroid[Y] = ((tri.imageVerts[0][Y] + tri.imageVerts[1][Y] + tri.imageVerts[2][Y]) / 3.0);
		centroid[Z] = ((tri.imageVerts[0][Z] + tri.imageVerts[1][Z] + tri.imageVerts[2][Z]) / 3.0);
		minus(centroid, ray.origin, ray.direction);
		normalize(ray.direction, ray.direction);

		//shoot ray through pixels in triangle
		int result = RayIntersection(tri);

		//compute colors if hit (we do flat for now)
		if (result == GZ_FAILURE)
		{
			return GZ_FAILURE;
		}
		else
		{
			Rasterize(tri);
		}

	}

	return GZ_SUCCESS;
}

int GzRender::PointAtTValue(float t, GzCoord coord)
{
	//origin + t * dir
	GzCoord scaled;
	scalarProduct(ray.direction, t, scaled);
	add(ray.origin, scaled, coord);

	return GZ_SUCCESS;
}

int GzRender::RayIntersection(GzTri triangle)
{
	//compute triangle (plane) normal using: crossProduct(1, 2, result)
	GzCoord edge1, edge2, edge3;
	minus(triangle.imageVerts[1], triangle.imageVerts[0], edge1);
	minus(triangle.imageVerts[2], triangle.imageVerts[1], edge2);
	minus(triangle.imageVerts[0], triangle.imageVerts[2], edge3);

	GzCoord planeNorm;
	crossProduct(edge1, edge2, planeNorm);
	normalize(planeNorm, planeNorm);

	//check if ray and triangle are parallel
	float dir = dotProduct(planeNorm, ray.direction);
	if (dir == 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	//compute D then t, ensure t >= 0
	//float t = (D - dot(N, orig)) / dot(N, dir); 
	float dCoeff = dotProduct(planeNorm, triangle.imageVerts[2]);
	float tValue = (dCoeff - dotProduct(planeNorm, ray.origin)) / (dotProduct(planeNorm, ray.direction));
	if (tValue < 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	//check if intersection point is inside triangle
	GzCoord Pvalue, coord1, coord2, coord3, cross1, cross2, cross3;
	PointAtTValue(tValue, Pvalue);

	minus(Pvalue, triangle.imageVerts[0], coord1);
	minus(Pvalue, triangle.imageVerts[1], coord2);
	minus(Pvalue, triangle.imageVerts[2], coord3);
	crossProduct(edge1, coord1, cross1);
	crossProduct(edge2, coord2, cross2);
	crossProduct(edge3, coord3, cross3);
	
	float dot1, dot2, dot3;
	dot1 = dotProduct(cross1, planeNorm);
	dot2 = dotProduct(cross2, planeNorm);
	dot3 = dotProduct(cross3, planeNorm);

	if (dot1 > 0 && dot2 > 0 && dot3 > 0)
	{
		return GZ_SUCCESS;	 //ray hits triangle :)
	}
	
	return GZ_FAILURE;
}

int GzRender::Rasterize(GzTri triangle)
{
	GzCoord triCoord[3], triNorm[3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			triCoord[i][j] = triangle.vertices[i][j];
			triNorm[i][j] = triangle.normals[i][j];
		}
	}

	//compute ABCD for interpolation
	float abcd[5][4];

	if (interp_mode == GZ_FLAT)
	{
		computeColor(triNorm[0], flatcolor);
	}
	else if (interp_mode == GZ_COLOR)
	{

		GzColor vertexRGB[3];
		for (int vol = 0; vol < 3; vol++)
		{
			Ks[vol] = 1.0;
			Kd[vol] = 1.0;
			Ka[vol] = 1.0;
		}
		for (int vertex = 0; vertex < 3; vertex++)// for every vertex
		{
			computeColor(triNorm[vertex], vertexRGB[vertex]);
		}


		for (int volume = 0; volume < 3; volume++)
		{
			GzCoord e01;
			GzCoord e12;
			minus(triCoord[0], triCoord[1], e01);
			e01[2] = vertexRGB[0][volume] - vertexRGB[1][volume];
			minus(triCoord[1], triCoord[2], e12);
			e12[2] = vertexRGB[1][volume] - vertexRGB[2][volume];
			crossProduct(e01, e12, abcd[volume]);

			abcd[volume][3] = (-1) * (abcd[volume][0] * triCoord[0][0] + abcd[volume][1] * triCoord[0][1]
				+ abcd[volume][2] * vertexRGB[0][volume]);

		}


	}
	else if (interp_mode == GZ_NORMALS)
	{

		for (int volume = 0; volume < 3; volume++)
		{
			GzCoord e01;
			GzCoord e12;
			minus(triCoord[0], triCoord[1], e01);
			e01[2] = triNorm[0][volume] - triNorm[1][volume];
			minus(triCoord[1], triCoord[2], e12);
			e12[2] = triNorm[1][volume] - triNorm[2][volume];
			crossProduct(e01, e12, abcd[volume]);

			abcd[volume][3] = (-1) * (abcd[volume][0] * triCoord[0][0] + abcd[volume][1] * triCoord[0][1]
				+ abcd[volume][2] * triNorm[0][volume]);

		}
	}

	for (int volume = 3; volume < 5; volume++)//u and v
	{

		GzCoord e01;
		GzCoord e12;
		minus(triCoord[0], triCoord[1], e01);
		e01[2] = uv[0][volume - 3] - uv[1][volume - 3];
		minus(triCoord[1], triCoord[2], e12);
		e12[2] = uv[1][volume - 3] - uv[2][volume - 3];
		crossProduct(e01, e12, abcd[volume]);

		abcd[volume][3] = (-1) * (abcd[volume][0] * triCoord[0][0] + abcd[volume][1] * triCoord[0][1]
			+ abcd[volume][2] * uv[0][volume - 3]);

	}

	DDA* dda01 = new DDA(triCoord[0], triCoord[1]);
	DDA* dda12 = new DDA(triCoord[1], triCoord[2]);
	DDA* dda02 = new DDA(triCoord[0], triCoord[2]);

	//decide left/right
	if (triCoord[0][Y] == triCoord[1][Y])//top horizon
	{
		//advance
		dda12->initEdge();
		dda02->initEdge();

		dda01->leftOrTop = true;
		dda12->leftOrTop = false;
		dda02->leftOrTop = true;
	}
	else if (triCoord[1][Y] == triCoord[2][Y])//bottom horizon
	{
		//advance
		dda01->initEdge();
		dda02->initEdge();

		dda01->leftOrTop = true;
		dda12->leftOrTop = false;
		dda02->leftOrTop = false;
	}
	else
	{
		//advance
		dda01->initEdge();
		dda12->initEdge();
		dda02->initEdge();

		if (dda01->slopeX > dda02->slopeX)
		{
			dda01->leftOrTop = false;
			dda12->leftOrTop = false;
			dda02->leftOrTop = true;
		}
		else
		{
			dda01->leftOrTop = true;
			dda12->leftOrTop = true;
			dda02->leftOrTop = false;
		}

	}

	std::vector<std::array<float, 3>> result;
	std::vector<std::array<float, 3>> resultColor;
	//first half

	while (dda01->current[Y] < dda01->end[Y])
	{
		DDA* ddaSpan;
		if (dda01->leftOrTop)
		{
			ddaSpan = new DDA(dda01->current, dda02->current);

		}
		else
		{
			ddaSpan = new DDA(dda02->current, dda01->current);
		}
		ddaSpan->initSpan();
		ddaSpan->iterateSpan(this, abcd);
		dda01->nextY();
		dda02->nextY();
		delete ddaSpan;

	}

	//second half

	while (dda02->current[Y] < dda02->end[Y])
	{
		DDA* ddaSpan;
		if (dda02->leftOrTop)
		{
			ddaSpan = new DDA(dda02->current, dda12->current);

		}
		else
		{
			ddaSpan = new DDA(dda12->current, dda02->current);
		}
		ddaSpan->initSpan();
		ddaSpan->iterateSpan(this, abcd);
		dda02->nextY();
		dda12->nextY();
		delete ddaSpan;

	}

	/*
	char* str = new char[80];
	sprintf(str, "%f, %f, %f\n", triCoord[0][1], triCoord[1][1], triCoord[2][2]);
	OutputDebugStringA((LPCSTR)str);
	*/

	delete dda01;
	delete dda02;
	delete dda12;

	return GZ_SUCCESS;
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	float radius = degree / 180 * PI;
	GzMatrix rotateX = {
	{1, 0, 0, 0},
	{0, cos(radius), -1 * sin(radius), 0},
	{0, sin(radius), cos(radius), 0},
	{0, 0, 0, 1} };
	std::copy(&rotateX[0][0], &rotateX[0][0] + 16, &mat[0][0]);

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	float radius = degree / 180 * PI;
	GzMatrix rotateY = {
	{cos(radius), 0, sin(radius), 0},
	{0, 1, 0, 0},
	{-1 * sin(radius), 0, cos(radius), 0},
	{0, 0, 0, 1} };
	std::copy(&rotateY[0][0], &rotateY[0][0] + 16, &mat[0][0]);

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	float radius = degree / 180 * PI;
	GzMatrix rotateZ = {
	{cos(radius), -1 * sin(radius), 0, 0},
	{sin(radius),cos(radius), 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1} };
	std::copy(&rotateZ[0][0], &rotateZ[0][0] + 16, &mat[0][0]);
	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/

	GzMatrix transMatrix = {
	{1, 0, 0, translate[0]},
	{0, 1, 0, translate[1]},
	{0, 0, 1, translate[2]},
	{0, 0, 0, 1} };
	std::copy(&transMatrix[0][0], &transMatrix[0][0] + 16, &mat[0][0]);
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	GzMatrix scaleMatrix = {
	{scale[0], 0, 0, 0},
	{0, scale[1], 0, 0},
	{0, 0, scale[2], 0},
	{0, 0, 0, 1} };
	std::copy(&scaleMatrix[0][0], &scaleMatrix[0][0] + 16, &mat[0][0]);
	return GZ_SUCCESS;
}






