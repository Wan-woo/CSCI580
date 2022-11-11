#include "stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include "raytracing.h"

int Ray::PointAtTValue(float t, GzCoord coord)
{
	//origin + t * dir
	GzCoord scaled;
	scalarProduct(direction, t, scaled);
	add(origin, scaled, coord);

	return GZ_SUCCESS;
}

int Triangle::RayIntersection(Ray ray)
{
	//compute triangle (plane) normal using: crossProduct(1, 2, result)
	GzCoord edge1, edge2, edge3;
	minus(vertices[1], vertices[0], edge1);
	minus(vertices[2], vertices[1], edge2);
	minus(vertices[0], vertices[2], edge3);

	GzCoord planeNorm;
	crossProduct(edge1, edge2, planeNorm);

	//check if ray and triangle are parallel
	float dir = dotProduct(planeNorm, ray.direction);
	if (dir == 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	//compute D then t, ensure t >= 0
	float dCoeff = -1 * dotProduct(planeNorm, vertices[0]);
	float tValue = -1 * (dotProduct(planeNorm, ray.origin) + dCoeff);
	if (tValue < 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	//check if intersection point is inside triangle
	GzCoord point, pEdge, result;
	ray.PointAtTValue(tValue, point);
	minus(point, vertices[0], pEdge);
	crossProduct(edge1, pEdge, result);
	if (dotProduct(planeNorm, result) < 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	minus(point, vertices[1], pEdge);
	crossProduct(edge2, pEdge, result);
	if (dotProduct(planeNorm, result) < 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	minus(point, vertices[2], pEdge);
	crossProduct(edge3, pEdge, result);
	if (dotProduct(planeNorm, result) < 0)
	{
		//no intersection
		return GZ_FAILURE;
	}

	return GZ_SUCCESS; //ray hits triangle :)
}