class Ray {
public:
	GzCoord origin;
	GzCoord direction;
	
	int PointAtTValue(float t, GzCoord coord);

	void Ray::add(GzCoord s1, GzCoord s2, GzCoord res)
	{
		for (int i = 0; i < 3; i++)
		{
			res[i] = s1[i] + s2[i];

		}

	}

	void Ray::scalarProduct(GzCoord s, float n, GzCoord res)
	{
		for (int i = 0; i < 3; i++)
		{
			res[i] = s[i] * n;
		}

	}
};


class Triangle {
public:
	GzCoord vertices[3]; //use interpolated z-values

	int RayIntersection(Ray ray);

	float Triangle::dotProduct(GzCoord s1, GzCoord s2)
	{
		float res = 0;
		for (int i = 0; i < 3; i++)
		{
			res += s1[i] * s2[i];
		}
		return res;
	}

	void Triangle::crossProduct(GzCoord s1, GzCoord s2, GzCoord res)
	{
		res[0] = s1[1] * s2[2] - s1[2] * s2[1];
		res[1] = s1[2] * s2[0] - s1[0] * s2[2];
		res[2] = s1[0] * s2[1] - s1[1] * s2[0];

	}
	void Triangle::normalize(GzCoord s, GzCoord res)
	{
		float value = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
		for (int i = 0; i < 3; i++)
		{
			res[i] = s[i] / value;

		}

	}
	void Triangle::minus(GzCoord s1, GzCoord s2, GzCoord res)
	{
		for (int i = 0; i < 3; i++)
		{
			res[i] = s1[i] - s2[i];

		}

	}
	void Triangle::add(GzCoord s1, GzCoord s2, GzCoord res)
	{
		for (int i = 0; i < 3; i++)
		{
			res[i] = s1[i] + s2[i];

		}

	}
	void Triangle::scalarProduct(GzCoord s, float n, GzCoord res)
	{
		for (int i = 0; i < 3; i++)
		{
			res[i] = s[i] * n;

		}

	}

};
