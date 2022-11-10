Student Name: Wanyu Zhang
Student ID: 6773445719
Email: wanyuzha@usc.edu
Visual Studio Version: Visual Studio 2019


At HW5, three functions are added to GzRender class in order that different modes could use same functions and avoid code being too long.

User help function is added additionally to GzRender without changing its original API.
// User add function compute color at a normal
void GzComputeColor(GzCoord norm, GzColor color);
//To compute color avoiding a bunch of params on param list


// User add function compute (R,G,B) plane coeff or (N1,N2,N3) coeff
void GzComputePlane(GzCoord* vertices, GzCoord* coeff);
void GzComputePlane(GzCoord* vertices, GzCoord* coeff, float (*plane)[4]);
float planes[3][4];
//To interpolate on plane used both for GZ_COLOR GZ_NORMALS
//To interpolate on uv_plane