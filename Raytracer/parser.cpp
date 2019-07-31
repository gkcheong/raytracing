//File parsing example
#include <omp.h>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>

#include "image.h"
#include "Vect.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;


//function prototype
bool InTriangle(Vect<float> point, Vect<float> v1, Vect<float> v2, Vect<float> v3);
bool helperIT(Vect<float> pointA, Vect<float> pointB, Vect<float> a, Vect<float> b);

class Camera {
public:
	Camera();
	Vect<float> pos;
	Vect<float> dir;
	float d;
	Vect<float> up;
	Vect<float> right;
	double hAngle;
	double wAngle;
	int width;
	int height;
	string output; //consider removing this

};

class Ray {
public:
	//Ray = p + td
	void normalize();
	Vect<float> d;
	Vect<float> p;
	float t;
};

class Intersection {
public:
	//only keep one intesection. Decided by the smallest non negative t;
	Intersection& operator=(const Intersection& rhs);
	bool hit; //if hit == true
	float t;//amount of travel until hit. p + t*d
	Vect<float> p; //intersection point
	Vect<float>  n; //normal of the intersection point
	int objIndex; //object index
	Intersection();

};

Intersection::Intersection() {
	hit = false;
	t = 100000;
	objIndex = -1;
}

Intersection& Intersection::operator=(const Intersection& rhs) {
	if (this != &rhs) {
		this->hit = rhs.hit;
		this->t = rhs.t;
		this->p = rhs.p;
		this->n = rhs.n;
		this->objIndex = rhs.objIndex;
	}

	return *this;
}

class Polygon {
public:
	virtual bool GetIntersection(Intersection& a, Ray& ray) = 0;
	virtual float GetMin(int axis);
	virtual float GetMax(int axis);
	Vect<float> c;
	int myindex;
};

float Polygon::GetMin(int axis) {
	return -1;
}

float Polygon::GetMax(int axis) {
	return -1;
}

typedef vector<Polygon*> Scene;

class Sphere : public Polygon{
	// a sphere has a center and a radius r
public:
	//todo write constructor and destructor
	float r;
	bool GetIntersection(Intersection& a, Ray& ray);
	float GetMin(int axis);
	float GetMax(int axis);
};

float Sphere::GetMin(int axis) {
	switch (axis) {
	case 0:
		return c[0] - r;
		break;
	case 1:
		return c[1] - r;
		break;
	case 2:
		return c[2] - r;
		break;
	default:
		return -1;
	}
}

float Sphere::GetMax(int axis) {
	switch (axis) {
	case 0:
		return c[0] + r;
		break;
	case 1:
		return c[1] + r;
		break;
	case 2:
		return c[2] + r;
		break;
	default:
		return -1;

	}
}
bool Sphere::GetIntersection(Intersection& xsec, Ray& ray) {
	Vect<float> shift = ray.p - c;
	float dis = dot(ray.d, shift);
	dis *= dis;
	dis -= dot(ray.d, ray.d) * (dot(shift, shift) - r * r);

	if (dis < -0.01) {
		//no solution
		return false; //no intersection, we move onto the next object
	}
	else if (dis < 0.01) { //one solution
		float tempF;

		tempF = -dot(ray.d, shift) + sqrt(dis);
		if (tempF < xsec.t && tempF > 0.01) {
			//this is a better solution. take this
			xsec.hit = true;
			xsec.t = tempF;
			xsec.p = ray.p + xsec.t * ray.d;
			xsec.n = xsec.p - c;
			xsec.n.normalize();
			xsec.objIndex = myindex;
			//add this after
			//xsec.objIndex = i;
			//cout << "what are you "<<i << endl;
			return true;
		}
		else {
			return false;
		}
	}
	else {
		float tempF;
		bool better_found = false;
		tempF = -dot(ray.d, shift) + sqrt(dis);
		if (tempF < xsec.t && tempF > 0.01) {
			xsec.hit = true;
			//this is a better solution. take this
			xsec.t = tempF;
			xsec.p = ray.p + xsec.t * ray.d;
			xsec.n = xsec.p - c;
			xsec.n.normalize();
			xsec.objIndex = myindex;
			//xsec.objIndex = i;
			better_found = true;
		}

		tempF = -dot(ray.d, shift) - sqrt(dis);
		if (tempF < xsec.t && tempF > 0.01) {
			//check second point
			//this is a better solution. take this
			xsec.hit = true;
			xsec.t = tempF;
			xsec.p = ray.p + xsec.t * ray.d;
			xsec.n = xsec.p - c;
			xsec.n.normalize();
			xsec.objIndex = myindex;
			//xsec.objIndex = i;
			better_found = true;
		}
		return better_found;
	}
}

class Triangle : public Polygon{
	//a triangle can be represented by 3 points and a normal/or two vectors
	//Barycentric coordinates needs two vector
	//This class contains keeps both 2 vectors and a normal
	//This class keeps one vertex only for ray plane intersection
public:
	Vect<float> p1;
	Vect<float> p2;
	Vect<float> p3;
	Vect<float> n1;
	Vect<float> n2;
	Vect<float> n3;
	bool GetIntersection(Intersection& a, Ray& ray);
	float GetMin(int axis);
	float GetMax(int axis);
};


bool Triangle::GetIntersection(Intersection& xsec, Ray& ray) {
	float tempt = 0;
	bool inside;
	//obtain ray-plane intersection
	//cout << "tempt " << -(dot(ray.p, n1) + -dot(p1, n1)) / dot(ray.d, n1) <<' '<<myindex<< endl;
	tempt = -(dot(ray.p - p1, n1) ) /
		dot(ray.d, n1);
	//cout << "ray hits " << ray.p[0] + tempt*ray.d[0] << ' ' << ray.p[1] + tempt*ray.d[1] << ' ' << ray.p[2] + tempt*ray.d[2] << endl;
	inside = InTriangle(ray.p + tempt*ray.d, p1, p2, p3);
	
	if (tempt < xsec.t && tempt > 0.00001 && inside) {
		//this is a better solution. take this
		xsec.hit = true;
		xsec.t = tempt;
		xsec.p = ray.p + xsec.t * ray.d;
		xsec.n = n1;
		xsec.n.normalize();
		xsec.objIndex = myindex;
		//xsec.objIndex = i + scene.size();
		return true;
	}
	return false;
}

float Triangle::GetMin(int axis) {
	float min;
	switch (axis) {
	case 0:
		min = p1[0];
		if (p2[0] < min) {
			min = p2[0];
		}
		if (p3[0] < min) {
			min = p3[0];
		}
		return min;
		break;
	case 1:
		min = p1[1];
		if (p2[1] < min) {
			min = p2[1];
		}
		if (p3[1] < min) {
			min = p3[1];
		}
		return min;
		break;
	case 2:
		min = p1[2];
		if (p2[2] < min) {
			min = p2[2];
		}
		if (p3[2] < min) {
			min = p3[2];
		}
		return min;
		break;
	default:
		return -1;
	}
}

float Triangle::GetMax(int axis) {
	float max;
	switch (axis) {
	case 0:
		max = p1[0];
		if (p2[0] > max) {
			max = p2[0];
		}
		if (p3[0] > max) {
			max = p3[0];
		}
		return max;
		break;
	case 1:
		max = p1[1];
		if (p2[1] > max) {
			max = p2[1];
		}
		if (p3[1] > max) {
			max = p3[1];
		}
		return max;
		break;
	case 2:
		max = p1[2];
		if (p2[2] > max) {
			max = p2[2];
		}
		if (p3[2] > max) {
			max = p3[2];
		}
		return max;
		break;
	default:
		return -1;
	}
}

class nTriangle : public Triangle {
public:
	bool GetIntersection(Intersection& a, Ray& ray);
};

bool nTriangle::GetIntersection(Intersection& xsec, Ray& ray) {
	float t1 = (p1[1] - p3[1])*ray.d[2] - ray.d[1] * (p1[2] - p3[2]); //ei-hf
	float t2 = (ray.d[0] * (p1[2] - p3[2])) - (p1[0] - p3[0]) *ray.d[2]; //gf-di
	float t3 = (p1[0] - p3[0]) * ray.d[1] - (p1[1] - p3[1])*ray.d[0]; //dh-eg
	float M = (p1[0] - p2[0]) * t1 + (p1[1] - p2[1])*t2 + (p1[2] - p2[2]) * t3;

	float t4 = (p1[0] - p2[0]) *(p1[1] - ray.p[1]) - (p1[0] - ray.p[0])*(p1[1] - p2[1]); //ak-jb
	float t5 = (p1[0] - ray.p[0])*(p1[2] - p2[2]) - (p1[0] - p2[0])*(p1[2] - ray.p[2]); //jc-al
	float t6 = (p1[1] - p2[1]) * (p1[2] - ray.p[2]) - (p1[1] - ray.p[1]) * (p1[2] - p2[2]); // bl -kc

	float t = -((p1[2] - p3[2]) * t4 + (p1[1] - p3[1]) * t5 + (p1[0] - p3[0]) * t6) / M;
	//cout << "this is t " << t <<" and my index "<<myindex<< endl;
	if (t > xsec.t || t < 0.0001) {
		return false;
	}
	float g = (ray.d[2] * t4 + ray.d[1] * t5 + ray.d[0] * t6) / M;
	//cout << "this is g " << g << endl;
	if (g < 0 || g > 1) {
		return false;
	}
	float b = ((p1[0] - ray.p[0]) * t1 + (p1[1] - ray.p[1])  * t2 + (p1[2] - ray.p[2]) * t3) / M;
	//cout << "this is b " << b << endl;

	if (b < 0 || b > 1 - g) {
		return false;
	}

	//found a better solution
	xsec.hit = true;
	xsec.t = t;
	xsec.p = (1- b - g) *p1 + b * (p2) + g * (p3);
	xsec.n = (1- b - g)* n1 + b * (n2) + g * (n3);
	//xsec.n = n1;
	xsec.n.normalize();
	xsec.objIndex = myindex;
	//xsec.objIndex = i + scene.size();
}

class Box {
public:
	float xmin;
	float xmax;
	float ymin;
	float ymax;
	float zmin;
	float zmax;
	bool alwaysHit; //for bounding box that contains leaves
	bool Intersect(Ray& ray);
	Box();
};

Box::Box() {
	xmin = 0;
	xmax = 0;
	ymin = 0;
	ymax = 0;
	zmin = 0;
	zmax = 0;
	alwaysHit = false;
}

bool Box::Intersect(Ray& ray) {
	if (alwaysHit) {
		return true;
	}
	float txmin;
	float txmax;
	float tymin;
	float tymax;
	float tzmin;
	float tzmax;

	float temp = 1.0 / ray.d[0];
	if (temp >= 0) {
		txmin =  temp * (xmin - ray.p[0]);
		txmax = temp * (xmax - ray.p[0]);
	}
	else {
		txmax = temp * (xmin - ray.p[0]);
		txmin = temp * (xmax - ray.p[0]);
	}
	temp = 1.0 / ray.d[1];
	if (temp >= 0) {
		tymin = temp * (ymin - ray.p[1]);
		tymax = temp * (ymax - ray.p[1]);
	}
	else {
		tymax = temp * (ymin - ray.p[1]);
		tymin = temp * (ymax - ray.p[1]);
	}
	
	temp = 1.0 / ray.d[2];
	if (temp >= 0) {
		tzmin = temp * (zmin - ray.p[2]);
		tzmax = temp * (zmax - ray.p[2]);
	}
	else {
		tzmax = temp * (zmin - ray.p[2]);
		tzmin = temp * (zmax - ray.p[2]);
	}
	
	if (txmin > tymax || txmin > tzmax || tymin > txmax || tymin > tzmax
		|| tzmin > txmax || tzmin > tymax) {
		return false;
	}
	else {
		return true;
	}

}

class BvhNode : public Polygon {
public:
	Polygon* left = NULL;
	Polygon* right = NULL;
	Box bound;
	bool GetIntersection(Intersection& a, Ray& ray);
	BvhNode(Scene& scene, int axis); //constructor
	~BvhNode(); //destructor
private:
	void const_helper(Scene& scene, int& axis, vector<Polygon*>& left, vector<Polygon*>& right);
	
};

BvhNode::~BvhNode() {
	if (left == right) {
		delete left;
	}
	else{
		delete left;
		delete right;
	}
}

bool BvhNode::GetIntersection(Intersection& a, Ray& ray) {
	if (bound.Intersect(ray)) {
		Intersection lefthit;
		Intersection righthit;
		//cout << "go left" << endl;
		left->GetIntersection(lefthit, ray);
		//cout << "go right" << endl;
		right->GetIntersection(righthit, ray);
		if (lefthit.hit && righthit.hit) {
			if (lefthit.t < righthit.t) {
				a = lefthit;
				//cout << "what is left " << a.objIndex << endl;
				//cout << "what is left " << a.t << endl;
			}
			else {
				a = righthit;
				//cout << "right does hit??? " << righthit.hit << endl;
				//cout << "what is right " << a.objIndex << endl;
				//cout << "what is right " << a.t << endl;
			}
			return true;
		}
		else if (lefthit.hit) {
			a = lefthit;
			return true;
		}
		else if (righthit.hit) {
			a = righthit;
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}
}

void BvhNode::const_helper(Scene& scene, int& axis, vector<Polygon*>& forLeft, vector<Polygon*>& forRight) {
	int n = scene.size();
	//the mid point is chosen by center of gravity
	float xmin = scene[0]->c[0];
	float xmax = scene[0]->c[0];
	float ymin = scene[0]->c[1];
	float ymax = scene[0]->c[1];
	float zmin = scene[0]->c[2];
	float zmax = scene[0]->c[2];
	
	bound.xmin = scene[0]->GetMin(0);
	bound.xmax = scene[0]->GetMax(0);
	bound.ymin = scene[0]->GetMin(1);
	bound.ymax = scene[0]->GetMax(1);
	bound.zmin = scene[0]->GetMin(2);
	bound.zmax = scene[0]->GetMax(2);
	for (int i = 1; i < n; i++) {
		if (scene[i]->GetMin(0) < bound.xmin) {
			bound.xmin = scene[i]->GetMin(0);
		}
		if (scene[i]->GetMax(0) > bound.xmax) {
			bound.xmax = scene[i]->GetMax(0);
		}
		if (scene[i]->GetMin(1) < bound.ymin) {
			bound.ymin = scene[i]->GetMin(1);
		}
		if (scene[i]->GetMax(1) > bound.ymax) {
			bound.ymax = scene[i]->GetMax(1);
		}
		if (scene[i]->GetMin(2) < bound.zmin) {
			bound.zmin = scene[i]->GetMin(2);
		}
		if (scene[i]->GetMax(2) > bound.zmax) {
			bound.zmax = scene[i]->GetMax(2);
		}
		if (scene[i]->c[0] < xmin) {
			xmin = scene[i]->c[0];
		}
		if (scene[i]->c[0] > xmax) {
			xmax = scene[i]->c[0];
		}
		if (scene[i]->c[1] < ymin) {
			ymin = scene[i]->c[1];
		}
		if (scene[i]->c[1] > ymax) {
			ymax = scene[i]->c[1];
		}
		if (scene[i]->c[2] < zmin) {
			zmin = scene[i]->c[2];
		}
		if (scene[i]->c[2] > zmax) {
			zmax = scene[i]->c[2];
		}
	}
	float mid;
	/*switch (axis) {
	case 0:
		mid = (bound.xmax + bound.xmin) / 2.0;
		break;
	case 1:
		mid = (bound.ymax + bound.ymin) / 2.0;
		break;
	case 2:
		mid = (bound.zmax + bound.zmin) / 2.0;
		break;
	default:
		cout << "error in switch" << endl;
		exit(1);
	}*/
	switch (axis) {
	case 0:
		mid = (xmax + xmin) / 2.0;
		break;
	case 1:
		mid = (ymax + ymin) / 2.0;
		break;
	case 2:
		mid = (zmax + zmin) / 2.0;
		break;
	default:
		cout << "error in switch" << endl;
		exit(1);
	}
	//cout << "This is my final bound size" << endl;
	//cout << "(" << bound.xmin << ',' << bound.xmax << ")" <<
	//	"(" << bound.ymin << ',' << bound.ymax << ")" <<
	//	"(" << bound.zmin << ',' << bound.zmax << ")" << endl;


	for (int i = 0; i < n; i++) {
		switch (axis) {
		case 0:
			if (scene[i]->c[0] < mid) {
				forLeft.push_back(scene[i]);
			}
			else {
				forRight.push_back(scene[i]);
			}
			break;
		case 1:
			if (scene[i]->c[1] < mid) {
				//cout << "this is called how many times?" << endl;
				//cout << "mid is ?" << mid << endl;
				//cout << "scene[i] " << scene[i]->c[1] << endl;
				forLeft.push_back(scene[i]);
			}
			else {
				forRight.push_back(scene[i]);
			}
			break;
		case 2:
			if (scene[i]->c[2] < mid) {
				//cout << "this is called how many times?" << endl;
				forLeft.push_back(scene[i]);
			}
			else {
				//cout << "Scene size and id " << scene[i]->c[2] << ' ' << scene[i]->myindex << endl;
				//cout << "right is called how many times?" << endl;
				forRight.push_back(scene[i]);
			}
			break;
		default:
			cout << "error is dividing\n";
			break;
		}
	}
}
BvhNode::BvhNode(Scene& scene, int axis) {
	int n = scene.size();
	//cout << "This is the start of constructor " << endl;
	//cout << "given size is " << n << endl;
	if (n == 1) {
		left = scene[0]; //left is a pointer to a member of scene which is a polygon ptr
		//cout << "what is left? " << scene[0]->myindex << endl;
		right = left; //this avoids check for null but is slower as left is checked twice
		//cout << "what is right? " << scene[0]->myindex << endl;
		//bound = BoundingBox(scene[0]);
		//make it always hits
		bound.alwaysHit = true;
	}
	else if (n == 2) {
		left = scene[0];
		//cout << "what is left? " << scene[0]->myindex << endl;
		right = scene[1];
		//cout << "what is right? " << scene[1]->myindex << endl;
		//bound = Combine(BoundingBox(scene[0]), BoundingBox(scene[1]);
		//make it always hits
		bound.alwaysHit = true;
	}
	else {
		vector<Polygon*> forLeft;
		vector<Polygon*> forRight;
		
		const_helper(scene, axis, forLeft, forRight);
		int count = 0;
		for (int i = 0; i < 3; i++) {
			if (forLeft.size() == 0) {
				cout << "left bounding box empty!" << endl;
				cout << "recalculating in new axis!" << endl;
				cout << "Size of right" << forRight.size() << endl;
				forLeft.clear();
				forRight.clear();
				axis = (axis + 1) % 3;
				count++;
				const_helper(scene, axis, forLeft, forRight);
			} else if (forRight.size() == 0) {
				cout << "right bounding box empty!" << endl;
				cout << "recalculating in new axis!" << endl;
				cout << "Size of left" << forLeft.size() << endl;
				forLeft.clear();
				forRight.clear();
				axis = (axis + 1) % 3; 
				count++;
				const_helper(scene, axis, forLeft, forRight);
			}
			else {
				break; //no need to check
			}
			if (count > 2) { //looped around
				cout << "Multiple objects in the same location" << endl;
				//can we divide by the center point of the object?
				exit(1);
			}
		}
		
		left = new BvhNode(forLeft, (axis + 1) % 3);
		right = new BvhNode(forRight, (axis + 1) % 3);
	}

	
}
class Material {
public:
	Vect<float> a; //ambient term
	Vect<float> d; //diffuse term
	Vect<float> s; //specular term
	float ns;
	Vect<float> t; //transmissive term
	float ior;
};

typedef vector<Material> Surface;

class Light {
public:
	virtual Vect<float> GetColor(Intersection& point, Surface& material, Ray& shadow, Ray& ret) = 0;
	virtual void GetDir(Intersection& point, Ray& ray) = 0;
};
class PLight : public Light{
public:
	Vect<float> col;
	Vect<float> pos; 
	Vect<float> GetColor(Intersection& point, Surface& material, Ray& ray, Ray& ret);
	void GetDir(Intersection& point, Ray& ray);
};

void PLight::GetDir(Intersection& point, Ray& ray) {
	ray.d = pos - point.p;
	ray.normalize();

	return;
}

Vect<float> PLight::GetColor(Intersection& point, Surface& material, Ray& shadow, Ray& ret) {
	Vect<float> color;
	color[0] = 0; color[1] = 0; color[2] = 0;

	float fallOff = 1;
	float dist1 = 1;
	float dist;
	dist = sqrt(dot((ret.p - point.p), (ret.p - point.p)));
	//cout << "This is the distance " << dist << endl;
	Vect<float> illu; //from light source to object
	illu = shadow.d;
	dist1 = shadow.t;
	//cout << "dist 1 is " << dist1 << endl;
	illu.normalize();


	//(ret)urn ray
	//from object to cam
	//dist2 = sqrt(dot((ret.p - point.p), (ret.p - point.p)));

	//(ref)lected ray
	Vect<float> Itemp;
	Itemp[0] = -illu[0]; Itemp[1] = -illu[1]; Itemp[2] = -illu[2];
	Vect<float> ref;
	ref = Itemp - 2 * dot(Itemp, point.n) * point.n;
	ref.normalize();

	//about (30 pixel)^2 is one unit of drop off
	//fallOff = 1 / (0.08*(dist1 + dist2) *(dist1 + dist2));
	//fallOff = 1 / (0.3*(dist1*dist1)); //good for test_reasonable
	fallOff = 1 / (1.2*(dist1)*(dist1)); //good for spheres1.scn

	
	//cout << " diffuse falloff is " << fallOff << endl;
	//cout << "The dist 1 is " << dist1 << endl;
	//cout << "angle is " << fmax(0, dot(illu, point.n)) << endl;
	//diffuse component
	color += fmin(1,fallOff) * fmax(0, dot(illu, point.n)) * material[point.objIndex].d * col;

	//specular component
	//fallOff = 1 / (0.2*(dist1 + dist2) *(dist1 + dist2));
	fallOff = 1 / (1*(dist1)*(dist1)); //good for spheres1.scn
	
	//cout << "falloff is " << fallOff << endl;
	color += fmin(1,fallOff) * material[point.objIndex].s * col * pow(fmax(0, dot(-ret.d, ref)), material[point.objIndex].ns);
	return color;
	//return color;
}

class DLight : public Light{
public:
	Vect<float> col;
	Vect<float> dir; //direction
	Vect<float> GetColor(Intersection& point, Surface& material, Ray& ray, Ray& ret);
	void GetDir(Intersection& point, Ray& ray);
};

void DLight::GetDir(Intersection& point, Ray& ray) {
	ray.d = -dir;
	ray.normalize();
	ray.t = 100000000; //rays comes from infinity;
	return;
}

Vect<float> DLight::GetColor(Intersection& point, Surface& material, Ray& shadow, Ray& ret) {
	Vect<float> color;
	color[0] = 0; color[1] = 0; color[2] = 0;

	float fallOff = 1;
	float dist1 = 1;
	float dist2 = 1;
	
	Vect<float> illu; //from light source to object
	illu = shadow.d;
	illu.normalize();


	//(ret)urn ray
	//automatically set by Ray ret

	//(ref)lected ray
	Vect<float> Itemp;
	Itemp[0] = -illu[0]; Itemp[1] = -illu[1]; Itemp[2] = -illu[2];
	Vect<float> ref;
	ref = Itemp - 2 * dot(Itemp, point.n) * point.n;
	ref.normalize();

	//diffuse component
	//cout << "the ouput here is " << material[point.objIndex].d[1]  << endl;
	color += fmax(0, dot(illu, point.n)) * material[point.objIndex].d * col;

	//specular component
	color += material[point.objIndex].s * col * pow(fmax(0, dot(-ret.d, ref)), material[point.objIndex].ns);
	return color;
}

class SLight : public Light{
public:
	Vect<float> col;
	Vect<float> pos;
	Vect<float> dir;
	float angle1;
	float angle2;
	Vect<float> GetColor(Intersection& point, Surface& material, Ray& ray, Ray& ret);
	void GetDir(Intersection& point, Ray& ray);
};

void SLight::GetDir(Intersection& point, Ray& ray) {
	ray.d = pos - point.p;
	ray.normalize();
	return;

}

Vect<float> SLight::GetColor(Intersection& point, Surface& material, Ray& shadow, Ray& ret) {
	Vect<float> color;
	color[0] = 0; color[1] = 0; color[2] = 0;

	float fallOff = 1;
	float dist1 = 1;
	float dist2 = 1;

	Vect<float> illu; //from light source to object
	illu = shadow.d;
	dist1 = shadow.t;
	//cout << "dist 1 is " << dist1 << endl;
	illu.normalize();

	//(ret)urn ray
	//from object to cam
	dist2 = sqrt(dot((ret.p - point.p), (ret.p - point.p)));

	//(ref)lected ray
	Vect<float> Itemp;
	Itemp[0] = -illu[0]; Itemp[1] = -illu[1]; Itemp[2] = -illu[2];
	Vect<float> ref;
	ref = Itemp - 2 * dot(Itemp, point.n) * point.n;
	ref.normalize();

	Itemp.normalize();
	dir.normalize();

	double angle = acos(dot(Itemp, dir)) / 3.1415 * 180;

	if (angle < angle1 && angle > 0) {

		//about (30 pixel)^2 is one unit of drop off
		fallOff = 1 / (0.08*(dist1 + dist2) *(dist1 + dist2));

		//diffuse component
		color += fmin(1,fallOff) * fmax(0, dot(illu, point.n)) * material[point.objIndex].d * col;

		//specular component
		fallOff = 1 / (0.2*(dist1 + dist2) *(dist1 + dist2));
		color += fmin(1,fallOff) * material[point.objIndex].s * col * pow(fmax(0, dot(ret.d, ref)), material[point.objIndex].ns);

	}
	else if (angle < angle2 && angle > 0) {

		//about (30 pixel)^2 is one unit of drop off
		fallOff = 1 / (0.08*(dist1 + dist2) *(dist1 + dist2));

		//diffuse component
		//modulate effect of light by the amount away from angle 2
		//value can take (angle1, angle2)
		color += (angle2 - angle) / (angle2 - angle1) *
			fmin(1,fallOff) * fmax(0, dot(illu, point.n)) * material[point.objIndex].d * col;

		//specular component
		fallOff = 1 / (0.2*(dist1 + dist2) *(dist1 + dist2));
		color += (angle2 - angle) / (angle2 - angle1) *
			fmin(1,fallOff) * material[point.objIndex].s * col * pow(fmax(0, dot(-ret.d, ref)), material[point.objIndex].ns);
	}
	else {
		//direction of ray larger than angle2. contribute nothing
		return color;
	}
	return color;

}

//use just a sphere scene right now

typedef vector<Light*> Lights;
typedef vector<Vect<float>> Vertex;
typedef vector<Vect<float>> Normal;


void Ray::normalize() {
	t = d.norm();
	d /= t;
}

Camera::Camera() {
	pos[0] = 0; pos[1] = 0; pos[2] = 0;
	dir[0] = 0; dir[1] = 0; dir[2] = 1;
	up[0] = 0; up[1] = 1; up[2] = 0;
	d = 1;
	hAngle = 45;
	width = 640; //640
	height = 480; //480
	output = "raytraced.bmp";
	right = cross(dir, up);
	wAngle = atan2(width * tan(hAngle / 180 * 3.1415), height) / 3.1415 * 180;
}



//function prototype
Vect<float> EvaluateRay(Ray& ray, Surface& material, 
						Vect<float>& background, Vect<float>& ambient, Lights& lights, 
						Scene& scene,BvhNode& head, int ii, int jj, int maxDepth, int depth);
Vect<float> ApplyLightingColor(Intersection& point, Ray& ray, Surface& material,
								Vect<float>& background, Vect<float>& ambient, Lights& lights,
								Scene& scene, BvhNode& head, int ii, int jj, int maxDepth, int depth);

/*Vect<float> topC;
Vect<float> botC;
Vect<float> leftC;
Vect<float> rightC;
*/
Vect<float> topLeftC;
Vect<float> topRightC;
Vect<float> botLeftC;
Vect<float> botRightC;

float focusLength;
bool CornersAreSet = false;

Ray GetRayFromFilm(Camera& cam, int x, int y) {
	if (!CornersAreSet) { //set top bot left right
		//focusLength = cam.height / 2.0 / tan(cam.hAngle / 180 * 3.1415);
		focusLength = cam.d;
		
		//transform cam.pos to new basis
		/*Vect<float> origin;
		origin[0] = dot(cam.pos, cam.right);
		origin[1] = dot(cam.pos, cam.up);
		origin[2] = dot(cam.pos, cam.dir);
		cout << "Check effect of origin " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << endl;*/

		//add an amount in the new basis from camera
		//cout << "cam dir " << cam.dir[0] << ' ' << cam.dir[1] << ' ' << cam.dir[2] << endl;
		
		topLeftC = cam.pos + focusLength*cam.dir - focusLength*tan(cam.wAngle / 180 *3.1415) *cam.right
					+ focusLength*tan(cam.hAngle / 180 *3.1415) * cam.up;
		topRightC = cam.pos + focusLength*cam.dir + focusLength*tan(cam.wAngle/180 *3.1415)*cam.right
					+ focusLength*tan(cam.hAngle/180*3.1415)*cam.up;
		botLeftC = cam.pos + focusLength*cam.dir - focusLength*tan(cam.wAngle/180*3.1415)*cam.right
					- focusLength*tan(cam.hAngle/180*3.1415)*cam.up;
		botRightC = cam.pos + focusLength*cam.dir + focusLength*tan(cam.wAngle/180*3.1415)*cam.right
						- focusLength*tan(cam.hAngle / 180 * 3.1415)*cam.up;
		
		/*leftC = cam.pos + focusLength*cam.dir - focusLength*tan(cam.wAngle / 180 * 3.1415)*cam.right;
		rightC = cam.pos + focusLength*cam.dir + focusLength*tan(cam.wAngle / 180 * 3.1415)*cam.right;
		topC = cam.pos + focusLength*cam.dir + focusLength*tan(cam.hAngle / 180 * 3.1415)*cam.up;
		botC = cam.pos + focusLength*cam.dir - focusLength*tan(cam.hAngle / 180 * 3.1415)*cam.up;
		
		Vect<float> temp = botC;
		temp.normalize();
		cout << " bottom vec " << temp[0] << ' ' << temp[1] << ' ' << temp[2] << endl;
		*/

		//transform the corners back in standard basis
		//to do make a matrix class
		/*Vect<float> basis1;
		basis1[0] = cam.right[0]; basis1[1] = cam.up[0]; basis1[2] = cam.dir[0];
		Vect<float> basis2;
		basis2[0] = cam.right[1]; basis2[1] = cam.up[1]; basis2[2] = cam.dir[1];
		Vect<float> basis3;
		basis3[0] = cam.right[2]; basis3[1] = cam.up[2]; basis3[2] = cam.dir[2];

		Vect<float> temp;
		temp[0] = dot(leftC, basis1);
		temp[1] = dot(leftC, basis2);
		temp[2] = dot(leftC, basis3);
		leftC = temp;

		temp[0] = dot(rightC, basis1);
		temp[1] = dot(rightC, basis2);
		temp[2] = dot(rightC, basis3);
		rightC = temp;

		temp[0] = dot(topC, basis1);
		temp[1] = dot(topC, basis2);
		temp[2] = dot(topC, basis3);
		topC = temp;

		temp[0] = dot(botC, basis1);
		temp[1] = dot(botC, basis2);
		temp[2] = dot(botC, basis3);
		botC = temp;*/

		//debug code for GetRayFromFilm
		/*cout << "Check left corners " << leftC[0] << ' ' << leftC[1] << ' ' << leftC[2] << endl;
		cout << "Check right corners " << rightC[0] << ' ' << rightC[1] << ' ' << rightC[2] << endl;
		cout << "Check top corners " << topC[0] << ' ' << topC[1] << ' ' << topC[2] << endl;
		cout << "Check bot corners " << botC[0] << ' ' << botC[1] << ' ' << botC[2] << endl;
		*/

		CornersAreSet = true;
	}
	//interpolate according to stb_image_write convention (top left 0,0)
	Ray ray;
	ray.p = cam.pos;
	
	float temp1;
	float temp2;
	temp1 = (topLeftC[0] + (topRightC[0] - topLeftC[0]) * (x + 0.5) / (cam.width));
	temp2 = (botLeftC[0] + (botRightC[0] - botLeftC[0]) * (x + 0.5) / (cam.width));
	ray.d[0] = (temp1 + (temp2 - temp1) * (y+0.5)/(cam.height)) - ray.p[0];
	temp1 = (topLeftC[1] + (topRightC[1] - topLeftC[1]) * (x + 0.5) / (cam.width));
	temp2 = (botLeftC[1] + (botRightC[1] - botLeftC[1]) * (x + 0.5) / (cam.width));
	ray.d[1] = (temp1 + (temp2 - temp1) * (y+0.5)/(cam.height)) - ray.p[1];
	temp1 = (topLeftC[2] + (topRightC[2] - topLeftC[2]) * (x + 0.5) / (cam.width));
	temp2 = (botLeftC[2] + (botRightC[2] - botLeftC[2]) * (x + 0.5) / (cam.width));
	ray.d[2] = (temp1 + (temp2 - temp1) * (y+0.5)/(cam.height)) - ray.p[2];
	

	ray.normalize();
	//cout << "This is the ray " << ray.d[0] << ' ' << ray.d[1] << ' ' << ray.d[2] << endl;
	
	return ray;

}

bool helperIT(Vect<float> pointA, Vect<float> pointB, Vect<float> a, Vect<float> b) {
	Vect<float> temp1 = cross(b - a, pointA - a);
	Vect<float> temp2 = cross(b - a, pointB - a);
	return (dot(temp1, temp2) >= 0);
}

bool InTriangle(Vect<float> point, Vect<float> v1, Vect<float> v2, Vect<float> v3) {
	return helperIT(point, v1, v2, v3) &&
		helperIT(point, v2, v1, v3) &&
		helperIT(point, v3, v1, v2);
}

void FindIntersection(Ray& ray, Scene& scene, Intersection& xsec) {

	xsec.hit = false;
	bool insec_found;
	for (int i = 0; i < scene.size(); i++) {
		insec_found = scene[i]->GetIntersection(xsec, ray);
		if (insec_found) {
			xsec.objIndex = i;
		}

	}
}

void newFindIntersection(Ray& ray, BvhNode& head, Intersection& xsec) {
	//cout << "do we ever get here?" << endl;
	xsec.hit = false;
	head.GetIntersection(xsec, ray);
}
Ray Refract(Ray& illu, Intersection& point, Surface& material) {
	//Ray illu is towards object
	Ray temp;
	temp.p = point.p;
	float thetaI;
	float thetaR;
	thetaI = acos(dot(-illu.d, point.n)); //in radians;
	if (dot(illu.d, point.n) < 0) {//if from air i.e entering
		thetaR = asin(1 * sin(thetaI) / material[point.objIndex].ior);
		temp.d = (1 / material[point.objIndex].ior * cos(thetaI) - cos(thetaR)) * point.n
			- 1 / material[point.objIndex].ior * -illu.d;
	}
	else {
		thetaI = 3.1415 - thetaI;
		thetaR = asin(material[point.objIndex].ior * sin(thetaI) / 1);
		temp.d = (material[point.objIndex].ior / 1 * cos(thetaI) - cos(thetaR)) * -point.n
			- material[point.objIndex].ior / 1 * -illu.d;
	}
	temp.d.normalize();
	return temp;
}


Ray Reflect(Ray& illu, Intersection& point ) {
	//i.e must point towards object
	//Ray illu is cam to object i.e point.p - cam.p
	Ray temp;
	temp.p = point.p;
	temp.d = illu.d - 2 * dot(illu.d, point.n) * point.n;
	temp.d.normalize();
	return temp;

	/*Vect<float> Itemp;
	Itemp[0] = -illu[0]; Itemp[1] = -illu[1]; Itemp[2] = -illu[2];
	Vect<float> ref;
	ref = Itemp - 2 * dot(Itemp, point.n) * point.n;
	ref.normalize();*/
}

Vect<float> ApplyLightingColor(Intersection& point, Ray& ret, Surface& material
							, Vect<float>& background, Vect<float>& ambient, Lights& lights
							, Scene& scene,BvhNode& head, int ii, int jj, int maxDepth, int depth) {
	Vect<float> color;
	
	color[0] = 0; color[1] = 0; color[2] = 0;
	
	for (int i = 0; i < lights.size(); ++i) {

		Ray shadow;
		shadow.p = point.p;
		lights[i]->GetDir(point, shadow);
		Intersection shadow_xsec;
		
		//FindIntersection(shadow, scene, shadow_xsec);
		newFindIntersection(shadow, head, shadow_xsec);
		
		if (shadow_xsec.hit && shadow_xsec.t < shadow.t
			&& shadow_xsec.t > 0.00001) {
			continue;
		}
		
		Vect<float> temp;
		temp = lights[i]->GetColor(point, material, shadow, ret);
		
		color += temp;
	}

	if (depth < maxDepth) {
		
		Ray reflect = Reflect(ret, point);
		Vect<float> temp;
		
		temp =  material[point.objIndex].s * EvaluateRay(reflect, material, background,
										ambient, lights, scene,head, ii, jj, maxDepth, depth + 1);
		
		color +=  temp;

		Ray refract = Refract(ret, point, material);
		temp = material[point.objIndex].t * EvaluateRay(refract, material, background,
										ambient, lights, scene,head, ii, jj, maxDepth, depth + 1);
		
		color += temp;
	}
	
	
	color += material[point.objIndex].a * ambient;
	
	return color;

}

Vect<float> EvaluateRay(Ray& ray, Surface& material, Vect<float>& background, 
						Vect<float>& ambient, Lights& lights, Scene& scene, BvhNode& head,
						int ii, int jj, int maxDepth, int depth) {

	// ii and jj for debugging purposes will remove in next iteration

	Intersection xsec;
	//reset hit record at every level
	xsec.t = 100000; //do this by setting a large distance
	
	//FindIntersection(ray, scene, xsec);
	
	newFindIntersection(ray, head, xsec);

	
	
	
	Vect<float> temp;
	if (xsec.hit) {
		return ApplyLightingColor(xsec, ray, material, background, ambient, lights,
								scene,head, ii, jj, maxDepth, depth);
	}
	else {
		temp[0] = background[0]; temp[1] = background[1]; temp[2] = background[2];
	}

	return temp;
}


int main(int argc, char* argv[]){
  string line;

  if (argc > 2) {
	  cout << "Too many command line argument" << endl;
	  exit(1);
  }

  // open the file containing the scene description
  ifstream input(argv[1]);

  // check for errors in opening the file
  if(input.fail()){
    cout << "Can't open file '" << argv[1] << "'" << endl;
    return 0;
  }
  
  // determine the file size (this is optional -- feel free to delete the 6 lines below)
  streampos begin,end;
  begin = input.tellg();
  input.seekg(0, ios::end);
  end = input.tellg();
  cout << "File '" << argv[1] << "' is: " << (end-begin) << " bytes long.\n\n";
  input.seekg(0, ios::beg);

  
  //Loop through reading each line
  bool max_vertices_set = false;
  int max_vertices = 0;
  bool max_normals_set = false;
  int max_normals = 0;
  Camera cam;
  string command;
  Scene scene;
  int myIndex = 0;
  Surface material;

  //setting up containers for triangles
  Vertex vertices;
  Normal normals;

  //setting defaults. TO do: make a default constructor
  Material currentMat;
  currentMat.a[0] = 0; currentMat.a[1] = 0; currentMat.a[2] = 0;
  currentMat.d[0] = 1; currentMat.d[1] = 1; currentMat.d[2] = 1;
  currentMat.s[0] = 0; currentMat.s[1] = 0; currentMat.s[2] = 0;
  currentMat.ns = 5;
  currentMat.t[0] = 0; currentMat.t[1] = 0; currentMat.t[2] = 0;
  currentMat.ior = 1;

  Lights lights;

  Vect<float> background;
  Vect<float> ambient;
  ambient[0] = 0; ambient[1] = 0; ambient[2] = 0;
  background[0] = 0; background[1] = 0; background[2] = 0;

  int maxDepth = 5;

  while(input >> command) { //Read first word in the line (i.e., the command type)
    
    if (command[0] == '#'){
      getline(input, line); //skip rest of line
      cout << "Skipping comment: " << command  << line <<  endl;
      continue;
    }
    
    
    //start of camera block 
    if (command == "camera") { 

      input >> cam.pos[0] >> cam.pos[1] >> cam.pos[2] 
			>> cam.dir[0] >> cam.dir[1] >> cam.dir[2] 
			>> cam.up[0] >> cam.up[1] >> cam.up[2] 
			>> cam.hAngle;
	  
	  cam.d = cam.dir.norm();
	  cam.dir.normalize();
	  cam.up.normalize(); //can normalize dont care

      printf("Camera set to position (%f,%f,%f) looking from (%f,%f,%f) and  \
        the up direction is (%f,%f,%f) with a half angle of %f\n", cam.pos[0], cam.pos[1],
        cam.pos[2], cam.dir[0], cam.dir[1], cam.dir[2], cam.up[0], cam.up[1], cam.up[2], cam.hAngle);

	  cam.right = cross(cam.dir, cam.up); //already normalize

    } else if (command == "film_resolution") {

      input >> cam.width >> cam.height;
      printf("The desired resolution of the output image is %i x %i\n", cam.width
        , cam.height);

    } else if (command == "output_image"){ //If the command is an output_image command

       input >> cam.output;
       printf("Render to file named: %s\n", cam.output.c_str());

    } //start of geometry block
    else if (command == "max_vertices") {

      max_vertices_set = true;
      input >> max_vertices;
      printf("Max vertices set to %i\n", max_vertices);
	  vertices.reserve(max_vertices);

    } else if (command == "max_normals") {

      input >> max_normals;
      max_normals_set = true;
      printf("Max normals set to %i\n", max_vertices);
	  normals.reserve(max_normals);

    } else if (command == "vertex") {

      if(!max_vertices_set) {
        cout << "Error: max vertices not set\n"<<endl;
        exit(1);
      }
	  Vect<float> temp;
      input >> temp[0] >> temp[1] >> temp[2];
	  if (vertices.size() >= max_vertices) {
		  cout << " Error: too many vertices\n" << endl;
	  }
	  vertices.push_back(temp);
      //printf("Instantiated a vertex at (%f,%f,%f)\n", temp[0], temp[1], temp[2]);


    } else if (command == "normal") {

      if(!max_normals_set) {
        cout << "Error: max normals no set\n"<<endl;
		exit(1);
      }
	  Vect<float> temp;
      input >> temp[0] >> temp[1] >> temp[2];
	  if (normals.size() >= max_normals) {
		  cout << " Error: too many normals \n" << endl;
		  exit(1);
	  }
	  normals.push_back(temp);
      printf("Instantiated a normal vector with (%f,%f,%f)\n", temp[0], temp[1], temp[2]);

    } else if (command == "triangle") {

      unsigned int v1, v2, v3;
      input >> v1 >> v2 >> v3;
	  if (v1 > vertices.capacity() - 1) {
		  cout << "Specified vertex " << v1 
			  << " is larger than vertices capacity of " << vertices.capacity() << endl;
		  exit(1);
	  }
	  if (v2 > vertices.capacity() - 1) {
		  cout << "Specified vertex " << v2
			  << " is larger than vertices capacity of " << vertices.capacity() << endl;
		  exit(1);
	  }
	  if (v3 > vertices.capacity() - 1) {
		  cout << "Specified vertex " << v3
			  << " is larger than vertices capacity of " << vertices.capacity() << endl;
		  exit(1);
	  }
	  Triangle* temp_ptr = new Triangle;
	  temp_ptr->p1 = vertices[v1];
	  temp_ptr->p2 = vertices[v2];
	  temp_ptr->p3 = vertices[v3];
	  Vect<float> vec1 = vertices[v2] - vertices[v1];
	  float temp1 = vec1.norm();
	  vec1.normalize();
	  Vect<float> vec2 = vertices[v3] - vertices[v1];
	  float temp2 = vec2.norm();
	  vec2.normalize();
	  Vect<float> vec3 = vertices[v3] - vertices[v2];
	  vec3.normalize();

	  temp_ptr->c = 1/3.0*(vertices[v1] + vertices[v2] + vertices[v3]); 
	  //quick to find a consistent `center'
	  //cout << " my triangle center is " << temp_ptr->c[0] << ' ' << temp_ptr->c[1] << ' ' << temp_ptr->c[2] << endl;
	  temp_ptr->myindex = myIndex;
	  myIndex++;
	  Vect<float> cam2tri = temp_ptr->c - cam.pos;
	  temp_ptr->n1 = cross(vec1, vec2);
	  temp_ptr->n1.normalize();
	  //cout << "temp_ptr->n1 " << temp_ptr->n1[0] << ' ' << temp_ptr->n1[1] << ' ' << temp_ptr->n1[2] << endl;
	  if (dot(cam2tri, temp_ptr->n1) > 0) {
		  temp_ptr->n1 = cross(vec2, vec1);
		  temp_ptr->n1.normalize();
	  }

	  temp_ptr->n2 = cross(vec3, -vec1);
	  temp_ptr->n2.normalize();
	  if (dot(cam2tri, temp_ptr->n2) > 0) {
		  temp_ptr->n2 = cross(vec3, -vec1);
		  temp_ptr->n2.normalize();
	  }

	  temp_ptr->n3 = cross(-vec2, -vec3);
	  temp_ptr->n3.normalize();
	  if (dot(cam2tri, temp_ptr->n3) > 0) {
		  temp_ptr->n3 = cross(-vec3, -vec2);
		  temp_ptr->n3.normalize();
	  }

	  scene.push_back(temp_ptr);
	  material.push_back(currentMat);
	  //cout << "temp_ptr->n2 " << temp_ptr->n2[0] << ' ' << temp_ptr->n2[1] << ' ' << temp_ptr->n2[2] << endl;
	  //cout << "temp_ptr->n3 " << temp_ptr->n3[0] << ' ' << temp_ptr->n3[1] << ' ' << temp_ptr->n3[2] << endl;
	  //cout << "my id " << temp_ptr->myindex << endl;
	  //cout << "This is n1 " << temp_ptr->n1[0] << ' ' << temp_ptr->n1[1] << ' ' << temp_ptr->n1[2] << endl;
      //printf("Instantiated a triangle with vertices labeled %i, %i, %i\n",
      //        v1, v2, v3);

    } else if (command == "normal_triangle") {

      unsigned int v1, v2, v3, n1, n2, n3;

      input >> v1 >> v2 >> v3 >> n1 >> n2 >> n3;
	  if (v1 > vertices.capacity() - 1) {
		  cout << "Specified vertex " << v1
			  << " is larger than vertices capacity of " << vertices.capacity() << endl;
		  exit(1);
	  }
	  if (v2 > vertices.capacity() - 1) {
		  cout << "Specified vertex " << v2
			  << " is larger than vertices capacity of " << vertices.capacity() << endl;
		  exit(1);
	  }
	  if (v3 > vertices.capacity() - 1) {
		  cout << "Specified vertex " << v3
			  << " is larger than vertices capacity of " << vertices.capacity() << endl;
		  exit(1);
	  }
	  if (n1 > normals.capacity() - 1) {
		  cout << "Specified normal " << n1
			  << " is larger than vertices capacity of " << normals.capacity() << endl;
		  exit(1);
	  }
	  if (n2 > normals.capacity() - 1) {
		  cout << "Specified normal " << n2
			  << " is larger than vertices capacity of " << normals.capacity() << endl;
		  exit(1);
	  }
	  if (n3 > normals.capacity() - 1) {
		  cout << "Specified normal " << n3
			  << " is larger than vertices capacity of " << normals.capacity() << endl;
		  exit(1);
	  }
	  nTriangle* temp_ptr = new nTriangle;
	  temp_ptr->p1 = vertices[v1];
	  temp_ptr->p2 = vertices[v2];
	  temp_ptr->p3 = vertices[v3];
	  temp_ptr->n1 = normals[n1];
	  temp_ptr->n1.normalize();
	  temp_ptr->n2 = normals[n2];
	  temp_ptr->n2.normalize();
	  temp_ptr->n3 = normals[n3];
	  temp_ptr->n3.normalize();

	  temp_ptr->c = 1 / 3.0*(vertices[v1] + vertices[v2] + vertices[v3]);
	  //quick hack to find a consistent `center'
	  //cout << " my triangle center is " << temp_ptr->c[0] << ' ' << temp_ptr->c[1] << ' ' << temp_ptr->c[2] << endl;
	  temp_ptr->myindex = myIndex;
	  myIndex++;

	  scene.push_back(temp_ptr);
	  material.push_back(currentMat);
      //printf("Instantiated a triangle with vertices %i, %i, %i, and normals \
       // %i, %i, %i\n", v1, v2, v3, n1, n2, n3);

    }else if (command == "sphere"){ //If the command is a sphere command

		Sphere* temp_ptr = new Sphere;
       input >> temp_ptr->c[0] >> temp_ptr->c[1] >> temp_ptr->c[2] >> temp_ptr->r;
	   temp_ptr->myindex = myIndex;
	   myIndex++;
       printf("Sphere as position (%f,%f,%f) with radius %f\n",
		   temp_ptr->c[0],temp_ptr->c[1],temp_ptr->c[2],temp_ptr->r);
	   scene.push_back(temp_ptr);
	   material.push_back(currentMat);
	   

    } else if (command == "background"){ //If the command is a background command

       input >> background[0] >> background[1] >> background[2];
       printf("Background color of (%f,%f,%f)\n",background[0],background[1],background[2]);

    } //start of material parameters 
    else if (command == "material") {

      input >> currentMat.a[0] >> currentMat.a[1] >> currentMat.a[2] 
			>> currentMat.d[0] >> currentMat.d[1] >> currentMat.d[2] 
            >> currentMat.s[0] >> currentMat.s[1] >> currentMat.s[2] 
			>> currentMat.ns >> currentMat.t[0] >> currentMat.t[1] >> currentMat.t[2] >> currentMat.ior;

      printf("material property is set as:\n");
      printf("Ambient : (%f, %f, %f)\n", currentMat.a[0], currentMat.a[1] ,currentMat.a[2]);
      printf("Diffuse : (%f, %f, %f)\n", currentMat.d[0], currentMat.d[1], currentMat.d[2]);
      printf("Specular: (%f, %f, %f)\n", currentMat.s[0], currentMat.s[1], currentMat.s[2]);
      printf("ns      : %f\n", currentMat.ns);
      printf("Transmit: (%f, %f, %f)\n", currentMat.t[0], currentMat.t[1], currentMat.t[2]);
      printf("IOR     : %f\n", currentMat.ior);

	 

    } //start of lightning  parameters
    else if (command == "directional_light") {
      
		DLight* temp_ptr = new DLight;
      input >> temp_ptr->col[0] >> temp_ptr->col[1] >> temp_ptr->col[2] 
			>> temp_ptr->dir[0] >> temp_ptr->dir[1] >> temp_ptr->dir[2];

      printf("Directional light with color (%f,%f,%f) from (%f,%f,%f)\n",
            temp_ptr->col[0], temp_ptr->col[1], temp_ptr->col[2], 
			temp_ptr->dir[0], temp_ptr->dir[1], temp_ptr->dir[2]);
	  //cout << "my direction is " << temp_ptr->dir[0] << endl;
	  lights.push_back(temp_ptr);

    } else if (command == "point_light") {

		PLight* temp_ptr = new PLight;
      input >> temp_ptr->col[0] >> temp_ptr->col[1] >> temp_ptr->col[2]
			>> temp_ptr->pos[0] >> temp_ptr->pos[1] >> temp_ptr->pos[2];

      printf("Point light with color (%f,%f,%f) from (%f,%f,%f)\n",
              temp_ptr->col[0], temp_ptr->col[1], temp_ptr->col[2],
			  temp_ptr->pos[0], temp_ptr->pos[1], temp_ptr->pos[2]);
	  lights.push_back(temp_ptr);

    } else if (command == "spot_light") {

	  SLight* temp_ptr = new SLight;

	  input >> temp_ptr->col[0] >> temp_ptr->col[1] >> temp_ptr->col[2]
		  >> temp_ptr->pos[0] >> temp_ptr->pos[1] >> temp_ptr->pos[2]
		  >> temp_ptr->dir[0] >> temp_ptr->dir[1] >> temp_ptr->dir[2]
		  >> temp_ptr->angle1 >> temp_ptr->angle2;

      printf("Spot light with color (%f,%f,%f) from (%f,%f,%f)\n \
              and shines toward (%f,%f,%f) with angle (%f,%f)\n",
              temp_ptr->col[0], temp_ptr->col[1], temp_ptr->col[2], 
			  temp_ptr->pos[0], temp_ptr->pos[1], temp_ptr->pos[2], 
			  temp_ptr->dir[0], temp_ptr->dir[1], temp_ptr->dir[2], 
			  temp_ptr->angle1, temp_ptr->angle2);
	  lights.push_back(temp_ptr);

    } else if (command == "ambient_light") {
      
      
      input >> ambient[0] >> ambient[1] >> ambient[2];
      printf("Ambient light with color (%f,%f,%f)\n",ambient[0],ambient[1],ambient[2]);

    } //start of Miscellaneous
    else if (command == "max_depth") {
      
      input >> maxDepth;
      printf("Max recursion depth %i\n", maxDepth);

    } else {

      getline(input, line); //skip rest of line
      cout << "WARNING. Do not know command: " << command << endl;

    }
  }
  //input completed
  //Deducing the rest of camera values;
  //cout << cam.width * tan(cam.hAngle / 180 * 3.1415) << endl;
  //cout << atan2(cam.width * tan(cam.hAngle / 180 * 3.1415), cam.height) << endl;
  //first solve for y
  float heightOfScreen = cam.d * tan(cam.hAngle / 180 * 3.1415);

  cam.wAngle = atan2(cam.width * tan(cam.hAngle / 180 * 3.1415), cam.height) / 3.1415 * 180;
  //procceding to draw on file

  //Generate image to be drawn
  Image img(cam.width, cam.height);
  //Image img(600, 800);

  //debug code for camera;
  //cout << "Checking cross product " << cam.right[0] << cam.right[1] << cam.right[2] << endl;
  //cout << "Checking wAngle " << cam.wAngle<< endl;

  //Draw on image
  float percent[11] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  int completion = 0;
  float max = cam.width * cam.height;
  clock_t timer1 = clock();
  BvhNode head(scene, 0);
  clock_t timer2 = clock();
  cout << "Time for bvh calculations " << float(timer2 - timer1) / CLOCKS_PER_SEC<<'s' << endl;
  cout << "Start Processing ---" << endl;
  
 #pragma omp parallel for schedule (dynamic)
  for (int i = 0; i < cam.width; i++) {
	  for (int j = 0; j < cam.height; j++) {

			/*race condition if openmp. 8 threads on laptop
		  if (float(i*cam.height + j) / max > percent[completion] && completion < 10) {
			  int np = omp_get_num_threads();
			  cout << "Processing: Completion rate " << percent[completion] * 100 << " % " << endl;
			  cout << "this is the amount of thread currently " << np << endl;
			  completion++;
		  }*/
		  

		  Ray ray = GetRayFromFilm(cam, i, j);

		  Vect<float> color = EvaluateRay(ray, material, background, ambient,
			  lights, scene, head, i, j, maxDepth, 0);
		  
		  img.GetPixel(i, j).SetClamp(color[0] * 255, color[1] * 255, color[2] * 255, 255);

	  }
  }
  timer2 = clock();
  cout << "Time for image calculations " << float(timer2 - timer1) / CLOCKS_PER_SEC << 's' << endl;
  for (int i = 0; i < scene.size(); ++i) {
	  delete scene[i];
  }
  for (int i = 0; i < lights.size(); ++i) {
	  delete lights[i];
  }
  //Output to file
  char* outFile = new char[cam.output.length() + 1];
  strcpy(outFile, cam.output.c_str());
  img.Write(outFile);

  delete[] outFile;

 
  
  
  return 0;

}
