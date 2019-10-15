// Raytracing Prj1 by Xipeng Wang
#include <GL\freeglut.h>
#include <omp.h>
#include "scene.h"
#include "objects.h"
#include "tinyxml.h"
#include "materials.h"
#include "lights.h"
#include "cyTimer.h"

#include <stdio.h>
#include <string.h>
#include <chrono>
#include <iostream>

// best bias for prj6
//#define bias 0.00095f
// temp bias for prj7
#define bias 0.00015f
#define longDis 10000.0f
#define e_cons 2.718281828f

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
Plane thePlane;
LightList lights;
MaterialList materials;
ItemFileList<Object> objList;
TimerStats timer;
TexturedColor background;
TexturedColor environment;
TextureList textureList;

//char prjName[] = "test7";
char prjName[] = "prj7";
char prjSource[30];
char prjRender[30];
char prjZRender[30];

int LoadScene(char const *filename);
void ShowViewport();

class CameraSpaceInfo {
public:
	Vec3f px, py, lc;

	CameraSpaceInfo() {}

	void Init(Camera c) {
		if(!c.dir.IsZero())
			c.dir.Normalize();
		if (!c.up.IsZero())
			c.up.Normalize();
		Vec3f right = c.dir.Cross(c.up);
		float ratio = (float)c.imgWidth / (float)c.imgHeight;
		float halfHeight = tanf(c.fov / 2 * Pi<float>() / 180);
		float halfWidth = halfHeight * ratio;
		lc = c.dir + c.up * halfHeight - right * halfWidth;
		px = right * halfWidth * 2 / c.imgWidth;
		py = c.up * -1 * halfHeight * 2 / c.imgHeight;
	}

	Vec3f GetPixelDir(int posX, int posY) {
		Vec3f d = lc + px * ((float)posX + 0.5f) + py * ((float)posY + 0.5f);
		return d;
	}
};
CameraSpaceInfo csInfo;

// Main Function
int main(int argc, char *argv[])
{
	omp_set_num_threads(16);

	strcpy_s(prjSource, prjName);
	strcat_s(prjSource, ".xml");

	strcpy_s(prjRender, "./Release/");
	strcat_s(prjRender, prjName);
	strcpy_s(prjZRender, prjRender);
	strcat_s(prjRender, ".png");
	strcat_s(prjZRender, "Z.png");

	LoadScene(prjSource);

	ShowViewport();
	return 0;
}

// Shadow Help function
bool ShadowRayToNode(Ray &ray, HitInfo &hInfo, Node *node) {
	Object *obj = node->GetNodeObj();

	bool hResult = false;

	if (obj) {
		hResult = obj->IntersectRay(ray, hInfo, HIT_FRONT_AND_BACK);
	}
	if (hResult) {
		if (hInfo.z < 1) {
			return true;
		}
	}

	for (int i = 0; i < node->GetNumChild(); i++) {
		Ray r = node->GetChild(i)->ToNodeCoords(ray);
		bool childResult = ShadowRayToNode(r, hInfo, node->GetChild(i));
		if (childResult) {
			if (hInfo.z < 1) {
				return true;
			}
		}
	}

	return false;
}

bool RayToNode(Ray &ray, HitInfo &hInfo, Node *node, int hitSide = HIT_FRONT) {
	Object *obj = node->GetNodeObj();

	bool hResult = false;

	if (obj) {
		hResult = obj->IntersectRay(ray, hInfo, hitSide);
	}
	if (hResult) {
		hInfo.node = node;
		node->FromNodeCoords(hInfo);
	}

	for (int i = 0; i < node->GetNumChild(); i++) {
		Ray r = node->GetChild(i)->ToNodeCoords(ray);
		bool childResult = RayToNode(r, hInfo, node->GetChild(i), hitSide);
		if (childResult) {
			hResult = childResult;
			node->FromNodeCoords(hInfo);
		}
	}

	if (hResult) {
		return true;
	}
	else {
		return false;
	}
}

void TransToNode(Ray &ray, Node* startNode, const Node* aimNode) {
	Node* currentParent = startNode;
	Node* currentChild = nullptr;
}

float IntersectBox(const Ray &ray, Box box) {
	float hit = longDis + 1;

	Vec3f dir = ray.dir;
	Vec3f p = ray.p;

	float x1 = hit, x2 = hit;
	float y1 = hit, y2 = hit;
	float z1 = hit, z2 = hit;

	if (dir.x == 0) {
		if (p.x < box.pmin.x || p.x > box.pmax.x) {
			return hit;
		}
		else {
			x1 = bias;
			x2 = longDis;
		}
	}
	else {
		x1 = (box.pmin.x - p.x) / dir.x;
		x2 = (box.pmax.x - p.x) / dir.x;
		if (x2 < x1) {
			x1 = x1 + x2;
			x2 = x1 - x2;
			x1 = x1 - x2;
		}
		if (x1 < bias) {
			if (x2 < bias) {
				return hit;
			}
			else {
				x1 = bias;
			}
		}
	}
	if (dir.y == 0) {
		if (p.y < box.pmin.y || p.y > box.pmax.y) {
			return hit;
		}
		else {
			y1 = bias;
			y2 = longDis;
		}
	}
	else {
		y1 = (box.pmin.y - p.y) / dir.y;
		y2 = (box.pmax.y - p.y) / dir.y;
		if (y2 < y1) {
			y1 = y1 + y2;
			y2 = y1 - y2;
			y1 = y1 - y2;
		}
		if (y1 < bias) {
			if (y2 < bias) {
				return hit;
			}
			else {
				y1 = bias;
			}
		}
	}
	if (dir.x == 0) {
		if (p.z < box.pmin.z || p.z > box.pmax.z) {
			return hit;
		}
		else {
			z1 = bias;
			z2 = longDis;
		}
	}
	else {
		z1 = (box.pmin.z - p.z) / dir.z;
		z2 = (box.pmax.z - p.z) / dir.z;
		if (z2 < z1) {
			z1 = z1 + z2;
			z2 = z1 - z2;
			z1 = z1 - z2;
		}
		if (z1 < bias) {
			if (z2 < bias) {
				return hit;
			}
			else {
				z1 = bias;
			}
		}
	}

	float t1 = x1, t2 = x2;
	if (t1 > y2 || t2 < y1) {
		return hit;
	}
	else {
		t1 = max(t1, y1);
		t2 = Min(t2, y2);
	}
	if (t1 > z2 || t2 < z1) {
		return hit;
	}
	else {
		t1 = max(t1, z1);
		t2 = Min(t2, z2);
	}

	hit = t1;

	return hit;
}

Color ShadePixel(Ray &ray, HitInfo &hInfo, Node *node, Vec2f relativePos) {
	Color c = Color(0, 0, 0);

	bool hResult = RayToNode(ray, hInfo, node);
	if (hResult) {
		const Node* hitNode = hInfo.node;
		c = hitNode->GetMaterial()->Shade(ray, hInfo, lights, 5);
	}
	else {
		// cylinder background mapping
		//Vec3f r = ray.dir;
		// float u = 1 / (2 * Pi<float>()) * atan2f(r.y, r.x) + 0.5f;
		// float v = 1 / (Pi<float>()) * asinf(r.z) + 0.5f;
		//Vec3f uvw = Vec3f(u, v, 0);

		// plane background mapping
		Vec3f uvw = Vec3f(relativePos, 0);
		c = background.Sample(uvw);
	}

	return c;
}

Color24 CalculatePixelHit(int posX, int posY) {
	Color c = Color();

	Vec3f rayp = camera.pos;
	Vec3f rayd = csInfo.GetPixelDir(posX, posY);

	Ray r = Ray(rayp, rayd);
	HitInfo h = HitInfo();

	bool result = RayToNode(r, h, &rootNode);

	if (result) {
		c = Color().White();
	}
	else {
		c = Color().Black();
	}

	Color24 color = Color24(c);
	return color;
}

Color24 CalculatePixelColor(int posX, int posY) {
	Color c = Color();

	Vec3f rayp = camera.pos;
	Vec3f rayd = csInfo.GetPixelDir(posX, posY);

	Ray r = Ray(rayp, rayd);
	HitInfo h = HitInfo();

	float relativeX = (float)posX / (float)camera.imgWidth;
	float relativeY = (float)posY / (float)camera.imgHeight;

	c = ShadePixel(r, h, &rootNode, Vec2f(relativeX, relativeY));

	Color24 color = Color24(c);
	return color;
}

float CalculatePixelZ(int posX, int posY) {
	float z = BIGFLOAT;

	Vec3f rayp = camera.pos;
	Vec3f rayd = csInfo.GetPixelDir(posX, posY);

	Ray r = Ray(rayp, rayd);
	HitInfo h = HitInfo();

	bool result = RayToNode(r, h, &rootNode);

	if (result) {
		z = h.z;
	}

	return z;
}

void BeginRender() {
	csInfo.Init(camera);

	int width = renderImage.GetWidth();
	int height = renderImage.GetHeight();
	Color24* img = renderImage.GetPixels();
	float* zBuffer = renderImage.GetZBuffer();

	
	timer.Start();
	// build bvh here

#pragma omp parallel for
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			Color24 c = CalculatePixelColor(i, j);
			img[j * width + i] = c;
			float z = CalculatePixelZ(i, j);
			zBuffer[j * width + i] = z;
		}
	}

	renderImage.ComputeZBufferImage();

	timer.Stop();
	auto end = timer.GetAverage();
	printf("Elapsed time: %f\n", end);
}

void StopRender() {

}

// Outside Functions
Color MtlBlinn::Shade(
	Ray const &ray,
	const HitInfo &hInfo,
	const LightList &lightsUsing,
	int bounceCount) const
{
	Color c = Color(0, 0, 0);
	if (bounceCount > 0) {

		Vec3f inDir = ray.dir * -1;
		Vec3f nDir = hInfo.N.GetNormalized();
		Vec3f reflecDir = 2 * nDir.Dot(inDir) * nDir - inDir;
		reflecDir.Normalize();

		if (hInfo.front == true) {
			if (reflection.GetColor().r > 0) {
				Color currentC = Color(0, 0, 0);

				Ray reflecRay = Ray(hInfo.p, reflecDir);
				HitInfo reflectHit = HitInfo();
				bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = reflectHit.node;
					currentC = hitNode->GetMaterial()
						->Shade(reflecRay, reflectHit, lights, bounceCount - 1);
					currentC *= reflection.GetColor();
				}
				else {
					currentC = environment.SampleEnvironment(reflecDir) * reflection.GetColor();
				}
				c += currentC;
			}
			if (refraction.GetColor().r > 0) {
				Color currentC = Color(0, 0, 0);

				float n1 = 1;
				float n2 = ior;

				Vec3f inDir = ray.dir * -1;
				inDir.Normalize();
				Vec3f nDir = hInfo.N.GetNormalized();
				
				float cosFi1 = abs(inDir.Dot(nDir));
				float sinFi1 = sqrt(1 - cosFi1 * cosFi1);
				float sinFi2 = sinFi1 * n1 / n2;
				float cosFi2 = sqrt(1 - sinFi2 * sinFi2);
				Vec3f surfaceDir = (inDir - (cosFi1 * nDir));
				surfaceDir.Normalize();
				Vec3f tl = -1 * sinFi2 * surfaceDir;
				Vec3f tn = -1 * nDir * cosFi2;
				Vec3f refracDir = tl + tn;

				Color reflecC = Color(0, 0, 0);
				Color refracC = Color(0, 0, 0);

				float r0 = pow(((n1 - n2) / (n1 + n2)), 2);
				float fresnel = r0 + (1 - r0) * pow((1-cosFi1), 5);
				float kt = refraction.GetColor().r;

				Ray reflecRay = Ray(hInfo.p, reflecDir);
				HitInfo reflectHit = HitInfo();
				bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = reflectHit.node;
					reflecC = hitNode->GetMaterial()
						->Shade(reflecRay, reflectHit, lights, bounceCount - 1);
					reflecC *= (fresnel * refraction.GetColor() + reflection.GetColor());
				}
				else {
					reflecC = environment.SampleEnvironment(reflecDir) * reflection.GetColor();
				}

				Ray refracRay = Ray(hInfo.p, refracDir);
				HitInfo refracHit = HitInfo();
				hResult = RayToNode(refracRay, refracHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = refracHit.node;
					refracC = hitNode->GetMaterial()
						->Shade(refracRay, refracHit, lights, bounceCount - 1);
					Vec3f refracRoute = refracHit.p - hInfo.p;
					float dis = refracRoute.Length();
					float rAbsorb = pow(e_cons, -absorption.r * dis);
					float gAbsorb = pow(e_cons, -absorption.g * dis);
					float bAbsorb = pow(e_cons, -absorption.b * dis);
					Color remainC = Color(rAbsorb, gAbsorb, bAbsorb);
					refracC *= ((1 - fresnel) * remainC * refraction.GetColor());
				}
				else {
					refracC = environment.SampleEnvironment(refracDir) * refraction.GetColor();
				}
				currentC += reflecC;
				currentC += refracC;

				c += currentC;
			}
			for (int i = 0; i < lightsUsing.size(); i++) {
				Color currentC = Color(0, 0, 0);
				Light* l = lightsUsing[i];
				if (l->IsAmbient()) {
					currentC += l->Illuminate(hInfo.p, hInfo.N) * diffuse.Sample(hInfo.uvw);
				}
				else {
					Color illu = l->Illuminate(hInfo.p, hInfo.N);
					Color diff = Color(0, 0, 0);
					Color spec = Color(0, 0, 0);

					Vec3f lDir = l->Direction(hInfo.p) * -1;
					Vec3f vDir = -ray.dir;
					Vec3f nDir = hInfo.N.GetNormalized();
					Vec3f hDir = (vDir + lDir);
					hDir = hDir / hDir.Length();

					float cosF = lDir.GetNormalized().Dot(nDir);
					if (cosF <= 0)
						cosF = 0;
					else {
						diff = illu * cosF;
						float specValue = hDir.Dot(nDir);
						if (specValue < 0)
							specValue = 0;
						else {
							specValue = pow(specValue, glossiness);
						}
						spec = illu * specValue;
					}

					currentC += diff * diffuse.Sample(hInfo.uvw);
					currentC += spec * specular.Sample(hInfo.uvw);
				}
				c += currentC;
			}
		}
		else {
			Color currentC = Color(0, 0, 0);

			float n1 = ior;
			float n2 = 1;

			Vec3f inDir = ray.dir * -1;
			inDir.Normalize();
			Vec3f nDir = hInfo.N.GetNormalized() * -1;

			Color reflecC = Color(0, 0, 0);
			Color refracC = Color(0, 0, 0);

			float cosFi1 = abs(inDir.Dot(nDir));
			float sinFi1 = sqrt(1 - cosFi1 * cosFi1);
			float sinFi2 = sinFi1 * n1 / n2;
			if (sinFi2 < 1) {
				float cosFi2 = sqrt(1 - sinFi2 * sinFi2);
				Vec3f surfaceDir = (inDir - (cosFi1 * nDir));
				surfaceDir.Normalize();
				Vec3f tl = -1 * sinFi2 * surfaceDir;
				Vec3f tn = -1 * nDir * cosFi2;
				Vec3f refracDir = tl + tn;

				Ray refracRay = Ray(hInfo.p, refracDir);
				HitInfo refracHit = HitInfo();
				bool hResult = RayToNode(refracRay, refracHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = refracHit.node;
					refracC = hitNode->GetMaterial()
						->Shade(refracRay, refracHit, lights, bounceCount - 1);
					currentC += refracC;
				}
				else {
					refracC = environment.SampleEnvironment(refracDir) * refraction.GetColor();
					currentC += refracC;
				}
			}
			else {
				Ray reflectRay = Ray(hInfo.p, reflecDir);
				HitInfo reflectHit = HitInfo();
				bool hResult = RayToNode(reflectRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = reflectHit.node;
					reflecC = hitNode->GetMaterial()
						->Shade(reflectRay, reflectHit, lights, bounceCount - 1);
					currentC += reflecC;
				}
				else {
					reflecC = environment.SampleEnvironment(reflecDir) * reflection.GetColor();
					currentC += reflecC;
				}
			}

			c += currentC;
		}
	}
	return c;
}

bool Sphere::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide) const
{
	if (hitSide == HIT_NONE) {
		return false;
	}
	else {
		float a = ray.dir.Dot(ray.dir);
		float b = 2 * ray.p.Dot(ray.dir);
		float c = ray.p.Dot(ray.p) - 1;
		float delta = b * b - 4 * a * c;
		if (delta < 0) {
			return false;
		}
		else {
			float k1 = (-b - sqrt(delta)) / (2 * a);
			float k2 = (-b + sqrt(delta)) / (2 * a);
			if (hitSide == HIT_FRONT) {
				if (k1 >= bias) {
					if (k1 < hInfo.z) {
						hInfo.z = k1;
						hInfo.p = ray.p + ray.dir * k1;
						hInfo.N = hInfo.p;
						float u = 1 / (2 * Pi<float>())*atan2f(hInfo.p.y, hInfo.p.x) + 0.5f;
						// cylinder mapping
						//float v = 0.5f * hInfo.p.z + 0.5f;
						// sphere mapping
						float v = 1 / Pi<float>() * asinf(hInfo.p.z) + 0.5f;
						hInfo.uvw = Vec3f(u, v, 0);
						hInfo.front = true;
						return true;
					}
					else {
						return false;
					}
					return true;
				}
				else {
					return false;
				}
			}
			else if (hitSide == HIT_BACK) {
				if (k2 >= bias) {
					if (k2 < hInfo.z) {
						hInfo.z = k2;
						hInfo.p = ray.p + ray.dir * k2;
						hInfo.N = hInfo.p;
						float u = 1 / (2 * Pi<float>())*atan2f(hInfo.p.y, hInfo.p.x) + 0.5f;
						// cylinder mapping
						//float v = 0.5f * hInfo.p.z + 0.5f;
						// sphere mapping
						float v = 1 / Pi<float>() * asinf(hInfo.p.z) + 0.5f;
						hInfo.uvw = Vec3f(u, v, 0);
						hInfo.front = false;
						return true;
					}
				}
				else {
					return false;
				}
			}
			else if (hitSide == HIT_FRONT_AND_BACK) {
				if (k1 >= bias) {
					if (k1 < hInfo.z) {
						hInfo.z = k1;
						hInfo.p = ray.p + ray.dir * k1;
						hInfo.N = hInfo.p;
						float u = 1 / (2 * Pi<float>())*atan2f(hInfo.p.y, hInfo.p.x) + 0.5f;
						// cylinder mapping
						//float v = 0.5f * hInfo.p.z + 0.5f;
						// sphere mapping
						float v = 1 / Pi<float>() * asinf(hInfo.p.z) + 0.5f;
						hInfo.uvw = Vec3f(u, v, 0);
						hInfo.front = true;
						return true;
					}
				}
				else if (k2 >= bias) {
					if (k2 < hInfo.z) {
						hInfo.z = k2;
						hInfo.p = ray.p + ray.dir * k2;
						hInfo.N = hInfo.p;
						float u = 1 / (2 * Pi<float>())*atan2f(hInfo.p.y, hInfo.p.x) + 0.5f;
						float v = 0.5f * hInfo.p.z + 0.5f;
						hInfo.uvw = Vec3f(u, v, 0);
						hInfo.front = false;
						return true;
					}
				}
				else {
					return false;
				}
			}
		}
	}
}

// Without bvh intersect 1
//bool TriObj::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide) const
//{
//	bool hit = false;
//	for (unsigned int i = 0; i < nf; i++) {
//		hit = hit | IntersectTriangle(ray, hInfo, hitSide, i);
//	}
//	return hit;
//}

// With bvh intersect
bool TriObj::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide) const
{
	bool hit = false;
	hit = TraceBVHNode(ray, hInfo, hitSide, bvh.GetRootNodeID());
	return hit;
}

bool TriObj::TraceBVHNode(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID) const {
	if (bvh.IsLeafNode(nodeID)) {
		unsigned const int* elements = bvh.GetNodeElements(nodeID);
		unsigned int elementCount = bvh.GetNodeElementCount(nodeID);

		bool hit = false;
		for (int i = 0; i < elementCount; i++) {
			unsigned int currentTriangle = elements[i];
			hit = hit | IntersectTriangle(ray, hInfo, hitSide, currentTriangle);
		}
		return hit;
	}
	else {
		bool hit = false;

		int child1 = bvh.GetFirstChildNode(nodeID);
		int child2 = bvh.GetSecondChildNode(nodeID);

		Box box1 = Box(bvh.GetNodeBounds(child1));
		float hit1 = IntersectBox(ray, box1);
		Box box2 = Box(bvh.GetNodeBounds(child2));
		float hit2 = IntersectBox(ray, box2);

		if (hit1 <= hit2) {
			if (hit1 < longDis) {
				hit = hit | TraceBVHNode(ray, hInfo, hitSide, child1);
			}
			if (hit2 < longDis) {
				if (hInfo.z > hit2) {
					hit = hit | TraceBVHNode(ray, hInfo, hitSide, child2);
				}
			}
		}
		else {
			if (hit2 < longDis) {
				hit = hit | TraceBVHNode(ray, hInfo, hitSide, child2);
			}
			if (hit1 < longDis) {
				if (hInfo.z > hit1) {
					hit = hit | TraceBVHNode(ray, hInfo, hitSide, child1);
				}
			}
		}
		return hit;
	}
}

bool TriObj::IntersectTriangle(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int faceID) const {
	if (hitSide == HIT_NONE) {
		return false;
	}
	else {
		bool hit = false;
		Vec3f p = ray.p;
		Vec3f d = ray.dir;

		Vec3f v0 = v[f[faceID].v[0]];
		Vec3f v1 = v[f[faceID].v[1]];
		Vec3f v2 = v[f[faceID].v[2]];

		Vec3f n = (v1 - v0).Cross(v2 - v0);

		float h = -n.Dot(v0);
		float t = -(p.Dot(n) + h) / (d.Dot(n));

		float NdotD = d.Dot(n);

		if (NdotD == 0) {
			return false;
		}

		if (t < bias) {
			return false;
		}

		Vec3f x = p + t * d;

		Vec2f v0d;
		Vec2f v1d;
		Vec2f v2d;
		Vec2f xd;

		//Vec3f drop = Vec3f(1, 1, 1);
		if (n.x > n.y) {
			if (n.x > n.z) {
				//drop = Vec3f(0, 1, 1);
				v0d = Vec2f(v0.y, v0.z);
				v1d = Vec2f(v1.y, v1.z);
				v2d = Vec2f(v2.y, v2.z);
				xd = Vec2f(x.y, x.z);
			}
			else {
				//drop = Vec3f(1, 1, 0);
				v0d = Vec2f(v0.x, v0.y);
				v1d = Vec2f(v1.x, v1.y);
				v2d = Vec2f(v2.x, v2.y);
				xd = Vec2f(x.x, x.y);
			}
		}
		else {
			if (n.y > n.z) {
				//drop = Vec3f(1, 0, 1);
				v0d = Vec2f(v0.x, v0.z);
				v1d = Vec2f(v1.x, v1.z);
				v2d = Vec2f(v2.x, v2.z);
				xd = Vec2f(x.x, x.z);
			}
			else {
				//drop = Vec3f(1, 1, 0);
				v0d = Vec2f(v0.x, v0.y);
				v1d = Vec2f(v1.x, v1.y);
				v2d = Vec2f(v2.x, v2.y);
				xd = Vec2f(x.x, x.y);
			}
		}

		float a0 = (v1d - xd).Cross(v2d - xd);
		float a1 = (v2d - xd).Cross(v0d - xd);
		float a2 = (v0d - xd).Cross(v1d - xd);

		if ((!(a0 < 0 && a1 < 0 && a2 < 0)) && (a0 < 0 || a1 < 0 || a2 < 0)) {
			return false;
		}

		float aAdd = a0 + a1 + a2;

		Vec3f bc = Vec3f(a0 / aAdd, a1 / aAdd, a2 / aAdd);
		Vec3f vertexN = GetNormal(faceID, bc);
		Vec3f uvw = GetTexCoord(faceID, bc);

		if (hitSide == HIT_FRONT) {
			if (NdotD > 0) {
				return false;
			}
			if (t < hInfo.z) {
				hInfo.z = t;
				hInfo.p = x;
				hInfo.N = vertexN;
				hInfo.uvw = uvw;
				hInfo.front = true;
				hit = true;
			}
		}
		else if (hitSide == HIT_BACK) {
			if (NdotD < 0) {
				return false;
			}
			if (t < hInfo.z) {
				hInfo.z = t;
				hInfo.p = x;
				hInfo.N = vertexN;
				hInfo.uvw = uvw;
				hInfo.front = false;
				hit = true;
			}
		}
		else if (hitSide == HIT_FRONT_AND_BACK) {
			if (t < hInfo.z) {
				hInfo.z = t;
				hInfo.p = x;
				hInfo.N = vertexN;
				hInfo.uvw = uvw;
				hInfo.front = NdotD < 0;
				hit = true;
			}
		}
		return hit;
	}
}

bool Plane::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide) const
{
	if (hitSide == HIT_NONE) {
		return false;
	}
	else {
		Vec3f n = Vec3f(0, 0, 1);
		float h = 0;
		Vec3f p = ray.p;
		Vec3f d = ray.dir;

		float NdotD = d.Dot(n);
		if (NdotD == 0) {
			return false;
		}

		float t = -(p.Dot(n) + h) / NdotD;
		if (t < bias) {
			return false;
		}

		Vec3f x = p + d * t;
		x *= Vec3f(1, 1, 0);
		Vec3f corner1 = Vec3f(-1, -1, 0);
		Vec3f corner2 = Vec3f(1, -1, 0);
		Vec3f corner3 = Vec3f(1, 1, 0);
		Vec3f corner4 = Vec3f(-1, 1, 0);

		bool hitOnPlane = true;
		if ((corner1 - x).Cross(corner2 - x).z < 0) {
			hitOnPlane = false;
		}
		if ((corner2 - x).Cross(corner3 - x).z < 0) {
			hitOnPlane = false;
		}
		if ((corner3 - x).Cross(corner4 - x).z < 0) {
			hitOnPlane = false;
		}
		if ((corner4 - x).Cross(corner1 - x).z < 0) {
			hitOnPlane = false;
		}

		if (!hitOnPlane) {
			return false;
		}

		if (hitSide == HIT_FRONT) {
			if (NdotD > 0) {
				return false;
			}
			else {
				if (t < hInfo.z) {
					hInfo.z = t;
					hInfo.p = x;
					hInfo.N = n;
					// up-left corner uv
					//float u = (x.x - corner4.x) / 2;
					//float v = (corner4.y - x.y) / 2;
					// bottom-left corner uv
					float u = abs(x.x - corner1.x) / 2;
					float v = abs(x.y - corner1.y) / 2;
					hInfo.uvw = Vec3f(u, v, 0);
					hInfo.front = true;
					return true;
				}
				return false;
			}
		}
		else if (hitSide == HIT_BACK) {
			if (NdotD < 0) {
				return false;
			}
			else {
				if (t < hInfo.z) {
					hInfo.z = t;
					hInfo.p = x;
					hInfo.N = n;
					// up-left corner uv
					// float u = (x.x - corner4.x) / 2;
					// float v = (corner4.y - x.y) / 2;
					// bottom-left corner uv
					float u = abs(x.x - corner1.x) / 2;
					float v = abs(x.y - corner1.y) / 2;
					hInfo.uvw = Vec3f(u, v, 0);
					hInfo.front = false;
					return true;
				}
				return false;
			}
		}
		else if (hitSide == HIT_FRONT_AND_BACK) {
			if (t < hInfo.z) {
				hInfo.z = t;
				hInfo.p = x;
				hInfo.N = n;
				// up-left corner uv
				// float u = (x.x - corner4.x) / 2;
				// float v = (corner4.y - x.y) / 2;
				// bottom-left corner uv
				float u = abs(x.x - corner1.x) / 2;
				float v = abs(x.y - corner1.y) / 2;
				hInfo.uvw = Vec3f(u, v, 0);
				hInfo.front = NdotD < 0;
				return true;
			}
			return false;
		}
	}
}

float GenLight::Shadow(Ray ray, float t_max) {
	HitInfo shadowHitInfo = HitInfo();
	float dis = 1;
	if (t_max != 1) {
		ray.Normalize();
		dis = longDis;
	}

	Ray shadowRay = Ray(ray.p, ray.dir * dis);
	bool hit = ShadowRayToNode(shadowRay, shadowHitInfo, &rootNode);

	if (hit)
		return 0;
	else
		return 1;
}