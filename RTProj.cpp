// Raytracing Prj1 by Xipeng Wang
#include <GL\freeglut.h>
#include "scene.h"
#include "objects.h"
#include "tinyxml.h"
#include "materials.h"
#include "lights.h"

#include <stdio.h>
#include <string.h>
#include <chrono>
#include <iostream>

#define bias 0.00095f
#define longDis 10000.0f
#define e_cons 2.718281828f

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
LightList lights;
MaterialList materials;

//char prjName[] = "test3";
char prjName[] = "prj4";
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

Color ShadePixel(Ray &ray, HitInfo &hInfo, Node *node) {
	Color c = Color(0, 0, 0);

	bool hResult = RayToNode(ray, hInfo, node);
	if (hResult) {
		const Node* hitNode = hInfo.node;
		c = hitNode->GetMaterial()->Shade(ray, hInfo, lights, 5);
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

	c = ShadePixel(r, h, &rootNode);

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
	auto start = std::chrono::system_clock::now();

	csInfo.Init(camera);

	int width = renderImage.GetWidth();
	int height = renderImage.GetHeight();
	Color24* img = renderImage.GetPixels();
	float* zBuffer = renderImage.GetZBuffer();

	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			Color24 c = CalculatePixelColor(i, j);
			img[j * width + i] = c;
			float z = CalculatePixelZ(i, j);
			zBuffer[j * width + i] = z;
		}
	}

	renderImage.ComputeZBufferImage();

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
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
			if (reflection.r > 0) {
				Color currentC = Color(0, 0, 0);

				Ray reflecRay = Ray(hInfo.p, reflecDir);
				HitInfo reflectHit = HitInfo();
				bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = reflectHit.node;
					currentC = hitNode->GetMaterial()
						->Shade(reflecRay, reflectHit, lights, bounceCount - 1);
					currentC *= reflection;
				}
				c += currentC;
			}
			if (refraction.r > 0) {
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
				float kt = refraction.r;

				Ray reflecRay = Ray(hInfo.p, reflecDir);
				HitInfo reflectHit = HitInfo();
				bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
				if (hResult) {
					const Node* hitNode = reflectHit.node;
					reflecC = hitNode->GetMaterial()
						->Shade(reflecRay, reflectHit, lights, bounceCount - 1);
					reflecC *= (fresnel * refraction + reflection);
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
					refracC *= ((1 - fresnel) * remainC * refraction);
				}
				currentC += reflecC;
				currentC += refracC;

				c += currentC;
			}
			for (int i = 0; i < lightsUsing.size(); i++) {
				Color currentC = Color(0, 0, 0);
				Light* l = lightsUsing[i];
				if (l->IsAmbient()) {
					currentC += l->Illuminate(hInfo.p, hInfo.N) * diffuse;
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

					float cosF = lDir.GetNormalized().Dot(nDir.GetNormalized());
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

					currentC += diff * diffuse;
					currentC += spec * specular;
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
						hInfo.front = true;
						return true;
					}
				}
				else if (k2 >= bias) {
					if (k2 < hInfo.z) {
						hInfo.z = k2;
						hInfo.p = ray.p + ray.dir * k2;
						hInfo.N = hInfo.p;
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