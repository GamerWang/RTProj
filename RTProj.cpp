﻿// Raytracing Prj1 by Xipeng Wang
#include <GL\freeglut.h>
#include "scene.h"
#include "objects.h"
#include "tinyxml.h"
#include "materials.h"

#include <stdio.h>
#include <chrono>
#include <iostream>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
LightList lights;
MaterialList materials;

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

int main(int argc, char *argv[])
{
	LoadScene(".\\prj2.xml");
	ShowViewport();
	return 0;
}

bool Sphere::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide) const
{
	if (hitSide == HIT_NONE) {
		return false;
	}
	else {
		float delta = 4 * ray.p.Dot(ray.dir) * ray.dir.Dot(ray.p) - 4 * ray.dir.Dot(ray.dir) * (ray.p.Dot(ray.p) - 1);
		if (delta < 0) {
			return false;
		}
		else {
			float k1 = (-2 * ray.p.Dot(ray.dir) - sqrt(delta)) / (ray.dir.Dot(ray.dir));
			float k2 = (-2 * ray.p.Dot(ray.dir) + sqrt(delta)) / (ray.dir.Dot(ray.dir));
			// Temp: only hit detect
			if (hitSide == HIT_FRONT) {
				if (k1 >= 0) {
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
				if (k2 >= 0) {
					hInfo.z = -k2;
					hInfo.front = false;
					return true;
				}
				else {
					return false;
				}
			}
			else if (hitSide == HIT_FRONT_AND_BACK) {
				// Temp: deprecated
				return false;
			}
		}
	}
}

Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights) const {
	Color c = Color(0, 0, 0);

	for (int i = 0; i < lights.size(); i++) {
		Color currentC = Color(0, 0, 0);
		Light* l = lights[i];
		if (l->IsAmbient()) {
			currentC += l->Illuminate(hInfo.p, hInfo.N);
		}
		else {
			Color illu = l->Illuminate(hInfo.p, hInfo.N);
			
		}
		c += currentC;
	}
	return c;
}

bool RayToNode(Ray &ray, HitInfo &hInfo, Node *node) {
	Object *obj = node->GetNodeObj();
	Ray r = node->ToNodeCoords(ray);

	bool hResult = false;

	if (obj) {
		hResult = obj->IntersectRay(r, hInfo);
	}

	if (hResult) {
		hInfo.node = node;
		ray.dir = r.dir;
		ray.p = r.p;
	}

	for (int i = 0; i < node->GetNumChild(); i++) {
		hResult = hResult || RayToNode(r, hInfo, node->GetChild(i));
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
		c = hitNode->GetMaterial()->Shade(ray, hInfo, lights);
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