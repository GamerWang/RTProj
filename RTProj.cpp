// Raytracing Prj1 by Xipeng Wang
#include <GL\freeglut.h>
#include "scene.h"
#include "objects.h"
#include "tinyxml.h"

#include <stdio.h>
#include <chrono>
#include <iostream>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;

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
	LoadScene(".\\prj1.xml");
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
					hInfo.z = k1 < hInfo.z ? k1 : hInfo.z;
					hInfo.front = true;
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

bool RayToNode(Ray const &ray, HitInfo &hInfo, Node *node) {
	Object *obj = node->GetNodeObj();
	Ray r = node->ToNodeCoords(ray);

	bool hResult = false;

	if (obj) {
		hResult = obj->IntersectRay(r, hInfo);
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

Color24 CalculatePixelHit(int posX, int posY) {
	Color c = Color();

	Vec3f rayp = camera.pos;
	Vec3f rayd = csInfo.GetPixelDir(posX, posY);

	Ray r = Ray(rayp, rayd);
	HitInfo h = HitInfo();

	//Temp: only hit detect
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
			Color24 c = CalculatePixelHit(i, j);
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