// Raytracing Prj1 by Xipeng Wang
#pragma once
#include <GL\freeglut.h>
#include <omp.h>
#include "scene.h"
#include "objects.h"
#include "tinyxml.h"
#include "materials.h"
#include "lights.h"
#include "cyTimer.h"
#include "myHelper.h"

#include <stdio.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <stack>

// best bias for prj6
//#define bias 0.00095f
// temp bias for prj7
#define bias 0.00005f
#define longDis 10000.0f
#define e_cons 2.718281828f
#define deltaOffset 0.005f
#define gamma_term 0.4545454545f
// normal using constants
//#define max_refLignt_bounce 3
//#define max_giLight_bounce 2
//#define max_variance 0.15f
//#define max_sampe_count 32
//#define min_halton_sample 8
//#define min_shadow_ray 4
//#define full_shadow_ray 16
//#define max_gi_sample 1
//#define min_surface_sample 1
//#define max_surface_sample 1

// pretty constants
//#define max_refLignt_bounce 5
//#define max_giLight_bounce 2
//#define max_variance 0.15f
//#define max_sampe_count 64
//#define min_halton_sample 8
//#define min_shadow_ray 8
//#define full_shadow_ray 32
//#define max_gi_sample 32
//#define min_surface_sample 1
//#define max_surface_sample 1

// path tracing constants
//#define max_light_count 10
#define max_refLignt_bounce 1
#define max_giLight_bounce 3
#define max_variance 0.15f
#define max_sampe_count 256

#define min_halton_sample 8
#define min_shadow_ray 1
#define full_shadow_ray 1
#define max_gi_sample 1
#define min_surface_sample 1
#define max_surface_sample 1


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

//char prjName[] = "test11";
char prjName[] = "prj12";
//char prjName[] = "prj11_pretty";
char prjSource[30];
char prjRender[30];
char prjZRender[30];
char prjCRender[30];

int LoadScene(char const *filename);
void ShowViewport();

float MyRandom(float max) {
	return (float)(rand()) / (float)(RAND_MAX)*max;
}

class CameraSpaceInfo {
public:
	Vec3f px, py, lc, right;

	void Init(Camera c) {
		if(!c.dir.IsZero())
			c.dir.Normalize();
		if (!c.up.IsZero())
			c.up.Normalize();
		right = c.dir.Cross(c.up);
		float ratio = (float)c.imgWidth / (float)c.imgHeight;
		float halfHeight = tanf(c.fov / 2 * Pi<float>() / 180);
		float halfWidth = halfHeight * ratio;
		lc = c.dir + c.up * halfHeight - right * halfWidth;
		px = right * halfWidth * 2 / c.imgWidth;
		py = c.up * -1 * halfHeight * 2 / c.imgHeight;
	}

	Vec3f GetPixelDir(float posX, float posY) {
		Vec3f d = lc + px * (posX + 0.5f) + py * (posY + 0.5f);
		return d;
	}
};
CameraSpaceInfo csInfo;

ImportanceLightSampler sampler;

struct JitterNode {
	float squaredWeight;
	Vec2f centerPos;
	Color variance;
	Color color;
	JitterNode* sibling;
	JitterNode* child;

	JitterNode() {
		squaredWeight = 1;
		centerPos = Vec2f(0, 0);
		variance = Color(1, 1, 1);
		color = Color(0, 0, 0);
		sibling = NULL;
		child = NULL;
	}
};

// Main Function
int main(int argc, char *argv[])
{
	omp_set_num_threads(24);

	strcpy_s(prjSource, prjName);
	strcat_s(prjSource, ".xml");

	strcpy_s(prjRender, "./Release/");
	strcat_s(prjRender, prjName);
	strcpy_s(prjZRender, prjRender);
	strcpy_s(prjCRender, prjRender);
	strcat_s(prjRender, ".png");
	strcat_s(prjZRender, "Z.png");
	strcat_s(prjCRender, "C.png");

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

Color GammaCorrection(Color sRGB) {
	float linearR = pow(sRGB.r, gamma_term);
	float linearG = pow(sRGB.g, gamma_term);
	float linearB = pow(sRGB.b, gamma_term);
	return Color(linearR, linearG, linearB);
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
		t1 = Max(t1, y1);
		t2 = Min(t2, y2);
	}
	if (t1 > z2 || t2 < z1) {
		return hit;
	}
	else {
		t1 = Max(t1, z1);
		t2 = Min(t2, z2);
	}

	hit = t1;

	return hit;
}

// Without Ray Differential
Color ShadePixel(Ray &ray, HitInfo &hInfo, Node *node, Vec2f relativePos) {
	Color c = Color(0, 0, 0);

	bool hResult = RayToNode(ray, hInfo, node);
	if (hResult) {
		const Node* hitNode = hInfo.node;
		//c = hitNode->GetMaterial()->Shade(ray, hInfo, lights, max_refLignt_bounce, max_giLight_bounce);
		c = hitNode->GetMaterial()->Shade(ray, hInfo, lights, sampler, max_refLignt_bounce, max_giLight_bounce);
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

// With 2 other Ray Differential
Color ShadePixelRayDiff(Ray &ray, Ray *drays , HitInfo &hInfo, Node *node, Vec2f relativePos) {
	Color c = Color(0, 0, 0);

	bool hResult = RayToNode(ray, hInfo, node);
	if (hResult) {
		HitInfo dInfoX = HitInfo();
		Vec3f uvwX = Vec3f();
		HitInfo dInfoY = HitInfo();
		Vec3f uvwY = Vec3f();
		bool xHit = RayToNode(drays[0], dInfoX, node);
		if (xHit) {
			if (dInfoX.node == hInfo.node) {
				uvwX = dInfoX.uvw - hInfo.uvw;
				uvwX /= deltaOffset;
			}
		}
		bool yHit = RayToNode(drays[1], dInfoY, node);
		if (yHit) {
			if (dInfoY.node == hInfo.node) {
				uvwY = dInfoY.uvw - hInfo.uvw;
				uvwY /= deltaOffset;
			}
		}
		hInfo.duvw[0] = uvwX;
		hInfo.duvw[1] = uvwY;
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

Color24 CalculatePixelHit(float posX, float posY) {
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

Color CalculatePixelColor(float posX, float posY) {
	Color c = Color();

	Vec3f rayp = camera.pos;
	Vec3f rayd = csInfo.GetPixelDir(posX, posY);

	Ray r = Ray(rayp, rayd);

	// diff by passing in right and up vectors
	r.diffRight = csInfo.px;
	r.diffUp = -csInfo.py;

	r.Normalize();

	HitInfo h = HitInfo();

	float relativeX = (float)posX / (float)camera.imgWidth;
	float relativeY = (float)posY / (float)camera.imgHeight;

	// without sending real rays ray differential
	c = ShadePixel(r, h, &rootNode, Vec2f(relativeX, relativeY));

	// with simple sending 2 new rays ray differential
	//Ray* dRays = new Ray[2];

	//Ray dRayX = Ray(rayp, csInfo.GetPixelDir(posX + deltaOffset, posY));
	//Ray dRayY = Ray(rayp, csInfo.GetPixelDir(posX, posY + deltaOffset));

	//dRayX.Normalize();
	//dRayY.Normalize();

	//dRays[0] = dRayX;
	//dRays[1] = dRayY;

	//c = ShadePixelRayDiff(r, dRays, h, &rootNode, Vec2f(relativeX, relativeY));
	return c;
}

Color CalculatePixelColorDepth(Vec2f inPos, Vec2f outPos) {
	Color c = Color();
	
	Vec3f rayp = camera.pos + csInfo.right * inPos.x + camera.up * inPos.y;
	Vec3f focalDir = csInfo.GetPixelDir(outPos.x, outPos.y);
	focalDir *= camera.focaldist;
	Vec3f rayd = focalDir - rayp;

	//Ray r = Ray(rayp, rayd);
	Ray r = Ray(rayp, focalDir - rayp + camera.pos);

	r.diffRight = csInfo.px * focalDir;
	r.diffUp = csInfo.py * focalDir;

	r.Normalize();

	HitInfo h = HitInfo();
	float relativeX = outPos.x / (float)camera.imgWidth;
	float relativeY = outPos.y / (float)camera.imgHeight;

	c = ShadePixel(r, h, &rootNode, Vec2f(relativeX, relativeY));

	return c;
}

Color24 CalculatePixelColor24(float posX, float posY) {
	Color c = CalculatePixelColor(posX, posY);
	Color24 color = Color24(c);
	return color;
}

// non recursive way of dynamic sampling (can use almost all memory)
uint8_t SamplePixelColorJitteredAdaptive(float posX, float posY, float sampleSize, Color24& color) {
	std::stack<JitterNode*> s;
	JitterNode* root = new JitterNode();
	root->centerPos = Vec2f(posX, posY);
	root->color = CalculatePixelColor(posX, posY);
	root->squaredWeight = 1;
	JitterNode* p = root;
	bool newNodes = false;

	uint8_t totalSample = 1;
	uint8_t lastSample = 1;

	while (root->variance.Sum() > max_variance && totalSample < max_sampe_count) {
		lastSample = totalSample;
		p = root;
		while (!s.empty() || p != NULL) {
			while (p != NULL) {
				if (p->child == NULL) {
					JitterNode* newChild = NULL;
					JitterNode* lastChild = NULL;
					float childWeight = p->squaredWeight / 2;
					float childSize = sampleSize * childWeight;

					Color newColorP = Color(0, 0, 0);
					Color newVarianceP = Color(0, 0, 0);

					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < 2; j++) {
							newChild = new JitterNode();
							Vec2f jitter = Vec2f(MyRandom(childSize), MyRandom(childSize));
							Vec2f newPos = p->centerPos + jitter * childSize + Vec2f(i - 1, j - 1) * childSize;
							Vec2f newCenter = p->centerPos + Vec2f(i - 0.5f, j - 0.5f) * childSize;
							
							newChild->centerPos = newCenter;
							newChild->color = CalculatePixelColor(newPos.x, newPos.y);
							newChild->squaredWeight = childWeight;

							newColorP += newChild->color;

							if (lastChild == NULL) {
								p->child = newChild;
							}
							else {
								lastChild->sibling = newChild;
							}

							lastChild = newChild;
						}
					}
					totalSample += 3;

					p = NULL;
				}
				else {
					s.push(p);
					p = p->child;
				}
			}
			if (!s.empty()) {
				p = s.top();
				s.pop();
				p = p->sibling;
			}
		}

		p = root;
		Color finalColor = Color(0, 0, 0);
		float totalWeight = 0;
		while (!s.empty() || p != NULL) {
			while (p != NULL) {
				if (p->child == NULL) {
					float currentWeight = (p->squaredWeight * p->squaredWeight);
					finalColor += p->color * currentWeight;
					totalWeight += currentWeight;
				}
				s.push(p);
				p = p->child;
			}
			if (!s.empty()) {
				p = s.top();
				s.pop();
				p = p->sibling;
			}
		}

		if (totalWeight != 0) {
			finalColor /= totalWeight;
		}

		root->color = finalColor;
		
		p = root;
		Color finalVariance = Color(0, 0, 0);
		while (!s.empty() || p != NULL) {
			while (p != NULL) {
				if (p->child == NULL) {
					Color currentVariance = root->color - p->color;
					currentVariance *= currentVariance;
					currentVariance *= (p->squaredWeight * p->squaredWeight);
					finalVariance += currentVariance;
				}
				s.push(p);
				p = p->child;
			}
			if (!s.empty()) {
				p = s.top();
				s.pop();
				p = p->sibling;
			}
		}
		finalVariance = Color(sqrt(finalVariance.r), sqrt(finalVariance.g), sqrt(finalVariance.b));
		root->variance = finalVariance;
	}

	color = Color24(root->color);

	p = root;
	while (!s.empty() || p != NULL) {
		while (p != NULL) {
			s.push(p);
			p = p->child;
		}
		if (!s.empty()) {
			p = s.top();
			s.pop();
			JitterNode* destroy = p;
			p = p->sibling;
			delete destroy;
		}
	}

	return totalSample;
}

uint8_t SamplePixelHaltonAdaptive(float posX, float posY, Color24& color) {
	Color finalColor = Color(0, 0, 0);
	Color previousColor = Color(0, 0, 0);
	uint8_t sampleCount = 0;
	float xRand = MyRandom(1);
	float yRand = MyRandom(1);
	Color variance = Color(0, 0, 0);
	while (sampleCount < min_halton_sample) {
		sampleCount += 1;
		float x = Halton(sampleCount, 2);
		float y = Halton(sampleCount, 3);
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		x += xRand;
		y += yRand;
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		Color currentColor = CalculatePixelColor(posX + x, posY + y);
		previousColor = finalColor;
		finalColor = previousColor * (sampleCount - 1) / sampleCount + currentColor / sampleCount;
		variance = variance + ((currentColor - previousColor)*(currentColor - finalColor) - variance) / sampleCount;
	}

	while (
		(sqrt(variance.r) + sqrt(variance.g) + sqrt(variance.b)) > 
		max_variance && sampleCount < max_sampe_count) {
		sampleCount += 1;
		float x = Halton(sampleCount, 2);
		float y = Halton(sampleCount, 3);
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		x += xRand;
		y += yRand;
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		Color currentColor = CalculatePixelColor(posX + x, posY + y);
		previousColor = finalColor;
		finalColor = previousColor * (sampleCount - 1) / sampleCount + currentColor / sampleCount;
		variance = variance + ((currentColor - previousColor) * (currentColor - finalColor) - variance) / sampleCount;
	}

	finalColor = GammaCorrection(finalColor);
	color = Color24(finalColor);

	return sampleCount;
}

uint8_t SamplePixelHalton(float posX, float posY, float sampleCount, Color24& color) {
	Color finalColor = Color(0, 0, 0);
	float xRand = MyRandom(1);
	float yRand = MyRandom(1);
	for (int i = 0; i < sampleCount; i++) {
		float x = Halton(i, 2);
		float y = Halton(i, 3);
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		x += xRand;
		y += yRand;
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		Color currentColor = CalculatePixelColor(posX + x, posY + y);
		finalColor += currentColor;
	}
	finalColor /= sampleCount;
	finalColor = GammaCorrection(finalColor);
	color = Color24(finalColor);
	return sampleCount;
}

uint8_t SamplePixelDepthofField(float posX, float posY, float apertureSampleCount, Color24& color) {
	Color finalColor = Color(0, 0, 0);

	float R = camera.dof;

	float xRand = MyRandom(1);
	float yRand = MyRandom(1);
	for (int i = 0; i < apertureSampleCount; i++) {
		float r = MyRandom(R * R);
		r = sqrt(r);
		float s = MyRandom(2 * Pi<float>());
		Vec2f currentIn = Vec2f(cosf(s), sinf(s));
		currentIn *= r;

		float x = Halton(i, 2);
		float y = Halton(i, 3);
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		x += xRand;
		y += yRand;
		if (x > 0.5f) x -= 1;
		if (y > 0.5f) y -= 1;
		Vec2f currentOut = Vec2f(posX + x, posY + y);

		Color currentColor = CalculatePixelColorDepth(currentIn, currentOut);
		finalColor += currentColor;
	}

	finalColor /= apertureSampleCount;
	color = Color24(finalColor);
	return apertureSampleCount;
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
	uint8_t* sample = renderImage.GetSampleCount();
	
	timer.Start();
	// build bvh here

	uint8_t sampleCount = 0;

	// initialization
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			sample[j * width + i] = 0;
			img[j * width + i] = Color24(0,0,0);
			zBuffer[j * width + i] = 0;
		}
	}

	// generate light source intensity map
	sampler.Init();
	for (int i = 0; i < lights.size(); i++) {
		if (lights[i]->IsPoint()) {
			sampler.PushValue(lights[i]->OverallIntensity(), i);
		}
	}
	sampler.RecalculateAll();

#pragma omp parallel for
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			// with out antialiasing
			//Color24 c = CalculatePixelColor(i, j);
			Color24 c = Color24();
			//uint8_t currentSampleCount = SamplePixelColorJitteredAdaptive(i, j, 1.0f, c);
			uint8_t currentSampleCount = SamplePixelHalton(i, j, max_sampe_count, c);
			//uint8_t currentSampleCount = SamplePixelHaltonAdaptive(i, j, c);
			//uint8_t currentSampleCount = SamplePixelDepthofField(i, j, 32, c);
			sample[j * width + i] = currentSampleCount;
			img[j * width + i] = c;
			float z = CalculatePixelZ(i, j);
			zBuffer[j * width + i] = z;
		}
	}

	renderImage.ComputeZBufferImage();
	renderImage.ComputeSampleCountImage();

	timer.Stop();
	auto end = timer.GetAverage();
	printf("Elapsed time: %f\n", end);
}

void StopRender() {

}

// Outside Functions

// not using path tracing
Color MtlBlinn::Shade(
	Ray const &ray,
	const HitInfo &hInfo,
	const LightList &lightsUsing,
	int refBounceCount, int giBounceCount) const
{
	Color c = Color(0, 0, 0);
	if (refBounceCount > 0) {
		int surfaceSampleCount = 1;

		Vec3f inDir = ray.dir * -1;
		inDir.Normalize();
		Vec3f baseNDir = hInfo.N.GetNormalized();
		Vec3f coordDir1 = Vec3f(0, 0, 0);
		Vec3f coordDir2 = Vec3f(0, 0, 0);
		float randRadius = 0;

		if (reflectionGlossiness != 0 || refractionGlossiness != 0) {
			surfaceSampleCount = max_surface_sample;
			Vec3f randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
			randomDir.Normalize();
			coordDir1 = baseNDir.Cross(randomDir);
			while (coordDir1.Length() < 0.75f) {
				randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
				randomDir.Normalize();
				coordDir1 = baseNDir.Cross(randomDir);
			}
			coordDir1.Normalize();
			coordDir2 = baseNDir.Cross(coordDir1);

			if (reflectionGlossiness != 0) {
				randRadius = tanf(reflectionGlossiness);
			}
			else if (refractionGlossiness) {
				randRadius = tanf(refractionGlossiness);
			}
			randRadius = randRadius * randRadius;
		}

		for (int j = 0; j < surfaceSampleCount; j++) {
			randRadius = MyRandom(randRadius);

			float randAngle = MyRandom(2 * Pi<float>());
			randRadius = sqrt(randRadius);

			float randX = sinf(randAngle);
			float randY = cosf(randAngle);

			Vec3f nDir = baseNDir + (randX * coordDir1 + randY * coordDir2) * randRadius;
			nDir.Normalize();

			Vec3f reflecDir = 2 * nDir.Dot(inDir) * nDir - inDir;
			reflecDir.Normalize();

			if (hInfo.front == true) {
				if (reflection.GetColor().r > 0) {
					Color currentC = Color(0, 0, 0);

					Ray reflecRay = Ray(hInfo.p, reflecDir);
					reflecRay.diffRight = ray.diffRight - 2 * (ray.diffRight.Dot(nDir) * nDir);
					reflecRay.diffUp = ray.diffUp - 2 * (ray.diffUp.Dot(nDir) * nDir);

					HitInfo reflectHit = HitInfo();
					bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
					if (hResult) {
						const Node* hitNode = reflectHit.node;
						currentC = hitNode->GetMaterial()
							->Shade(reflecRay, reflectHit, lights, refBounceCount - 1);
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
					float fresnel = r0 + (1 - r0) * pow((1 - cosFi1), 5);
					float kt = refraction.GetColor().r;

					Ray reflecRay = Ray(hInfo.p, reflecDir);
					reflecRay.diffRight = ray.diffRight - 2 * (ray.diffRight.Dot(nDir) * nDir);
					reflecRay.diffUp = ray.diffUp - 2 * (ray.diffUp.Dot(nDir) * nDir);

					HitInfo reflectHit = HitInfo();
					bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
					if (hResult) {
						const Node* hitNode = reflectHit.node;
						reflecC = hitNode->GetMaterial()
							->Shade(reflecRay, reflectHit, lights, refBounceCount - 1);
						reflecC = fresnel * reflecC * refraction.GetColor();
					}
					else {
						reflecC = fresnel * environment.SampleEnvironment(reflecDir) * refraction.GetColor();
					}

					Ray refracRay = Ray(hInfo.p, refracDir);
					refracRay.diffRight = ray.diffRight;
					refracRay.diffUp = ray.diffUp;

					HitInfo refracHit = HitInfo();

					hResult = RayToNode(refracRay, refracHit, &rootNode, HIT_FRONT_AND_BACK);
					if (hResult) {
						const Node* hitNode = refracHit.node;
						refracC = hitNode->GetMaterial()
							->Shade(refracRay, refracHit, lights, refBounceCount - 1);
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
			}
			else {
				Color currentC = Color(0, 0, 0);

				float n1 = ior;
				float n2 = 1;

				Vec3f inDir = ray.dir * -1;
				inDir.Normalize();
				nDir = nDir * -1;

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
					refracRay.diffRight = ray.diffRight;
					refracRay.diffUp = ray.diffUp;

					HitInfo refracHit = HitInfo();

					bool hResult = RayToNode(refracRay, refracHit, &rootNode, HIT_FRONT_AND_BACK);
					if (hResult) {
						const Node* hitNode = refracHit.node;
						refracC = hitNode->GetMaterial()
							->Shade(refracRay, refracHit, lights, refBounceCount - 1);
						currentC += refracC;
					}
					else {
						refracC = environment.SampleEnvironment(refracDir) * refraction.GetColor();
						currentC += refracC;
					}
				}
				else {
					Ray reflecRay = Ray(hInfo.p, reflecDir);
					reflecRay.diffRight = ray.diffRight - 2 * (ray.diffRight.Dot(nDir) * nDir);
					reflecRay.diffUp = ray.diffUp - 2 * (ray.diffUp.Dot(nDir) * nDir);

					HitInfo reflectHit = HitInfo();

					bool hResult = RayToNode(reflecRay, reflectHit, &rootNode, HIT_FRONT_AND_BACK);
					if (hResult) {
						const Node* hitNode = reflectHit.node;
						reflecC = hitNode->GetMaterial()
							->Shade(reflecRay, reflectHit, lights, refBounceCount - 1);
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
		c /= surfaceSampleCount;

		if (hInfo.front == true) {
			// Direct illuminations
			for (int i = 0; i < lightsUsing.size(); i++) {
				Color currentC = Color(0, 0, 0);
				Light* l = lightsUsing[i];
				Color illu = l->Illuminate(hInfo.p, hInfo.N);
				Color diff = Color(0, 0, 0);
				Color spec = Color(0, 0, 0);

				Vec3f lDir = l->Direction(hInfo.p) * -1;
				Vec3f vDir = -ray.dir;
				Vec3f nDir = baseNDir;
				Vec3f hDir = (vDir + lDir);
				hDir = hDir / hDir.Length();
				float geoTerm = lDir.GetNormalized().Dot(nDir);

				if (l->IsAmbient()) {
					// discard when using GI
					//currentC += l->Illuminate(hInfo.p, hInfo.N) * diffuse.Sample(hInfo.uvw, hInfo.duvw);
				}
				else {
					if (geoTerm <= 0)
						geoTerm = 0;
					else {
						diff = illu * geoTerm;
						float specValue = hDir.Dot(nDir);
						if (specValue < 0)
							specValue = 0;
						else {
							specValue = pow(specValue, glossiness);
						}
						spec = illu * specValue;
					}

					currentC += diff * diffuse.Sample(hInfo.uvw, hInfo.duvw);
					currentC += spec * specular.Sample(hInfo.uvw, hInfo.duvw);
				}
				c += currentC;
			}
			// Indirect illuminations
			if (giBounceCount > 0) {
				// generate coordinates
				Vec3f zdir = baseNDir;
				Vec3f randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
				randomDir.Normalize();
				Vec3f xdir = zdir.Cross(randomDir);
				while (xdir.Length() < 0.75f) {
					randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
					randomDir.Normalize();
					xdir = zdir.Cross(randomDir);
				}
				xdir.Normalize();
				Vec3f ydir = zdir.Cross(xdir);

				// indirect diffuse
				Color giDiffuse = Color(0, 0, 0);
				Color giSpecular = Color(0, 0, 0);

				// uniform sampling
				{
					//for (int i = 0; i < max_gi_sample; i++) {
					//	Color currentC = Color(0, 0, 0);

					//	float x = MyRandom(1);
					//	float cosTheta = 1 - x;
					//	float sinTheta = sqrt(1 - cosTheta * cosTheta);
					//	float phi = MyRandom(2 * Pi<float>());

					//	float xWeight = sinTheta * cosf(phi);
					//	float yWeight = sinTheta * sinf(phi);
					//	float zWeight = cosTheta;

					//	Vec3f omegaDir = xWeight * xdir + yWeight * ydir + zWeight * zdir;
					//	Color gi = Color(0, 0, 0);

					//	Ray giDiffRay = Ray(hInfo.p, omegaDir);
					//	HitInfo giDiffHInfo = HitInfo();
					//	bool hResult = RayToNode(giDiffRay, giDiffHInfo, &rootNode, HIT_FRONT);
					//	if (hResult) {
					//		const Node* hitNode = giDiffHInfo.node;
					//		gi = hitNode->GetMaterial()->Shade(giDiffRay, giDiffHInfo, lights, refBounceCount - 1, giBounceCount - 1);
					//	}
					//	else {
					//		gi = environment.SampleEnvironment(omegaDir);
					//	}

					//	float geoTerm = omegaDir.Dot(baseNDir);
					//	if (geoTerm > 0.01f) {
					//		currentC = gi * diffuse.Sample(hInfo.uvw, hInfo.duvw);
					//	}
					//	giDiffuse += currentC;
					//}
					//giDiffuse = giDiffuse / max_gi_sample;
					//////giDiffuse *= (Pi<float>() * 2);
				}
				// cosine-weighted sampling
				{
					for (int i = 0; i < max_gi_sample; i++) {
						Color currentC = Color(0, 0, 0);
						float x = MyRandom(1);
						float cosTheta = sqrt(x);
						float phi = MyRandom(2 * Pi<float>());
						float zWeight = sqrt(1 - x);
						float xWeight = cosf(phi) * cosTheta;
						float yWeight = sinf(phi) * cosTheta;

						Vec3f omegaDir = xWeight * xdir + yWeight * ydir + zWeight * zdir;

						Color gi = Color(0, 0, 0);

						Ray giDiffRay = Ray(hInfo.p, omegaDir);
						HitInfo giDiffHInfo = HitInfo();
						bool hResult = RayToNode(giDiffRay, giDiffHInfo, &rootNode, HIT_FRONT);
						if (hResult) {
							const Node* hitNode = giDiffHInfo.node;
							gi = hitNode->GetMaterial()->Shade(giDiffRay, giDiffHInfo, lights, refBounceCount - 1, giBounceCount - 1);
						}
						else {
							gi = environment.SampleEnvironment(omegaDir);
						}

						float geoTerm = omegaDir.Dot(baseNDir);
						if (geoTerm > 0) {
							currentC += (gi * diffuse.Sample(hInfo.uvw, hInfo.duvw));
						}

						giDiffuse += currentC;
					}
					giDiffuse = giDiffuse / max_gi_sample;
					//giDiffuse *= (Pi<float>());
				}
				// indirect specular
				
				{

				}

				c = c + giDiffuse + giSpecular;
			}

			// emmisive
		}
	}
	return c;
}

// using path tracing
Color MtlBlinn::Shade(
	Ray const& ray,
	const HitInfo& hInfo,
	const LightList & lightsUsing,
	ImportanceLightSampler theSampler,
	int refBounceCount, int giBounceCount) const {
	Color c = Color(0, 0, 0);
	if (giBounceCount > 0) {
		Vec3f inDir = ray.dir * -1;
		inDir.Normalize();
		Vec3f nDir = hInfo.N.GetNormalized();
		Vec3f pos = hInfo.p;
		
		// not divide light intensity by probablity correctly
		int currentLightIndex = theSampler.RandomGet();
		PointLight* currentLight = (PointLight *)lightsUsing[currentLightIndex];

		// generate gi ray
		Vec3f diffSampleDir = Vec3f(0, 0, 0);
		Vec3f specSampleDir = Vec3f(0, 0, 0);
		Vec3f directLightDir = Vec3f(0, 0, 0);

		// sample based on cosine

		// sample based on brdf
		float phiIndirect = MyRandom(2 * Pi<float>());
		float xIndirect = MyRandom(1);
		float cosThetaIndirect = pow(xIndirect, 1 / (2 * glossiness + 1));
		//float p2 = pow(cosThetaIndirect, 2 * glossiness);
		//float p2 = pow(cosThetaIndirect, 2 * glossiness) * (glossiness + 2) / (2 * Pi<float>());

		bool samplingLightSource = true;
		float sampleWeighter = MyRandom(1);

		Vec3f reflecDir = 2 * nDir.Dot(inDir) * nDir - inDir;
		reflecDir.Normalize();
		Vec3f randomDir = Vec3f(0, 0, 0);

		// cosine weighted random dir
		{
			Vec3f coszdir = nDir;
			randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
			randomDir.Normalize();
			Vec3f cosxdir = coszdir.Cross(randomDir);
			while (cosxdir.Length() < 0.75f) {
				randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
				randomDir.Normalize();
				cosxdir = coszdir.Cross(randomDir);
			}
			cosxdir.Normalize();
			Vec3f cosydir = coszdir.Cross(cosxdir);

			Color currentC = Color(0, 0, 0);
			float x = MyRandom(1);
			float cosTheta = sqrt(x);
			float phi = MyRandom(2 * Pi<float>());
			float zWeight = sqrt(1 - x);
			float xWeight = cosf(phi) * cosTheta;
			float yWeight = sinf(phi) * cosTheta;

			diffSampleDir = xWeight * cosxdir + yWeight * cosydir + zWeight * coszdir;
		}

		// using brdf
		// generate direction
		{
			Vec3f zdir = reflecDir;
			randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
			randomDir.Normalize();
			Vec3f xdir = zdir.Cross(randomDir);
			while (xdir.Length() < 0.75f) {
				randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
				randomDir.Normalize();
				xdir = zdir.Cross(randomDir);
			}
			xdir.Normalize();
			Vec3f ydir = zdir.Cross(xdir);

			float sinThetaIndirect = sqrt(1 - cosThetaIndirect * cosThetaIndirect);
			specSampleDir = zdir * cosThetaIndirect + sinThetaIndirect * (xdir * cosf(phiIndirect) + ydir * sinf(phiIndirect));
		}

		Color DiffIllumination = Color(0, 0, 0);
		Color SpecIllumination = Color(0, 0, 0);

		Ray giRay = Ray(hInfo.p, diffSampleDir);
		HitInfo giHit = HitInfo();
		bool hResult = RayToNode(giRay, giHit, &rootNode, HIT_FRONT);

		if (hResult) {
			const Node* hitNode = giHit.node;
			DiffIllumination = hitNode->GetMaterial()->Shade(giRay, giHit, lights, theSampler, refBounceCount, giBounceCount - 1);
		}
		else {
			DiffIllumination = environment.SampleEnvironment(diffSampleDir);
		}

		giRay = Ray(hInfo.p, specSampleDir);
		giHit = HitInfo();
		hResult = RayToNode(giRay, giHit, &rootNode, HIT_FRONT);
		if (hResult) {
			const Node* hitNode = giHit.node;
			SpecIllumination = hitNode->GetMaterial()->Shade(giRay, giHit, lights, theSampler, refBounceCount, giBounceCount - 1);
		}
		else {
			SpecIllumination = environment.SampleEnvironment(specSampleDir);
		}
		//SpecIllumination /= cosThetaIndirect;

		float geoTerm = diffSampleDir.Dot(nDir);
		if (geoTerm < 0) {
			geoTerm = 0;
		}

		Vec3f hDir = inDir + specSampleDir;
		hDir.Normalize();
		float specValue = hDir.Dot(nDir);
		if (specValue < 0) {
			specValue = 0;
		}
		specValue = pow(specValue, glossiness);

		Color diffuseBrdf = diffuse.Sample(hInfo.uvw, hInfo.duvw);
		Color specularBrdf = specular.Sample(hInfo.uvw, hInfo.duvw) * specValue;

		Color indirectIllumination = DiffIllumination * diffuseBrdf +
			SpecIllumination * specularBrdf;

		Color directIllumination = Color(0, 0, 0);
		{
			directIllumination = currentLight->Illuminate(hInfo.p, hInfo.N);
			directLightDir = currentLight->Direction(hInfo.p) * -1;
			directLightDir.Normalize();

			float directGeo = directLightDir.Dot(nDir);
			if (directGeo < 0) {
				directGeo = 0;
			}
			Vec3f directHDir = inDir + directLightDir;
			directHDir.Normalize();
			float directSpec = directHDir.Dot(nDir);
			if (directSpec < 0) {
				directSpec = 0;
			}
			directSpec = pow(directSpec, glossiness);
			//Color directBrdf = diffuse.Sample(hInfo.uvw, hInfo.duvw) * directGeo / Pi<float>() +
			//	(directGeo + 2) / (Pi<float>() * 2) * directSpec * specular.Sample(hInfo.uvw, hInfo.duvw);
			Color directBrdf = diffuse.Sample(hInfo.uvw, hInfo.duvw) * directGeo +
				directSpec * specular.Sample(hInfo.uvw, hInfo.duvw);
			directIllumination = directIllumination * directBrdf;
		}

		// calculate base on the sample direction
		Color finalColor = directIllumination + indirectIllumination;
		//Color finalColor = specValue * Color(1, 1, 1);

		c = finalColor;
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



						hInfo.duvw[0] = Vec3f(0, 0, 0);
						hInfo.duvw[1] = Vec3f(0, 0, 0);

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

						hInfo.duvw[0] = Vec3f(0, 0, 0);
						hInfo.duvw[1] = Vec3f(0, 0, 0);

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

						hInfo.duvw[0] = Vec3f(0, 0, 0);
						hInfo.duvw[1] = Vec3f(0, 0, 0);

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

						hInfo.duvw[0] = Vec3f(0, 0, 0);
						hInfo.duvw[1] = Vec3f(0, 0, 0);

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

		if (t > hInfo.z) {
			return false;
		}

		if (hitSide == HIT_FRONT) {
			if (NdotD > 0) {
				return false;
			}
		}
		else if (hitSide == HIT_BACK) {
			if (NdotD < 0) {
				return false;
			}
		}

		hInfo.z = t;
		hInfo.p = x;
		hInfo.N = vertexN;
		hInfo.uvw = uvw;
		hInfo.front = NdotD < 0;

		hInfo.duvw[0] = Vec3f(0, 0, 0);
		hInfo.duvw[1] = Vec3f(0, 0, 0);

		return true;
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

		if (t > hInfo.z) {
			return false;
		}

		if (hitSide == HIT_FRONT) {
			if (NdotD > 0) {
				return false;
			}
		}
		else if (hitSide == HIT_BACK) {
			if (NdotD < 0) {
				return false;
			}
		}

		hInfo.z = t;
		hInfo.p = x;
		hInfo.N = n;
		// up-left corner uv
		// float u = (x.x - corner4.x) / 2;
		// float v = (corner4.y - x.y) / 2;
		// bottom-left corner uv
		float u = (x.x - corner1.x) / 2;
		float v = (x.y - corner1.y) / 2;
		hInfo.uvw = Vec3f(u, v, 0);
		hInfo.front = NdotD < 0;

		Vec3f D = Vec3f(d);
		D.Normalize();
		float T = -(p.Dot(n)) / D.Dot(n);

		// calculate px & py by using passed in right & up vectors
		Vec3f deltaDX = (d.Dot(d) * (ray.diffRight) - d.Dot(ray.diffRight) * d) / pow(d.Dot(d), 1.5f);
		Vec3f deltaDY = (d.Dot(d) * (ray.diffUp) - D.Dot(ray.diffUp) * d) / pow(d.Dot(d), 1.5f);

		//printf("diff right: %f\n", deltaDX.x * 10000);

		float deltaTX = -((T * deltaDX).Dot(n)) / (D.Dot(n));
		float deltaTY = -((T * deltaDY).Dot(n)) / (D.Dot(n));

		//printf("diff right: %f\n", deltaTX * 100000);
		//printf("diff right: %f\n", T);

		Vec3f deltaPX = T * deltaDX + D * deltaTX;
		Vec3f deltaPY = T * deltaDY + D * deltaTY;

		hInfo.duvw[0] = (deltaPX) / 2.0f;
		hInfo.duvw[1] = (deltaPY) / 2.0f;

		return true;
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

Color PointLight::Illuminate(Vec3f const& p, Vec3f const& N) const {
	Vec3f centerDir = position - p;
	Vec3f nCenterDir = centerDir / centerDir.Length();
	Vec3f randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
	randomDir.Normalize();
	Vec3f coordDir1 = nCenterDir.Cross(randomDir);
	while (coordDir1.Length() < 0.75f) {
		randomDir = Vec3f(MyRandom(1), MyRandom(1), MyRandom(1));
		randomDir.Normalize();
		coordDir1 = nCenterDir.Cross(randomDir);
	}
	coordDir1.Normalize();
	Vec3f coordDir2 = nCenterDir.Cross(coordDir1);

	bool colorSame = true;
	float firstShadow = 0;
	float sizeS = size * size;
	//sizeS = 0;
	float finalShadow = 0;

	for (int i = 0; i < min_shadow_ray; i++) {
		float randAngle = MyRandom(2 * Pi<float>());
		float randRadius = MyRandom(sizeS);
		randRadius = sqrt(randRadius);

		float randX = sinf(randAngle);
		float randY = cosf(randAngle);

		Vec3f currentDir = centerDir + (randX * coordDir1 + randY * coordDir2) * randRadius;
		float currentShadow = Shadow(Ray(p, currentDir), 1);
		if (currentDir.Dot(N) < 0) {
			currentShadow = 0;
		}

		if (size == 0) {
			return currentShadow * intensity;
		}

		if (i == 0) {
			firstShadow = currentShadow;
		}

		if (currentShadow != firstShadow) {
			colorSame = false;
		}

		finalShadow += currentShadow;
	}

	if (colorSame) {
		finalShadow /= min_shadow_ray;
	}
	else {
		for (int i = 0; i < full_shadow_ray - min_shadow_ray; i++) {
			float randAngle = MyRandom(2 * Pi<float>());
			float randRadius = MyRandom(sizeS);
			randRadius = sqrt(randRadius);

			float randX = sinf(randAngle);
			float randY = cosf(randAngle);

			Vec3f currentDir = centerDir + (randX * coordDir1 + randY * coordDir2) * randRadius;
			float currentShadow = Shadow(Ray(p, currentDir), 1);

			finalShadow += currentShadow;
		}

		finalShadow /= full_shadow_ray;
	}

	return finalShadow * intensity;
}