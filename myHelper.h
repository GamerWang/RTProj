#pragma once
#include <stdio.h>

#define max_light_count 10

class ImportanceLightSampler {
public:
	void Init() {
		int lightsCount = 0;

		for (int i = 0; i < max_light_count; i++) {
			importanceList[i] = 0;
			lightIndex[i] = -1;
		}
	}

	void PushValue(float intensity, int index) {
		if (lightsCount > 0) {
			importanceList[lightsCount] = importanceList[lightsCount - 1] + intensity;
		}
		else {
			importanceList[lightsCount] = intensity;
		}
		lightIndex[lightsCount] = index;
		lightsCount++;
	}

	void RecalculateAll() {
		if (lightsCount > 0) {
			float sumImportance = importanceList[lightsCount - 1];
			for (int i = 0; i < lightsCount; i++) {
				importanceList[i] /= sumImportance;
			}
		}
	}

	int RandomGet() {
		float randomValue = (float)(rand()) / (float)(RAND_MAX);
		for (int i = 0; i < lightsCount; i++) {
			if (randomValue < importanceList[i]) {
				return lightIndex[i];
			}
		}
	}

private:
	float importanceList[max_light_count];
	int lightIndex[max_light_count];
	int lightsCount;
};