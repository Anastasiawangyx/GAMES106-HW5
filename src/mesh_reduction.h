#pragma once

#include "assimp_helper.h"
#include <unordered_map>
#include <vector>
using namespace std;
class MeshReduction
{

public:
    


	// Task 1: Mesh reduction
	static Mesh Reduction(string file, const Mesh& mesh, float ratio);

	// Task 2: Mesh reduction with appearance presentation
	static Mesh ReductionWithAppearancePresentation(string file, const Mesh& mesh, float ratio);

	// Task3: Mesh reduction with UV map
	static Mesh ReductionWithUVMap(string file, const Mesh& mesh, float ratio);
};

