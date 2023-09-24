#include "mesh_reduction.h"

#include <Eigen/Sparse>
#include <igl/collapse_edge.h>
#include <igl/min_heap.h>
#include <igl/parallel_for.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/edge_flaps.h>
#include <igl/shortest_edge_and_midpoint.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory_resource>
#include <set>
#include <thread>
#include <vector>

#include "assimp_helper.h"

#ifdef QEM_USE_UHEMESH
#include <QEMUHEMesh.h>
#endif

#ifdef QEM_USE_OPENMESH
#include <QEMOpenMesh.h>
#endif



// vert idx, neighbor faces idx
std::unordered_map<int, std::vector<int>> findVertNeighborFaces(const Mesh &mesh)
{
	int len = mesh.F.rows();
	std::unordered_map<int, std::vector<int>> map;
	for (int i = 0; i < len; ++i)
	{
		int idx1 = mesh.F(i, 0);
		int idx2 = mesh.F(i, 1);
		int idx3 = mesh.F(i, 2);
		map[idx1].push_back(i);
		map[idx2].push_back(i);
		map[idx3].push_back(i);
	}
	return map;
};


//计算每一个面对应的Kp矩阵
Eigen::Matrix4d countFaceKPMatrix(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3) {
	Eigen::Vector3d n1(v1 - v2);
	Eigen::Vector3d n2(v2 - v3);
	Eigen::Vector3d normal = n1.cross(n2);
	normal.normalize();
	double d = -(normal.x() * v1.x() + normal.y() * v1.y() + normal.z() * v1.z());
	Eigen::Vector4d p(normal.x(), normal.y(), normal.z(), d);
	Eigen::Matrix4d mat = p*p.transpose();
	return mat;
}


Mesh MeshReduction::Reduction(const Mesh &mesh, float ratio)
{	
	int verNum = mesh.V.rows();
	std::unordered_map<int, vector<int>> vertNeighborFaces = findVertNeighborFaces(mesh);

	for (int i = 0; i < verNum; ++i)
	{
		Eigen::Matrix4d vertAddMatrix;
		vertAddMatrix.setZero();
		for (int idx : vertNeighborFaces[i])
		{
			Eigen::Vector3d v1(mesh.V(mesh.F(idx, 0), 0), mesh.V(mesh.F(idx, 0), 1), mesh.V(mesh.F(idx, 0), 2));
			Eigen::Vector3d v2(mesh.V(mesh.F(idx, 1), 0), mesh.V(mesh.F(idx, 1), 1), mesh.V(mesh.F(idx, 1), 2));
			Eigen::Vector3d v3(mesh.V(mesh.F(idx, 2), 0), mesh.V(mesh.F(idx, 2), 1), mesh.V(mesh.F(idx, 2), 2));
			vertAddMatrix += countFaceKPMatrix(v1,v2,v3);
		}



	}




	//-------------
	QEM qem;
	qem.ImportMesh(mesh);
	// TODO: implement me
	int targetCount = std::ceil(mesh.F.rows() * ratio);
	qem.DoSimplification(QEM::SimplificationMode::Position, targetCount);
	Mesh outputMesh = qem.ExportMesh();
	outputMesh.texture = mesh.texture;
	return outputMesh;
}

Mesh MeshReduction::ReductionWithAppearancePresentation(const Mesh &mesh, float ratio)
{
	QEM qem;
	qem.ImportMesh(mesh);
	// TODO: implement me
	int targetCount = std::ceil(mesh.F.rows() * ratio);
	qem.DoSimplification(QEM::SimplificationMode::PositionNormal, targetCount);

	Mesh outputMesh = qem.ExportMesh();
	outputMesh.texture = mesh.texture;
	return outputMesh;
}

Mesh MeshReduction::ReductionWithUVMap(const Mesh &mesh, float ratio)
{
	QEM qem;
	qem.ImportMesh(mesh);
	// TODO: implement me
	int targetCount = std::ceil(mesh.F.rows() * ratio);
	qem.DoSimplification(QEM::SimplificationMode::PositionNormalUV, targetCount);

	Mesh outputMesh = qem.ExportMesh();
	outputMesh.texture = mesh.texture;
	return outputMesh;
}
