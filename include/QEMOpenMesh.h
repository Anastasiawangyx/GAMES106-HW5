#pragma once

#include "QEMDebug.h"
#include "assimp_helper.h"
#include <vector>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriConnectivity.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>


class QEM
{
  public:
	enum SimplificationMode
	{
		Position,
		PositionNormal,
		PositionNormalUV
	};

	struct MeshTraits : public OpenMesh::DefaultTraits
	{
		using Point      = OpenMesh::Vec3d;
		using Normal     = OpenMesh::Vec3d;
		using TexCoord2D = OpenMesh::Vec2d;
		VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::TexCoord2D);
		FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
		EdgeAttributes(OpenMesh::Attributes::Status);
		HalfedgeAttributes(OpenMesh::Attributes::Status);
	};

	using QEMMesh = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;
  private:
	QEMMesh heMesh;

  public:
	inline void ImportMesh(const Mesh &mesh)
	{
		std::vector<QEMMesh::VertexHandle> vhs;

		for (int i = 0; i < mesh.V.rows(); i++)
		{
			auto vh = heMesh.add_vertex(OpenMesh::Vec3d(mesh.V(i, 0), mesh.V(i, 1), mesh.V(i, 2)));
			heMesh.set_normal(vh, OpenMesh::Vec3d(mesh.N(i, 0), mesh.N(i, 1), mesh.N(i, 2)));
			heMesh.set_texcoord2D(vh, OpenMesh::Vec2d(mesh.UV(i, 0), mesh.UV(i, 1)));
			vhs.push_back(vh);
		}

		for (int i = 0; i < mesh.F.rows(); i++)
		{
			heMesh.add_face(vhs[mesh.F(i, 0)], vhs[mesh.F(i, 1)], vhs[mesh.F(i, 2)]);
		}
	}

	inline Mesh ExportMesh()
	{
		Mesh outputMesh;
		outputMesh.V.resize(heMesh.n_vertices(), 3);
		outputMesh.F.resize(heMesh.n_faces(), 3);
		outputMesh.N.resize(heMesh.n_vertices(), 3);
		outputMesh.UV.resize(heMesh.n_vertices(), 2);

		for (int i = 0; i < heMesh.n_vertices(); i++)
		{
			auto vh            = heMesh.vertex_handle(i);
			outputMesh.V(i, 0) = heMesh.point(vh)[0];
			outputMesh.V(i, 1) = heMesh.point(vh)[1];
			outputMesh.V(i, 2) = heMesh.point(vh)[2];
			outputMesh.N(i, 0) = heMesh.normal(vh)[0];
			outputMesh.N(i, 1) = heMesh.normal(vh)[1];
			outputMesh.N(i, 2) = heMesh.normal(vh)[2];
			outputMesh.UV(i, 0) = heMesh.texcoord2D(vh)[0];
			outputMesh.UV(i, 1) = heMesh.texcoord2D(vh)[1];
		}

		for (int i = 0; i < heMesh.n_faces(); i++)
		{
			auto fh  = heMesh.face_handle(i);
			auto it  = heMesh.fv_cwbegin(fh);
			int  idx = 0;
			while (it != heMesh.fv_cwend(fh))
			{
				assert(idx < 3);
				outputMesh.F(i, idx) = it->idx();
				idx++;
				it++;
			}
		}

		return outputMesh;
	}

	// 计算每一个面对应的Kp矩阵
	Eigen::Matrix4d countFaceKPMatrix(QEMMesh::FaceHandle fh)
	{
		auto normal = heMesh.normal(fh);
		normal.normalize();
		QEMMesh::FaceVertexIter fv_it  = heMesh.fv_iter(fh);
		QEMMesh::VertexHandle   vert   = *fv_it;
		QEMMesh::Point          position = heMesh.point(vert);
		double                  d        = -(normal[0] * position[0] + normal[1] * position[1] + normal[2] * position[2]);
		Eigen::Vector4d         p(normal[0], normal[1], normal[2], d);
		Eigen::Matrix4d         mat = p * p.transpose();
		return mat;
	}

	inline bool DoSimplification(SimplificationMode mode, int targetFaceCount)
	{
		QEM_DEBUG("DoSimplification(mode=%d, targetFaceCount=%d)", mode, targetFaceCount);
		// TODO: check validity
		std::cout << "hello openmesh qem" << endl;
		heMesh.request_face_status();
		heMesh.request_edge_status();
		heMesh.request_vertex_status();
		OpenMesh::VPropHandleT<Eigen::Matrix4d> vertMatrixProp;        
		heMesh.add_property(vertMatrixProp);
		//遍历顶点求kp矩阵
		for (int i = 0; i < heMesh.n_vertices(); i++)
		{
			Eigen::Matrix4d mat;
			mat.setZero();
			auto vh = heMesh.vertex_handle(i);
			for (QEMMesh::VertexFaceIter vf_it = heMesh.vf_iter(vh); vf_it.is_valid(); ++vf_it)
			{
				auto fh = heMesh.face_handle(vf_it);
				mat += countFaceKPMatrix(fh);
			}
			heMesh.property(vertMatrixProp, vh) = mat;
		}

		// error new-vertex-position
		OpenMesh::EPropHandleT<std::pair<double,Eigen::Vector4d>> edgeMatrixProp;
		heMesh.add_property(edgeMatrixProp);
		std::priority_queue<std::pair<double, QEMMesh::EdgeHandle>, vector<std::pair<double, QEMMesh::EdgeHandle>>, greater<std::pair<double, QEMMesh::EdgeHandle>>> heap;
		// 遍历每一条边
		for (QEMMesh::EdgeIter e_it = heMesh.edges_begin(); e_it != heMesh.edges_end(); ++e_it)
		{
			QEMMesh::EdgeHandle edge = *e_it;

			// 获取边的两个顶点
			QEMMesh::HalfedgeHandle halfedge1 = heMesh.halfedge_handle(edge, 0);
			

			QEMMesh::VertexHandle vertex1 = heMesh.to_vertex_handle(halfedge1);
			QEMMesh::VertexHandle vertex2 = heMesh.from_vertex_handle(halfedge1);

			// 在这里使用 vertex1 和 vertex2，它们分别是边的两个顶点
			Eigen::Matrix4d edgeMatrix = heMesh.property(vertMatrixProp, vertex1) + heMesh.property(vertMatrixProp, vertex2);
			Eigen::Matrix4d mat;
			mat << edgeMatrix(0, 0), edgeMatrix(0, 1), edgeMatrix(0, 2), edgeMatrix(0, 3),
			       edgeMatrix(0, 1), edgeMatrix(1, 1), edgeMatrix(1, 2), edgeMatrix(1, 3),
			       edgeMatrix(0, 2), edgeMatrix(1, 2), edgeMatrix(2, 2), edgeMatrix(2, 3),
			       0, 0, 0, 1;
			Eigen::Vector4d b(0, 0, 0, 1);
			Eigen::Vector4d x;
			x            = mat.colPivHouseholderQr().solve(b);
			double error = x.transpose() * edgeMatrix * x;
			heMesh.property(edgeMatrixProp, edge) = make_pair(error, x);
			heap.push(make_pair(error, edge));
		}

		std::cout << "size : " << heap.size()<< endl;
		
		//while (!heap.empty())
		//{
		//	auto targetEdge = heap.top();
		//	// TODO 有效性检测
		//	heap.pop();
		//	std::cout << targetEdge.first << endl
		//	          << targetEdge.second.is_valid() << endl;
		//}
		
		auto targetEdge = heap.top();
		// TODO 有效性检测
		heap.pop();
		// 坍缩边
		QEMMesh::HalfedgeHandle heh1 = heMesh.halfedge_handle(targetEdge.second, 0);

		heMesh.collapse(heh1);

		//while (!heap.empty() && heMesh.n_vertices() > targetFaceCount)
		//{
		//	
		//	// 新顶点
		//	
		//	break;


		//}

		return true;
	}
};