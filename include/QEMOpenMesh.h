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
	inline void ImportMesh(string file, const Mesh &mesh)
	{
		//直接将输入模型的文件读入进行处理
		cout << "file dir : " << file << endl;
		if (!OpenMesh::IO::read_mesh(heMesh, file))
		{
			std::cerr << "read error\n";
			exit(1);
		}
		else
		{
			cout << "read file success!" << endl;
		}
		//直接从顶点构造simple的四个三角形
		//std::vector<QEMMesh::VertexHandle> vhandle(6);

		//vhandle[0] = heMesh.add_vertex(QEMMesh::Point(0, 0, 1));
		//vhandle[1] = heMesh.add_vertex(QEMMesh::Point(1, 0, 1));
		//vhandle[2] = heMesh.add_vertex(QEMMesh::Point(2, 0, 1));
		//vhandle[3] = heMesh.add_vertex(QEMMesh::Point(0, 0, 0));
		//vhandle[4] = heMesh.add_vertex(QEMMesh::Point(1, 0, 0));
		//vhandle[5] = heMesh.add_vertex(QEMMesh::Point(2, 0, 0));

		//auto fh1 = heMesh.add_face(vhandle[1], vhandle[3], vhandle[0]);
		//auto fh2 = heMesh.add_face(vhandle[2], vhandle[4], vhandle[1]);
		//auto fh3 = heMesh.add_face(vhandle[1], vhandle[4], vhandle[3]);
		//auto fh4 = heMesh.add_face(vhandle[2], vhandle[5], vhandle[4]);

		// origin method
		//std::vector<QEMMesh::VertexHandle> vhs;

		//for (int i = 0; i < mesh.V.rows(); i++)
		//{
		//	auto vh = heMesh.add_vertex(OpenMesh::Vec3d(mesh.V(i, 0), mesh.V(i, 1), mesh.V(i, 2)));
		//	heMesh.set_normal(vh, OpenMesh::Vec3d(mesh.N(i, 0), mesh.N(i, 1), mesh.N(i, 2)));
		//	heMesh.set_texcoord2D(vh, OpenMesh::Vec2d(mesh.UV(i, 0), mesh.UV(i, 1)));
		//	vhs.push_back(vh);
		//}

		//for (int i = 0; i < mesh.F.rows(); i++)
		//{
		//	heMesh.add_face(vhs[mesh.F(i, 0)], vhs[mesh.F(i, 1)], vhs[mesh.F(i, 2)]);
		//}
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

	void updataVertexNormal(QEMMesh::VertexHandle vertex)
	{
		// 更新顶点周围的法线
		QEMMesh::Normal updated_normal(0.0, 0.0, 0.0);        // 初始化为零向量

		// 遍历与顶点相邻的面片
		for (QEMMesh::VertexFaceIter vf_it = heMesh.vf_iter(vertex); vf_it.is_valid(); ++vf_it)
		{
			QEMMesh::FaceHandle face = *vf_it;

			// 计算每个面片的法线
			QEMMesh::Normal face_normal = heMesh.normal(face);

			// 累加法线以计算平均法线
			updated_normal += face_normal;
		}

		// 标准化法线以确保它是单位向量
		if (updated_normal.norm() > 0.0)
		{
			updated_normal.normalize();
		}

		// 更新顶点的法线
		heMesh.set_normal(vertex, updated_normal);
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
		for (auto& vh : heMesh.vertices())
		{
			Eigen::Matrix4d mat;
			mat.setZero();
			for (auto& fh : vh.faces())
			{
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
			QEMMesh::HalfedgeHandle halfedge2 = heMesh.halfedge_handle(edge, 1);

			QEMMesh::VertexHandle vertex1 = heMesh.to_vertex_handle(halfedge1);
			QEMMesh::VertexHandle vertex2 = heMesh.from_vertex_handle(halfedge1);

			// 在这里使用 vertex1 和 vertex2，它们分别是边的两个顶点
			Eigen::Matrix4d edgeMatrix = heMesh.property(vertMatrixProp, vertex1) + heMesh.property(vertMatrixProp, vertex2);
			//cout << edgeMatrix << endl;
			Eigen::Matrix4d mat;
			mat << edgeMatrix(0, 0), edgeMatrix(0, 1), edgeMatrix(0, 2), edgeMatrix(0, 3),
			       edgeMatrix(0, 1), edgeMatrix(1, 1), edgeMatrix(1, 2), edgeMatrix(1, 3),
			       edgeMatrix(0, 2), edgeMatrix(1, 2), edgeMatrix(2, 2), edgeMatrix(2, 3),
			       0, 0, 0, 1;
			Eigen::Vector4d b(0, 0, 0, 1);
			Eigen::Vector4d x;
			x            = mat.colPivHouseholderQr().solve(b);
			double error = x.transpose() * edgeMatrix * x;
			if (heMesh.is_boundary(halfedge1) || heMesh.is_boundary(halfedge2))
			{
				//for edge that is boundarys, increase its error to a large value
				error += 10000.0;
			}
			heMesh.property(edgeMatrixProp, edge) = make_pair(error, x);
			heap.push(make_pair(error, edge));
		}

		
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
		QEMMesh::HalfedgeHandle heh2 = heMesh.halfedge_handle(targetEdge.second, 1);
		/*heMesh.delete_edge(targetEdge.second,false);
		heMesh.delete_vertex(heMesh.to_vertex_handle(heh1),false);
		heMesh.delete_vertex(heMesh.from_vertex_handle(heh1),false);*/
		//cout << "x : " << heMesh.point(heMesh.to_vertex_handle(heh1))[0] << " y: " << heMesh.point(heMesh.to_vertex_handle(heh1))[1] << " z : " << heMesh.point(heMesh.to_vertex_handle(heh1))[2] << endl;
		//cout << "x : " << heMesh.point(heMesh.from_vertex_handle(heh1))[0] << " y: " << heMesh.point(heMesh.from_vertex_handle(heh1))[1] << " z : " << heMesh.point(heMesh.from_vertex_handle(heh1))[2] << endl;
		//for (QEMMesh::VertexIter v_it = heMesh.vertices_begin(); v_it != heMesh.vertices_end(); ++v_it)
		//{
		//	cout << "vertex : " << heMesh.status(*v_it).deleted() << endl;
		//}
		//std::cout << "size : " << heap.size() << endl;
		//int cnt = 0;
		//for (QEMMesh::HalfedgeIter h_it = heMesh.halfedges_begin(); h_it != heMesh.halfedges_end(); ++h_it)
		//{
		//	cout << "half edge : " << heMesh.status(*h_it).deleted() << endl;
		//	++cnt;
		//}
		//cout << "cnt:" << cnt << endl;
		//for (QEMMesh::VertexIter bv_it = heMesh.vertices_begin(); v_it != heMesh.vertices_end(); ++v_it)
		//{
		//	cout << "vertex : " << heMesh.status(*v_it).deleted() << endl;
		//}
		//for (QEMMesh::EdgeIter e_it = heMesh.edges_begin(); e_it != heMesh.edges_end(); ++e_it)
		//{
		//	cout << "edge : " << heMesh.status(*e_it).deleted() << endl;
		//}
		//for (QEMMesh::FaceIter f_it = heMesh.faces_begin(); f_it != heMesh.faces_end(); ++f_it)
		//{
		//	cout << "face : " << heMesh.status(*f_it).deleted() << endl;
		//}
		heMesh.collapse(heh1);
		//QEMMesh::VertexHandle newVertex = heMesh.to_vertex_handle(heh1);
		//Eigen::Vector3d       pos       = heMesh.property(edgeMatrixProp, heMesh.edge_handle(heh1)).second.segment<3>(0);
		//QEMMesh::Point        p(pos.x(), pos.y(), pos.z());
		//heMesh.set_point(newVertex, p);
		
		bool flag = true;
		cout << "to 1 : " << heMesh.status(heMesh.to_vertex_handle(heh1)).deleted() << endl;
		cout << "from 1 : " << heMesh.status(heMesh.from_vertex_handle(heh1)).deleted()<<endl;
		cout << "to 2 : " << heMesh.status(heMesh.to_vertex_handle(heh2)).deleted() << endl;
		cout << "from 2 : " << heMesh.status(heMesh.from_vertex_handle(heh2)).deleted() << endl;
		cout << "half edge 1: " << heMesh.status(heh1).deleted() << endl;
		cout << "half edge 2: " << heMesh.status(heh2).deleted() << endl;
		cout << "edge : " << heMesh.status(targetEdge.second).deleted() << endl;
		cout << "flag : " << flag << endl;
		heMesh.garbage_collection();
		//bool valid = heMesh.is_valid_handle(newVertex);
		//if (valid)
		//{
		//	updataVertexNormal(newVertex);
		//	Eigen::Matrix4d mat;
		//	mat.setZero();
		//	for (QEMMesh::VertexFaceIter vf_it = heMesh.vf_iter(newVertex); vf_it.is_valid(); ++vf_it)
		//	{
		//		mat += countFaceKPMatrix(*vf_it);
		//	}
		//	heMesh.property(vertMatrixProp, newVertex) = mat;


		//	cout << "valid" << endl;
		//}

		



		
		//while (!heap.empty() && heMesh.n_vertices() > targetFaceCount)
		//{
		//	
		//	// 新顶点
		//	
		//	break;


		//}

		// write mesh to output.obj
		try
		{
			if (!OpenMesh::IO::write_mesh(heMesh, "output.obj"))
			{
				std::cerr << "Cannot write mesh to file 'output.obj'" << std::endl;
				return 1;
			}
		}
		catch (std::exception &x)
		{
			std::cerr << x.what() << std::endl;
			return 1;
		}

		return true;
	}
};