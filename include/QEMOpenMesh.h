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
		const auto &n = heMesh.normal(fh);

		const double            a        = n[0];
		const double            b        = n[1];
		const double            c        = n[2];
		
		QEMMesh::FaceVertexIter fv_it  = heMesh.fv_iter(fh);
		QEMMesh::VertexHandle   vert   = *fv_it;
		QEMMesh::Point          position = heMesh.point(vert);
		const double            d        = -(heMesh.point(vert)| n);
		Eigen::Vector4d         p(a,b,c,d);
		Eigen::Matrix4d         mat = p * p.transpose();
		return mat;
	}

	Eigen::Matrix4d updateVertexMatrix(QEMMesh::VertexHandle vertex)
	{
		Eigen::Matrix4d mat;
		mat.setZero();
		for (QEMMesh::VertexFaceIter vf_it = heMesh.vf_iter(vertex); vf_it.is_valid(); ++vf_it)
		{
			if (heMesh.status(*vf_it).deleted())
				continue; 
			mat += countFaceKPMatrix(*vf_it);
		}
		return mat;
	}
	//检车矩阵是否可逆
	bool checkInvertible(Eigen::Matrix4d mat)
	{
		Eigen::FullPivLU<Eigen::Matrix4d> lu(mat);
		if (lu.isInvertible())
			return true;
		else
			return false;
	}

	bool isFlipped(QEMMesh::VertexHandle v0, QEMMesh::VertexHandle v1, const QEMMesh::Point &ptTarget)
	{
		for (auto fIt = heMesh.vf_iter(v0); fIt.is_valid(); ++fIt)
		{
			if (heMesh.status(*fIt).deleted())
			{
				continue;
			}

			auto                   vIt = heMesh.fv_iter(*fIt);
			OpenMesh::VertexHandle fv[3];
			fv[0] = *vIt;
			++vIt;
			fv[1] = *vIt;
			++vIt;
			fv[2] = *vIt;

			// ignore the face will be deleted
			if (fv[0] == v1 || fv[1] == v1 || fv[2] == v1)
			{
				continue;
			}

			int idxV0 = 0;
			for (int i = 0; i < 3; ++i)
			{
				if (fv[i] == v0)
				{
					idxV0 = i;
				}
			}

			QEMMesh::Point pt0 = heMesh.point(v0);
			QEMMesh::Point pt1 = heMesh.point(fv[(idxV0 + 1) % 3]);
			QEMMesh::Point pt2 = heMesh.point(fv[(idxV0 + 2) % 3]);

			QEMMesh::Point dir1 = pt1 - ptTarget;
			dir1.normalize();
			QEMMesh::Point dir2 = pt2 - ptTarget;
			dir2.normalize();

			// The angle below can be adjusted, but now is enough
			// if the angle between dir1 and dir2 small than 2.6 angle, return true
			if (fabs(dir1 | dir2) > 0.999)
				return true;

			QEMMesh::Point normold = heMesh.normal(*fIt);
			QEMMesh::Point norm    = dir1 % (dir2);
			norm.normalize();

			// if the angle between normold and norm large than 78.5 angle, return true
			if ((normold | norm) < 0.2)
				return true;
		}
		return false;
	}
	void UpdateFaceNormal(OpenMesh::VertexHandle v0)
	{
		if (!heMesh.is_valid_handle(v0) || heMesh.is_isolated(v0))
		{
			return;
		}
		for (auto fIt = heMesh.vf_iter(v0); fIt.is_valid(); ++fIt)
		{
			if (!heMesh.status(*fIt).deleted())
			{
				heMesh.set_normal(*fIt, heMesh.calc_face_normal(*fIt));
			}
		}
	}
	
	inline bool DoSimplification(SimplificationMode mode, float ratio)
	{
		/*------------------init-----------------*/
		QEM_DEBUG("DoSimplification(mode=%d, targetFaceCount=%d)", mode, targetFaceCount);
		cout << "hello qem" << endl;
		double dAgressiveness = 7;
		if (!heMesh.has_vertex_status())
			heMesh.request_vertex_status();
		if (!heMesh.has_face_status())
			heMesh.request_face_status();
		if (!heMesh.has_edge_status())
			heMesh.request_edge_status();
		if (!heMesh.has_face_normals())
		{
			heMesh.request_face_normals();
		}
		heMesh.update_face_normals();
		OpenMesh::VPropHandleT<Eigen::Matrix4d> vertMatrixProp;        
		heMesh.add_property(vertMatrixProp);
		//遍历顶点求kp矩阵
		for (auto& vh : heMesh.vertices())
		{
			heMesh.property(vertMatrixProp, vh) = updateVertexMatrix(vh);
		}

		// error new-vertex-position
		OpenMesh::EPropHandleT<std::pair<double,Eigen::Vector4d>> edgeMatrixProp;
		heMesh.add_property(edgeMatrixProp);
		std::priority_queue<std::pair<double, QEMMesh::EdgeHandle>, vector<std::pair<double, QEMMesh::EdgeHandle>>, greater<std::pair<double, QEMMesh::EdgeHandle>>> heap;
		// 遍历每一条边
		for (QEMMesh::EdgeIter e_it = heMesh.edges_sbegin(); e_it != heMesh.edges_end(); ++e_it)
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
			Eigen::Vector4d x(0,0,0,0);
			double          error = 10000;
			if (checkInvertible(mat) && !heMesh.is_boundary(edge))
			{
				Eigen::Vector4d b(0, 0, 0, 1);
				x     = mat.inverse() * b;
				error = x.transpose() * edgeMatrix * x;
			}
			else
			{
				//如果矩阵不可逆，则取边的两个端点和中点中误差最小的点
				Eigen::Vector4d v1(heMesh.point(vertex1)[0], heMesh.point(vertex1)[1], heMesh.point(vertex1)[2], 1.0);
				Eigen::Vector4d v2(heMesh.point(vertex2)[0], heMesh.point(vertex2)[1], heMesh.point(vertex2)[2], 1.0);
				Eigen::Vector4d vmid = (v1 + v2) / 2;
				double          e1   = v1.transpose() * edgeMatrix * v1;
				double          e2   = v2.transpose() * edgeMatrix * v2;
				double          e3   = vmid.transpose() * edgeMatrix * vmid;
				if (e1 < error)
				{
					error = e1;
					x     = v1;
				}
				if (e2 < error)
				{
					error = e2;
					x     = v2;
				}
				if (e3 < error)
				{
					error = e3;
					x     = vmid;
				}
			}
			heMesh.property(edgeMatrixProp, edge) = make_pair(error, x);
			heap.push(make_pair(error, edge));
		}
		/*------------------init over-----------------*/

		/*----------------- version 2 -------------------*/
		// 需要简化的次数
		int cnt = heMesh.n_vertices()*(1-ratio);
		int idx = 0;
		while (!heap.empty() && heMesh.n_vertices() > cnt && idx < cnt)
		{
			auto targetEdge = heap.top();
			heap.pop();
			//有效性检测
			if (heMesh.status(targetEdge.second).deleted())
				continue;
			if (!heMesh.is_valid_handle(targetEdge.second))
				continue;

			QEMMesh::HalfedgeHandle heh1 = heMesh.halfedge_handle(targetEdge.second, 0);
			OpenMesh::VertexHandle   v0 = heMesh.from_vertex_handle(heh1);
			OpenMesh::VertexHandle   v1 = heMesh.to_vertex_handle(heh1);

			if (heMesh.is_boundary(v0) != heMesh.is_boundary(v1))
			{
				continue;
			}
			Eigen::Vector4d x = heMesh.property(edgeMatrixProp, targetEdge.second).second;
			QEMMesh::Point  ptTarget(x[0], x[1], x[2]);
			if (isFlipped(v0, v1, ptTarget))
				continue;
			if (isFlipped(v1, v0, ptTarget))
				continue;

			auto heh2 = heMesh.opposite_halfedge_handle(heh1);
			if (heMesh.is_collapse_ok(heh1))
			{
				heMesh.collapse(heh1);
				heMesh.point(v1) = ptTarget;
				++idx;
			}
			else if (heMesh.is_collapse_ok(heh2))
			{
				heMesh.collapse(heh2);
				heMesh.point(v0) = ptTarget;
				++idx;
			}
			else
			{
				continue;
			}

			//更新新顶点周围的元素
			UpdateFaceNormal(v0);
			UpdateFaceNormal(v1);
			auto &v0Quadric = heMesh.property(vertMatrixProp, v0);
			auto &v1Quadric = heMesh.property(vertMatrixProp, v1);
			heMesh.property(vertMatrixProp, v1) += v0Quadric;
			heMesh.property(vertMatrixProp, v0) += v1Quadric;
			// update edge matrix
			if (heMesh.is_valid_handle(v0) || !heMesh.is_isolated(v0))
			{
			    double         dError = 0;
			    QEMMesh::Point ptResult;
			    for (auto hIt = heMesh.voh_iter(v0); hIt.is_valid(); ++hIt)
			    {
			        OpenMesh::EdgeHandle edge = heMesh.edge_handle(*hIt);

			        if (!heMesh.is_valid_handle(*hIt) || heMesh.status(edge).deleted())
			        {
			        	continue;
			        }
			        // 获取边的两个顶点
			        QEMMesh::HalfedgeHandle halfedge1 = heMesh.halfedge_handle(edge, 0);
			        QEMMesh::HalfedgeHandle halfedge2 = heMesh.halfedge_handle(edge, 1);

			        QEMMesh::VertexHandle vertex1 = heMesh.to_vertex_handle(halfedge1);
			        QEMMesh::VertexHandle vertex2 = heMesh.from_vertex_handle(halfedge1);

			        // 在这里使用 vertex1 和 vertex2，它们分别是边的两个顶点
			        Eigen::Matrix4d edgeMatrix = heMesh.property(vertMatrixProp, vertex1) + heMesh.property(vertMatrixProp, vertex2);
			        Eigen::Matrix4d mat;
			        mat << edgeMatrix(0, 0), edgeMatrix(0, 1), edgeMatrix(0, 2), edgeMatrix(0, 3),
			        	edgeMatrix(0, 1), edgeMatrix(1, 1), edgeMatrix(1, 2), edgeMatrix(1, 3),
			        	edgeMatrix(0, 2), edgeMatrix(1, 2), edgeMatrix(2, 2), edgeMatrix(2, 3),
			        	0, 0, 0, 1;
			        Eigen::Vector4d x(0, 0, 0, 0);
			        double          error = 10000;
			        if (checkInvertible(mat) && !heMesh.is_boundary(edge))
			        {
			        	Eigen::Vector4d b(0, 0, 0, 1);
			        	x     = mat.inverse() * b;
			        	error = x.transpose() * edgeMatrix * x;
			        }
			        else
			        {
			        	// 如果矩阵不可逆，则取边的两个端点和中点中误差最小的点
			        	Eigen::Vector4d v1(heMesh.point(vertex1)[0], heMesh.point(vertex1)[1], heMesh.point(vertex1)[2], 1.0);
			        	Eigen::Vector4d v2(heMesh.point(vertex2)[0], heMesh.point(vertex2)[1], heMesh.point(vertex2)[2], 1.0);
			        	Eigen::Vector4d vmid = (v1 + v2) / 2;
			        	double          e1   = v1.transpose() * edgeMatrix * v1;
			        	double          e2   = v2.transpose() * edgeMatrix * v2;
			        	double          e3   = vmid.transpose() * edgeMatrix * vmid;
			        	if (e1 < error)
			        	{
			        		error = e1;
			        		x     = v1;
			        	}
			        	if (e2 < error)
			        	{
			        		error = e2;
			        		x     = v2;
			        	}
			        	if (e3 < error)
			        	{
			        		error = e3;
			        		x     = vmid;
			        	}
			        }
			        heMesh.property(edgeMatrixProp, edge) = make_pair(error, x);
			    }
			}
			if (heMesh.is_valid_handle(v1) || !heMesh.is_isolated(v1))
			{
			    double         dError = 0;
			    QEMMesh::Point ptResult;
			    for (auto hIt = heMesh.voh_iter(v1); hIt.is_valid(); ++hIt)
			    {
			        OpenMesh::EdgeHandle edge = heMesh.edge_handle(*hIt);

			        if (!heMesh.is_valid_handle(*hIt) || heMesh.status(edge).deleted())
			        {
			        	continue;
			        }
			        // 获取边的两个顶点
			        QEMMesh::HalfedgeHandle halfedge1 = heMesh.halfedge_handle(edge, 0);
			        QEMMesh::HalfedgeHandle halfedge2 = heMesh.halfedge_handle(edge, 1);

			        QEMMesh::VertexHandle vertex1 = heMesh.to_vertex_handle(halfedge1);
			        QEMMesh::VertexHandle vertex2 = heMesh.from_vertex_handle(halfedge1);

			        // 在这里使用 vertex1 和 vertex2，它们分别是边的两个顶点
			        Eigen::Matrix4d edgeMatrix = heMesh.property(vertMatrixProp, vertex1) + heMesh.property(vertMatrixProp, vertex2);
			        Eigen::Matrix4d mat;
			        mat << edgeMatrix(0, 0), edgeMatrix(0, 1), edgeMatrix(0, 2), edgeMatrix(0, 3),
			        	edgeMatrix(0, 1), edgeMatrix(1, 1), edgeMatrix(1, 2), edgeMatrix(1, 3),
			        	edgeMatrix(0, 2), edgeMatrix(1, 2), edgeMatrix(2, 2), edgeMatrix(2, 3),
			        	0, 0, 0, 1;
			        Eigen::Vector4d x(0, 0, 0, 0);
			        double          error = 10000;
			        if (checkInvertible(mat) && !heMesh.is_boundary(edge))
			        {
			        	Eigen::Vector4d b(0, 0, 0, 1);
			        	x     = mat.inverse() * b;
			        	error = x.transpose() * edgeMatrix * x;
			        }
			        else
			        {
			        
			        	// 如果矩阵不可逆，则取边的两个端点和中点中误差最小的点
			        	Eigen::Vector4d v1(heMesh.point(vertex1)[0], heMesh.point(vertex1)[1], heMesh.point(vertex1)[2], 1.0);
			        	Eigen::Vector4d v2(heMesh.point(vertex2)[0], heMesh.point(vertex2)[1], heMesh.point(vertex2)[2], 1.0);
			        	Eigen::Vector4d vmid = (v1 + v2) / 2;
			        	double          e1   = v1.transpose() * edgeMatrix * v1;
			        	double          e2   = v2.transpose() * edgeMatrix * v2;
			        	double          e3   = vmid.transpose() * edgeMatrix * vmid;
			        	if (e1 < error)
			        	{
			        		error = e1;
			        		x     = v1;
			        	}
			        	if (e2 < error)
			        	{
			        		error = e2;
			        		x     = v2;
			        	}
			        	if (e3 < error)
			        	{
			        		error = e3;
			        		x     = vmid;
			        	}
			        }
			        heMesh.property(edgeMatrixProp, edge) = make_pair(error, x);
			    }
			}
		}

		heMesh.garbage_collection();
		if (heMesh.has_vertex_status())
			heMesh.release_vertex_status();
		if (heMesh.has_face_status())
			heMesh.release_face_status();
		if (heMesh.has_edge_status())
			heMesh.release_edge_status();
		cout << "face count : " << heMesh.n_faces() << endl;
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