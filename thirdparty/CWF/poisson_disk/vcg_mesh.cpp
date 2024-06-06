#include "vcg_mesh.hpp"

#include <spdlog/spdlog.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/append.h>

namespace vcg {

MyMesh::MyMesh() { bbox = Box3<ScalarType>(); }

MyMesh::MyMesh(const MyMesh &m) {
  Clear();
  // Since append is not using const iterators for m I'm forced to use
  // the ugly const_cast(). shame I know ...
  tri::Append<MyMesh, MyMesh>::Mesh(*this, const_cast<MyMesh &>(m));
}

void MyMesh::concat(const float *verts, const int *triangles, int n_vert, int n_face) {
  vcg::tri::Allocator<MyMesh>::AddVertices(*this, n_vert);
  vcg::tri::Allocator<MyMesh>::AddFaces(*this, n_face);

  std::vector<VertexPointer> ivp(n_vert);
  // 设置顶点
  VertexIterator vi = vert.end();
  for (int i = (n_vert - 1); i >= 0; --i) {
    --vi;
    ivp[i] = &*vi;
    (*vi).P() = CoordType(verts[i * 3 + 0], verts[i * 3 + 1], verts[i * 3 + 2]);
  }

  // 设置面
  FaceIterator fi = face.end();
  for (int i = (n_face - 1); i >= 0; --i) {
    --fi;
    (*fi).V(0) = ivp[triangles[i * 3 + 0]];
    (*fi).V(1) = ivp[triangles[i * 3 + 1]];
    (*fi).V(2) = ivp[triangles[i * 3 + 2]];
  }
}

void MyMesh::concat(const Eigen::MatrixXf &V, const Eigen::MatrixXi &F) {
  int n_vert = V.rows();
  int n_face = F.rows();
  vcg::tri::Allocator<MyMesh>::AddVertices(*this, n_vert);
  vcg::tri::Allocator<MyMesh>::AddFaces(*this, n_face);

  std::vector<VertexPointer> ivp(n_vert);
  // 设置顶点
  VertexIterator vi = vert.end();
  for (int i = (n_vert - 1); i >= 0; --i) {
    --vi;
    ivp[i] = &*vi;
    (*vi).P() = CoordType(V(i, 0), V(i, 1), V(i, 2));
  }

  // 设置面
  FaceIterator fi = face.end();
  for (int i = (n_face - 1); i >= 0; --i) {
    --fi;
    (*fi).V(0) = ivp[F(i, 0)];
    (*fi).V(1) = ivp[F(i, 1)];
    (*fi).V(2) = ivp[F(i, 2)];
  }
}

void MyMesh::set_normals(const float *normals) {
  vcg::MyMesh::VertexIterator vi = vert.begin();
  for (int i = 0; i < vert.size(); ++i, ++vi) {
    (*vi).N() = NormalType(normals[i * 3 + 0], normals[i * 3 + 1], normals[i * 3 + 2]);
    //  spdlog::info("{}, {}, {}", normals[i * 3 + 0], normals[i * 3 + 1], normals[i * 3 + 2]);
  }
}

void MyMesh::set_normals(const Eigen::MatrixXf &N) {
  vcg::MyMesh::VertexIterator vi = vert.begin();
  for (int i = 0; i < vert.size(); ++i, ++vi) {
    (*vi).N() = NormalType(N(i, 0), N(i, 1), N(i, 2));
  }
}

void MyMesh::update_bb() { vcg::tri::UpdateBounding<MyMesh>::Box(*this); }

void MyMesh::fill_mesh(const float *verts, const int *triangles, int n_vert, int n_face) {
  Clear();
  concat(verts, triangles, n_vert, n_face);
}

void MyMesh::fill_mesh(const Eigen::MatrixXf &V, const Eigen::MatrixXi &F) {
  Clear();
  concat(V, F);
}

namespace MyAlgorithms {
void poison_disk_sampling(MyMesh &vcg_mesh, const PoisonParams &params, MyMesh &samples) {
  // 参数处理
  float rad = std::max(params._radius, 0.0f);
  int n_samp = std::max(params._n_samples, 1);
  int mc_rate = params._montecarlo_rate;

  bool sub_sample = params._sub_sampling;
  bool approx_geodesic = params._approx_geodesic_dist;
  bool refine = params._refine;

  if (rad <= 0) {
    // 根据样本数计算半径
    rad = SurfaceSampling::ComputePoissonDiskRadius(vcg_mesh, n_samp);
  } else {
    // 根据半径计算样本数
    n_samp = SurfaceSampling::ComputePoissonSampleNum(vcg_mesh, rad);
  }
  spdlog::info("计算泊松采样: r = {}, n = {}", rad, n_samp);

  // 首先生成 蒙特卡罗采样 用于快速查找
  MyMesh *presampledMesh;

  // this mesh is used only if we need real poisson sampling (and therefore we need to choose points
  // different from the starting mesh vertices)
  // 仅当我们需要真实的泊松采样时使用(因此我们需要选择与起始网格顶点不同的点)
  MyMesh montecarloMesh;

  if (sub_sample) {
    // 降采样
    presampledMesh = &vcg_mesh;
  } else {
    // 蒙特卡罗采样
    presampledMesh = &montecarloMesh;
    MySampler sampler(presampledMesh);

    SurfaceSampling::Montecarlo(vcg_mesh, sampler, n_samp * mc_rate);

    montecarloMesh.bbox = vcg_mesh.bbox; // 使用相同的bounding box

    spdlog::info("生成 {} 蒙特卡洛采样", montecarloMesh.vn);
  }

  SurfaceSampling::PoissonDiskParam pp;

  if (refine) {
    pp.preGenFlag = true;
    pp.preGenMesh = &vcg_mesh;
  }

  pp.geodesicDistanceFlag = approx_geodesic;
  pp.bestSampleChoiceFlag = true;
  pp.bestSamplePoolSize = 10;

  {
    auto timer = Utils::ScopeTimer("计算泊松采样");
    MySampler sampler(&samples);
    // SurfaceSampling::PoissonDiskPruningByNumber(sampler, *presampledMesh, n_samp, rad, pp);
    SurfaceSampling::PoissonDiskPruning(sampler, *presampledMesh, rad, pp);
    vcg_mesh.update_bb();
    // auto &g = pp.pds.gridSize;
    // spdlog::info("网格大小 {} x {} x {}", g[0], g[1], g[2]);

    spdlog::info("生成泊松采样: 顶点数 {}", samples.vn);
  }
}
} // namespace MyAlgorithms

} // namespace vcg
