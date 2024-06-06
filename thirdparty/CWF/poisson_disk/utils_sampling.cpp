#include "utils_sampling.hpp"
#include "vcg_mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/write_xyz_points.h>

#include <spdlog/spdlog.h>
#include <utility>

namespace UtilsSampling {
void poisson_disk(float radius,
                  int n_samples,
                  const std::vector<Vec3> &verts,
                  const std::vector<Vec3> &nors,
                  const std::vector<int> &tris,
                  std::vector<Vec3> &samples_pos,
                  std::vector<Vec3> &samples_nor) {
  assert(verts.size() > 0);
  assert(tris.size() > 0);

  vcg::MyMesh vcg_mesh, sampler;
  vcg_mesh.fill_mesh((float *)verts.data(), (int *)tris.data(), verts.size(), tris.size() / 3);
  // 当法线信息正确时
  if (verts.size() == nors.size()) {
    vcg_mesh.set_normals((float *)nors.data());
  }
  vcg_mesh.update_bb();

  vcg::MyAlgorithms::PoisonParams pp;
  pp._radius = radius;
  pp._n_samples = n_samples;
  pp._approx_geodesic_dist = true;
  vcg::MyAlgorithms::poison_disk_sampling(vcg_mesh, pp, sampler);

  const int n_vert = sampler.vert.size();
  samples_pos.clear();
  samples_nor.clear();
  samples_pos.resize(n_vert);
  samples_nor.resize(n_vert);
  vcg::MyMesh::VertexIterator vi = sampler.vert.begin();
  for (int i = 0; i < n_vert; ++i, ++vi) {
    vcg::MyMesh::CoordType p = (*vi).P();
    vcg::MyMesh::NormalType n = (*vi).N();
    samples_pos[i] = {p.X(), p.Y(), p.Z()};
    samples_nor[i] = {n.X(), n.Y(), n.Z()};
  }
}

void poisson_disk(float radius, int n_samples, const Eigen::MatrixXf &V, const Eigen::MatrixXf &N, const Eigen::MatrixXi &F, std::vector<Vec3> &samples_pos, std::vector<Vec3> &samples_nor) {
  assert(V.rows() == N.rows());

  vcg::MyMesh vcg_mesh, sampler;

  vcg_mesh.fill_mesh(V, F);
  // 当法线信息正确时
  if (V.rows() == N.rows()) {
    vcg_mesh.set_normals(N);
  }
  vcg_mesh.update_bb();

  vcg::MyAlgorithms::PoisonParams pp;
  pp._radius = radius;
  pp._n_samples = n_samples;
  pp._approx_geodesic_dist = false;
  vcg::MyAlgorithms::poison_disk_sampling(vcg_mesh, pp, sampler);

  const int n_vert = sampler.vert.size();
  samples_pos.clear();
  samples_nor.clear();
  samples_pos.resize(n_vert);
  samples_nor.resize(n_vert);
  vcg::MyMesh::VertexIterator vi = sampler.vert.begin();
  for (int i = 0; i < n_vert; ++i, ++vi) {
    vcg::MyMesh::CoordType p = (*vi).P();
    vcg::MyMesh::NormalType n = (*vi).N();
    samples_pos[i] = {p.X(), p.Y(), p.Z()};
    samples_nor[i] = {n.X(), n.Y(), n.Z()};
  }
}

void poisson_disk_to_file(float radius, int &n_samples, const Eigen::MatrixXf &V, const Eigen::MatrixXf &N, const Eigen::MatrixXi &F, std::function<std::string(int)> get_name) {
  assert(V.rows() == N.rows());

  vcg::MyMesh vcg_mesh, sampler;

  vcg_mesh.fill_mesh(V, F);
  // 当法线信息正确时
  if (V.rows() == N.rows()) {
    vcg_mesh.set_normals(N);
  }
  vcg_mesh.update_bb();

  vcg::MyAlgorithms::PoisonParams pp;
  pp._radius = radius;
  pp._n_samples = n_samples;
  pp._approx_geodesic_dist = false;
  vcg::MyAlgorithms::poison_disk_sampling(vcg_mesh, pp, sampler);
  n_samples = sampler.vert.size();
  
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point = Kernel::Point_3;
  using Vector = Kernel::Vector_3;
  using Pwn = std::pair<Point, Vector>;

  std::vector<Pwn> points;
  points.reserve(sampler.vert.size());

  for (auto &v : sampler.vert) {
    auto &p = v.P();
    auto &n = v.N();
    points.emplace_back(Point(p.X(), p.Y(), p.Z()), Vector(n.X(), n.Y(), n.Z()));
  }

  // 描述xyz写入时如何处理数据, point从哪取, normal从哪取, 精度是多少
  auto params = CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()).normal_map(CGAL::Second_of_pair_property_map<Pwn>()).stream_precision(6);
  CGAL::IO::write_XYZ(get_name(n_samples), points, params);
}
} // namespace UtilsSampling
