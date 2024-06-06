#include "Utils/Common.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <spdlog/formatter.h>
#include <time.h>

#include <BGAL/Algorithm/BOC/BOC.h>
#include <BGAL/CVTLike/CVT.h>
#include <BGAL/Integral/Integral.h>
#include <BGAL/Optimization/LinearSystem/LinearSystem.h>

#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/readOFF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>

using Path = std::filesystem::path;

typedef CGAL::Simple_cartesian<double> K_T;
typedef K_T::FT FT;
typedef K_T::Point_3 Point_T;

typedef K_T::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K_T> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K_T, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
double gamma = 0.00000000000001;

struct MyPoint {
  MyPoint(Eigen::Vector3d a) { p = a; }

  MyPoint(double a, double b, double c) {
    p.x() = a;
    p.y() = b;
    p.z() = c;
  }
  Eigen::Vector3d p;

  bool operator<(const MyPoint &a) const {

    double dis = (p - a.p).norm();
    if (dis < gamma) {
      return false;
    }

    if ((p.x() - a.p.x()) < 0.00000000001 && (p.x() - a.p.x()) > -0.00000000001) {
      if ((p.y() - a.p.y()) < 0.00000000001 && (p.y() - a.p.y()) > -0.00000000001) {
        return (p.z() < a.p.z());
      }
      return (p.y() < a.p.y());
    }
    return (p.x() < a.p.x());
  }

  bool operator==(const MyPoint &a) const {
    if ((p.x() - a.p.x()) < 0.00000000001 && (p.x() - a.p.x()) > -0.00000000001) {
      if ((p.y() - a.p.y()) < 0.00000000001 && (p.y() - a.p.y()) > -0.00000000001) {
        if ((p.z() - a.p.z()) < 0.00000000001 && (p.z() - a.p.z()) > -0.00000000001) {
          return 1;
        }
      }
    }
    return 0;
  }
};

struct MyFace {
  MyFace(Eigen::Vector3i a) { p = a; }

  MyFace(int a, int b, int c) {
    p.x() = a;
    p.y() = b;
    p.z() = c;
  }

  Eigen::Vector3i p;
  bool operator<(const MyFace &a) const {
    if (p.x() == a.p.x()) {
      if (p.y() == a.p.y()) {
        return p.z() > a.p.z();
      }
      return p.y() > a.p.y();
    }
    return p.x() > a.p.x();
  }
};

namespace BGAL {
_CVT3D::_CVT3D(const _ManifoldModel &model) : _model(model), _RVD(model), _RVD2(model), _para() {
  _rho = [](BGAL::_Point3 &p) { return 1; };
  _para.is_show = true;
  _para.epsilon = 1e-30;
  _para.max_linearsearch = 20;
}

_CVT3D::_CVT3D(const _ManifoldModel &model, std::function<double(_Point3 &p)> &rho, _LBFGS::_Parameter para) : _model(model), _RVD(model), _RVD2(model), _rho(rho), _para(para) {}

// 输出计算结果
void _CVT3D::DoOutput(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step) {
  if (_remesh_params.bOutputXyz)
    OutputXYZ(sites, RVD, num, modelname, step);
  if (_remesh_params.bOutputRemesh)
    OutputRemesh(sites, RVD, num, modelname, step);
  if (_remesh_params.bOutputRVD)
    OutputRVD(sites, RVD, num, modelname, step);
}

// --------------------------------------------写入XYZ--------------------------------------------
void _CVT3D::OutputXYZ(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step) {
  auto outpath = Path(_remesh_params.workDir) / "_LBFGSOUT";
  auto xyzpath = outpath / fmt::format("Ours_{}_{}_Points.xyz", num, modelname);

  if (step > 2) {
    xyzpath = outpath / fmt::format("Ours_{}_{}_Iter{}_Points.xyz", num, modelname, step - 3);
  }

  std::ofstream outP(xyzpath);

  int outnum = step == 1 ? sites.size() / 3 : sites.size();

  for (int i = 0; i < outnum; ++i) {
    outP << sites[i] << std::endl;
  }

  outP.close();
}

// --------------------------------------------写入Remesh--------------------------------------------
void _CVT3D::OutputRemesh(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step) {
  if (step < 2)
    return;
  auto outpath = Path(_remesh_params.workDir) / "_LBFGSOUT";
  auto remeshpath = outpath / fmt::format("Ours_{}_{}_Iter{}_Remesh.obj", num, modelname, to_string(step - 3));

  std::ofstream outRDT(remeshpath);

  auto Vs = sites;
  auto Edges = RVD.get_edges_();
  set<pair<int, int>> RDT_Edges;
  vector<set<int>> neibors;
  neibors.resize(Vs.size());
  for (int i = 0; i < Edges.size(); i++) {
    for (auto ee : Edges[i]) {
      RDT_Edges.insert({min(i, ee.first), max(i, ee.first)});
      neibors[i].insert(ee.first);
      neibors[ee.first].insert(i);
    }
  }

  for (auto v : Vs) {
    outRDT << "v " << v << endl;
  }

  // set<MyFace> rdtFaces;

  for (auto e : RDT_Edges) {
    for (int pid : neibors[e.first]) {
      if (RDT_Edges.find({min(pid, e.first), max(pid, e.first)}) == RDT_Edges.end())
        continue;
      if (RDT_Edges.find({min(pid, e.second), max(pid, e.second)}) == RDT_Edges.end())
        continue;

      int f1 = pid, f2 = e.first, f3 = e.second;
      auto maxf = max(f1, max(f2, f3));
      auto minf = min(f1, min(f2, f3));
      int mid;

      if (f1 != maxf && f1 != minf) {
        mid = f1;
      }

      if (f2 != maxf && f2 != minf) {
        mid = f2;
      }

      if (f3 != maxf && f3 != minf) {
        mid = f3;
      }
      outRDT << "f " << maxf + 1 << " " << mid + 1 << " " << minf + 1 << endl;
      // rdtFaces.insert(MyFace(maxf, mid, minf));
    }
  }

  // for (auto f : rdtFaces) {
  //   outRDT << "f " << f.p.x() + 1 << " " << f.p.y() + 1 << " " << f.p.z() + 1 << endl;
  // }

  outRDT.close();
  if (step == 2)
    _remesh_params.finishedCallback(step, remeshpath.string());
  else
    _remesh_params.batchRemeshCallback(step, remeshpath.string());
}

// --------------------------------------------写入RVD--------------------------------------------
void _CVT3D::OutputRVD(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step) {
  const auto &cells = RVD.get_cells_();
  auto lbfgoutpath = Path(_remesh_params.workDir) / "_LBFGSOUT";
  auto rvdpath = lbfgoutpath / fmt::format("Ours_{}_{}_RVD.obj", num, modelname); // step <= 2

  if (step > 2) {
    rvdpath = lbfgoutpath / fmt::format("Ours_{}_{}_Iter{}_RVD.obj", num, modelname, to_string(step - 3));
  }

  std::ofstream out(rvdpath);
  out << "g 3D_Object\nmtllib BKLineColorBar.mtl\nusemtl BKLineColorBar" << std::endl;

  for (int i = 0; i < RVD.number_vertices_(); ++i) {
    out << "v " << RVD.vertex_(i) << std::endl;
  }

  for (int i = 0; i < cells.size(); ++i) {
    // auto &cell = cells[i];

    // 面积计算
    // for (int j = 0; j < cell.size(); ++j) {
    //   BGAL::_Point3 p1 = RVD.vertex_(cell[j][0]);
    //   BGAL::_Point3 p2 = RVD.vertex_(cell[j][1]);
    //   BGAL::_Point3 p3 = RVD.vertex_(cell[j][2]);
    //   area += (p2 - p1).cross_(p3 - p1).length_() / 2;
    // }

    // color计算
    // double color = BGAL::_BOC::rand_();
    // if (step == 1 && i > cells.size() / 3) {
    //   color = 0;
    // }
    // out << "vt " << color << " 0" << std::endl;

    // for (int j = 0; j < cell.size(); ++j) {
    //   out << "f " << cell[j][0] + 1 << "/" << i + 1 << " " << cell[j][1] + 1 << "/" << i + 1 << " " << cell[j][2] + 1 << "/" << i + 1 << std::endl;
    // }

    for (auto &p : cells[i]) {
      out << fmt::format("f {}/{} {}/{} {}/{}\n", p[0] + 1, i + 1, p[1] + 1, i + 1, p[2] + 1, i + 1);
    }
  }

  out.close();
}

void _CVT3D::loadpoints(string name) {
  if (loaded)
    return;
  auto xyzpath = _remesh_params.OutputDirPath() / name;
  ifstream inPoints(xyzpath);
  double x, y, z, nx, ny, nz; // if xyz file has normal

  while (inPoints >> x >> y >> z >> nx >> ny >> nz) {
    Pts.push_back(Eigen::Vector3d(x, y, z));
    Nors.push_back(Eigen::Vector3d(nx, ny, nz)); // Nors here is useless, if do not have normal, just set it to (1,0,0)
  }

  inPoints.close();
  loaded = true;
}

void _CVT3D::calculate_(int num_sites, char *modelName) {
  double allTime = 0, RVDtime = 0;
  clock_t start, end;
  clock_t startRVD, endRVD;

  double PI = 3.14159265358;
  string modelname = modelName;
  Polyhedron polyhedron;
  std::ifstream input(_remesh_params.OutputDirPath() / "Temp.off");
  input >> polyhedron;
  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);

  // string xyzpath = "../../../../data/n" + to_string(num_sites) + "_" + modelname +
  // "_inputPoints.xyz"; load(xyzpath);

  // begin step 1.
  int num = Pts.size();

  spdlog::debug("\nBegin CWF.");

  auto Fnum = _remesh_params.fnum;
  auto alpha = _remesh_params.alpha;
  auto eplison = _remesh_params.eplison;
  auto lambda = _remesh_params.lambda;
  auto decay = _remesh_params.decay;

  auto fgm2 = [&](const Eigen::VectorXd &X, Eigen::VectorXd &g) {
    eplison = eplison * decay;
    double lossCVT = 0, lossQE = 0, loss = 0;

    startRVD = clock();

    for (int i = 0; i < num; ++i) {
      Point_T query(X(i * 3 + 0), X(i * 3 + 1), X(i * 3 + 2)); // project to base surface
      Point_T closest = tree.closest_point(query);
      auto tri = tree.closest_point_and_primitive(query);

      Polyhedron::Face_handle f = tri.second;
      auto p1 = f->halfedge()->vertex()->point();
      auto p2 = f->halfedge()->next()->vertex()->point();
      auto p3 = f->halfedge()->next()->next()->vertex()->point();
      Eigen::Vector3d v1(p1.x(), p1.y(), p1.z());
      Eigen::Vector3d v2(p2.x(), p2.y(), p2.z());
      Eigen::Vector3d v3(p3.x(), p3.y(), p3.z());
      Eigen::Vector3d N = (v2 - v1).cross(v3 - v1);
      N.normalize();
      Nors[i] = N;
      BGAL::_Point3 p(closest.x(), closest.y(), closest.z());
      _sites[i] = p;
    }

    _RVD.calculate_(_sites);
    Fnum++;
    DoOutput(_sites, _RVD, num_sites, modelname, Fnum); // output process

    endRVD = clock();
    RVDtime += (double)(endRVD - startRVD) / CLOCKS_PER_SEC;

    const auto &cells = _RVD.get_cells_();
    const auto &edges = _RVD.get_edges_();
    double energy = 0.0;
    g.setZero();
    vector<Eigen::Vector3d> gi;
    gi.resize(num);

    for (int i = 0; i < num; ++i) {
      gi[i] = Eigen::Vector3d(0, 0, 0);
    }
    int thread_num;
    thread_num = omp_get_max_threads();
    spdlog::debug("max threads: {}", thread_num);
    omp_set_num_threads(thread_num); // change to your CPU core numbers
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
      for (int j = 0; j < cells[i].size(); ++j) {
        auto &cell = cells[i][j];
        Eigen::VectorXd inte = BGAL::_Integral::integral_triangle3D(
            [&](BGAL::_Point3 p) {
              Eigen::VectorXd r(5);

              BGAL::_Point3 NorTriM = (_RVD.vertex_(cell[1]) - _RVD.vertex_(cell[0])).cross_(_RVD.vertex_(cell[2]) - _RVD.vertex_(cell[0]));
              NorTriM.normalized_();

              BGAL::_Point3 Nori(Nors[i].x(), Nors[i].y(), Nors[i].z());

              r[0] = (eplison * _rho(p) * ((_sites[i] - p).sqlength_())); // CVT

              r[1] = lambda * (NorTriM.dot_(p - _sites[i])) * (NorTriM.dot_(p - _sites[i])) + eplison * ((p - _sites[i]).sqlength_()); // qe+CVT

              r[2] = lambda * -2 * NorTriM.x() * (NorTriM.dot_(p - _sites[i])) + eplison * -2 * (p - _sites[i]).x(); // g
              r[3] = lambda * -2 * NorTriM.y() * (NorTriM.dot_(p - _sites[i])) + eplison * -2 * (p - _sites[i]).y(); // g
              r[4] = lambda * -2 * NorTriM.z() * (NorTriM.dot_(p - _sites[i])) + eplison * -2 * (p - _sites[i]).z(); // g

              return r;
            },
            _RVD.vertex_(cell[0]), _RVD.vertex_(cell[1]), _RVD.vertex_(cell[2]));
        // energy += alpha * inte(1);
        lossCVT += alpha * inte[0];
        loss += alpha * inte[1];
        gi[i].x() += alpha * inte[2];
        gi[i].y() += alpha * inte[3];
        gi[i].z() += alpha * inte[4];
      }
    }

    for (int i = 0; i < num; i++) {
      gi[i] = gi[i] - Nors[i] * (gi[i].dot(Nors[i]) / Nors[i].dot(Nors[i]));
      g(i * 3 + 0) += gi[i].x();
      g(i * 3 + 1) += gi[i].y();
      g(i * 3 + 2) += gi[i].z();
    }
    energy += loss;

    spdlog::info("energy: {:.7f} LossCVT: {:.7f} LossQE: {:.7f} Lambda_CVT: {:.7f}", energy, lossCVT / eplison, loss - lossCVT, eplison);

    return energy;
  };

  vector<Eigen::Vector3d> Pts2;

  Pts2 = Pts;
  num = Pts2.size();
  spdlog::info("Points: {}", num);
  _sites.resize(num);
  _para.max_linearsearch = 20;
  _para.max_iteration = 50;
  BGAL::_LBFGS lbfgs2(_para);
  Eigen::VectorXd iterX2(num * 3);

  for (int i = 0; i < num; ++i) {
    iterX2(i * 3 + 0) = Pts2[i].x();
    iterX2(i * 3 + 1) = Pts2[i].y();
    iterX2(i * 3 + 2) = Pts2[i].z();
    _sites[i] = BGAL::_Point3(Pts2[i].x(), Pts2[i].y(), Pts2[i].z());
  }

  _RVD.calculate_(_sites);
  start = clock();
  lbfgs2.minimize(fgm2, iterX2);
  end = clock();
  allTime += (double)(end - start) / CLOCKS_PER_SEC;

  spdlog::info("allTime: {:.7f} RVDtime: {:.7f} L-BFGS time: {:.7f}", allTime, RVDtime, allTime - RVDtime);

  for (int i = 0; i < num; ++i) {
    // Point_T query(x0[i * 3], x0[i * 3+1], x0[i * 3+2]);
    Point_T query(iterX2(i * 3 + 0), iterX2(i * 3 + 1), iterX2(i * 3 + 2));
    Point_T closest = tree.closest_point(query);
    auto tri = tree.closest_point_and_primitive(query);

    Polyhedron::Face_handle f = tri.second;
    auto p1 = f->halfedge()->vertex()->point();
    auto p2 = f->halfedge()->next()->vertex()->point();
    auto p3 = f->halfedge()->next()->next()->vertex()->point();
    Eigen::Vector3d v1(p1.x(), p1.y(), p1.z());
    Eigen::Vector3d v2(p2.x(), p2.y(), p2.z());
    Eigen::Vector3d v3(p3.x(), p3.y(), p3.z());
    Eigen::Vector3d N = (v2 - v1).cross(v3 - v1);
    N.normalize();
    Nors[i] = N;

    _sites[i] = BGAL::_Point3(closest.x(), closest.y(), closest.z());
  }

  _RVD.calculate_(_sites);

  DoOutput(_sites, _RVD, num_sites, modelname, 2);
}

void _CVT3D::calculatet_(int num_sites, string modelname) {
  double allTime = 0, RVDtime = 0;
  clock_t start, end;
  clock_t startRVD, endRVD;

  double PI = 3.14159265358;
  Polyhedron polyhedron;
  auto offpath = _remesh_params.OutputDirPath() / "Temp.off";
  std::ifstream input(offpath);
  input >> polyhedron;
  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);

  // string xyzpath = "../../../../data/n" + to_string(num_sites) + "_" + modelname +
  // "_inputPoints.xyz"; load(xyzpath);

  // begin step 1.
  int num = Pts.size();

  spdlog::debug("\nBegin CWF.");
  auto Fnum = _remesh_params.fnum;
  auto alpha = _remesh_params.alpha;
  auto eplison = _remesh_params.eplison;
  auto lambda = _remesh_params.lambda;
  auto decay = _remesh_params.decay;

  auto fgm2 = [&](const Eigen::VectorXd &X, Eigen::VectorXd &g) {
    auto timer = Utils::ScopeTimer("迭代");
    eplison = eplison * decay;
    double lossCVT = 0, lossQE = 0, loss = 0;

    startRVD = clock();

    for (int i = 0; i < num; ++i) {
      Point_T query(X(i * 3 + 0), X(i * 3 + 1), X(i * 3 + 2)); // project to base surface
      Point_T closest = tree.closest_point(query);
      auto tri = tree.closest_point_and_primitive(query);

      Polyhedron::Face_handle f = tri.second;
      auto p1 = f->halfedge()->vertex()->point();
      auto p2 = f->halfedge()->next()->vertex()->point();
      auto p3 = f->halfedge()->next()->next()->vertex()->point();
      Eigen::Vector3d v1(p1.x(), p1.y(), p1.z());
      Eigen::Vector3d v2(p2.x(), p2.y(), p2.z());
      Eigen::Vector3d v3(p3.x(), p3.y(), p3.z());
      Eigen::Vector3d N = (v2 - v1).cross(v3 - v1);
      N.normalize();
      Nors[i] = N;
      BGAL::_Point3 p(closest.x(), closest.y(), closest.z());
      _sites[i] = p;
    }

    {
      // auto timer = Utils::ScopeTimer("_RVD.calculate_");
      _RVD.calculate_(_sites); // 速度限制步骤
    }
    Fnum++;
    if (_remesh_params.bOutputOnlyEnd == false)
      DoOutput(_sites, _RVD, num_sites, modelname, Fnum); // 1.3s

    endRVD = clock();
    RVDtime += (double)(endRVD - startRVD) / CLOCKS_PER_SEC;

    const auto &cells = _RVD.get_cells_();
    const auto &edges = _RVD.get_edges_();
    double energy = 0.0;
    g.setZero();
    vector<Eigen::Vector3d> gi;
    gi.resize(num);

    for (int i = 0; i < num; ++i) {
      gi[i] = Eigen::Vector3d(0, 0, 0);
    }
    int thread_num;
    thread_num = omp_get_max_threads();
    spdlog::debug("max threads: {}", thread_num);
    omp_set_num_threads(thread_num); // change to your CPU core numbers
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
      // auto timer = Utils::ScopeTimer("迭代 " + to_string(i));
      for (int j = 0; j < cells[i].size(); ++j) {
        auto &cell = cells[i][j];
        Eigen::VectorXd inte = BGAL::_Integral::integral_triangle3D(
            [&](BGAL::_Point3 p) {
              Eigen::VectorXd r(5);

              BGAL::_Point3 NorTriM = (_RVD.vertex_(cell[1]) - _RVD.vertex_(cell[0])).cross_(_RVD.vertex_(cell[2]) - _RVD.vertex_(cell[0]));
              NorTriM.normalized_();

              BGAL::_Point3 Nori(Nors[i].x(), Nors[i].y(), Nors[i].z());

              r[0] = (eplison * _rho(p) * ((_sites[i] - p).sqlength_())); // CVT

              r[1] = lambda * (NorTriM.dot_(p - _sites[i])) * (NorTriM.dot_(p - _sites[i])) + eplison * ((p - _sites[i]).sqlength_()); // qe+CVT

              r[2] = lambda * -2 * NorTriM.x() * (NorTriM.dot_(p - _sites[i])) + eplison * -2 * (p - _sites[i]).x(); // g
              r[3] = lambda * -2 * NorTriM.y() * (NorTriM.dot_(p - _sites[i])) + eplison * -2 * (p - _sites[i]).y(); // g
              r[4] = lambda * -2 * NorTriM.z() * (NorTriM.dot_(p - _sites[i])) + eplison * -2 * (p - _sites[i]).z(); // g

              return r;
            },
            _RVD.vertex_(cell[0]), _RVD.vertex_(cell[1]), _RVD.vertex_(cell[2]));
        // energy += alpha * inte(1);
        lossCVT += alpha * inte[0];
        loss += alpha * inte[1];
        gi[i].x() += alpha * inte[2];
        gi[i].y() += alpha * inte[3];
        gi[i].z() += alpha * inte[4];
      }
    }

    for (int i = 0; i < num; i++) {
      gi[i] = gi[i] - Nors[i] * (gi[i].dot(Nors[i]) / Nors[i].dot(Nors[i]));
      g(i * 3 + 0) += gi[i].x();
      g(i * 3 + 1) += gi[i].y();
      g(i * 3 + 2) += gi[i].z();
    }
    energy += loss;
    spdlog::info("迭代 energy: {:.7f} LossCVT: {:.7f} LossQE: {:.7f} Lambda_CVT: {:.7f}", energy, lossCVT / eplison, loss - lossCVT, eplison);
    return energy;
  };

  vector<Eigen::Vector3d> Pts2;

  Pts2 = Pts;
  num = Pts2.size();
  spdlog::info("顶点数: {}", num);
  _sites.resize(num);
  _para.max_linearsearch = 20;
  _para.max_iteration = 50;
  BGAL::_LBFGS lbfgs2(_para);
  Eigen::VectorXd iterX2(num * 3);

  for (int i = 0; i < num; ++i) {
    iterX2(i * 3 + 0) = Pts2[i].x();
    iterX2(i * 3 + 1) = Pts2[i].y();
    iterX2(i * 3 + 2) = Pts2[i].z();
    _sites[i] = BGAL::_Point3(Pts2[i].x(), Pts2[i].y(), Pts2[i].z());
  }

  _RVD.calculate_(_sites);
  start = clock();
  lbfgs2.minimize(fgm2, iterX2);
  end = clock();
  allTime += (double)(end - start) / CLOCKS_PER_SEC;

  spdlog::info("总耗时: {:.2f}s\tRVD耗时: {:.2f}s\tL-BFGS耗时: {:.2f}s", allTime, RVDtime, allTime - RVDtime);

  auto timer = Utils::ScopeTimer("写入模型");

  for (int i = 0; i < num; ++i) {
    // Point_T query(x0[i * 3], x0[i * 3+1], x0[i * 3+2]);
    Point_T query(iterX2(i * 3 + 0), iterX2(i * 3 + 1), iterX2(i * 3 + 2));
    Point_T closest = tree.closest_point(query);
    auto tri = tree.closest_point_and_primitive(query);

    Polyhedron::Face_handle f = tri.second;
    auto p1 = f->halfedge()->vertex()->point();
    auto p2 = f->halfedge()->next()->vertex()->point();
    auto p3 = f->halfedge()->next()->next()->vertex()->point();
    Eigen::Vector3d v1(p1.x(), p1.y(), p1.z());
    Eigen::Vector3d v2(p2.x(), p2.y(), p2.z());
    Eigen::Vector3d v3(p3.x(), p3.y(), p3.z());
    Eigen::Vector3d N = (v2 - v1).cross(v3 - v1);
    N.normalize();
    Nors[i] = N;

    _sites[i] = BGAL::_Point3(closest.x(), closest.y(), closest.z());
  }

  _RVD.calculate_(_sites);

  DoOutput(_sites, _RVD, num_sites, modelname, 2);
}
} // namespace BGAL
