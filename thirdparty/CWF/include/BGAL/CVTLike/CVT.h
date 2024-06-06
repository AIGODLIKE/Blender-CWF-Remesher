#pragma once
#include <functional>
#ifdef _HAS_STD_BYTE
#undef _HAS_STD_BYTE
#endif
#define _HAS_STD_BYTE 0

#include "BGAL/BaseShape/Line.h"
#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
#include "BGAL/BaseShape/Triangle.h"
#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"
#include "BGAL/Optimization/LBFGS/LBFGS.h"
#include "BGAL/Tessellation3D/Tessellation3D.h"

using namespace std;

namespace BGAL {
struct RemeshParams {
  string workDir = "./";
  string meshName = "";
  string outputDir = "";
  int samples = 1000;
  int fnum = 4;
  double alpha = 1.0;
  double eplison = 1.0;
  double lambda = 1.0; // eplison is CVT weight, lambda is qe weight.
  double decay = 0.95;
  bool bOutputOnlyEnd = false;
  bool bOutputXyz = true;
  bool bOutputRemesh = true;
  bool bOutputRVD = true;

  std::function<void(int, string)> batchRemeshCallback = [](int step, string out) { spdlog::info("Remeshing {} => {}", step, out); };

  std::function<void(int, string)> finishedCallback = [](int step, string out) { spdlog::info("Remeshing {} => {}", step, out); };

  void Prepare() { outputDir = outputDir == "" ? workDir : outputDir; }

  std::filesystem::path WorkDirPath() { return workDir; }
  std::filesystem::path ObjPath() { return WorkDirPath().append(meshName).replace_extension(".obj"); }
  std::filesystem::path OutputDirPath() { return outputDir; }
};

class _CVT3D {
public:
  _CVT3D(const _ManifoldModel &model);
  _CVT3D(const _ManifoldModel &model, std::function<double(_Point3 &p)> &rho, _LBFGS::_Parameter para);
  void calculate_(int site_num, char *modelName);
  void calculatet_(int site_num, string modelname);
  void DoOutput(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step);
  void OutputXYZ(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step);
  void OutputRVD(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step);
  void OutputRemesh(std::vector<_Point3> &sites, _Restricted_Tessellation3D RVD, int num, string modelname, int step);
  void loadpoints(string xyzpath = "");

  const std::vector<_Point3> &get_sites() const { return _sites; }
  const _Restricted_Tessellation3D &get_RVD() const { return _RVD; }
  const _Restricted_Tessellation3D &get_RVD2() const { return _RVD2; }

public:
  const _ManifoldModel &_model;
  vector<Eigen::Vector3d> Pts;
  vector<Eigen::Vector3d> Nors;
  _Restricted_Tessellation3D _RVD;
  _Restricted_Tessellation3D _RVD2;
  std::vector<_Point3> _sites{};
  std::function<double(_Point3 &p)> _rho;
  _LBFGS::_Parameter _para;
  RemeshParams _remesh_params;

private:
  bool loaded = false;
};
} // namespace BGAL
