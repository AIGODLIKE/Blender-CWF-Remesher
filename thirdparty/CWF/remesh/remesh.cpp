#include <filesystem>
#include <functional>
#include <string>
#include <thread>

#include <omp.h>

#include <BGAL/CVTLike/CPD.h>
#include <BGAL/CVTLike/CVT.h>
#include <BGAL/Integral/Integral.h>
#include <BGAL/Model/ManifoldModel.h>
#include <BGAL/Model/Model_Iterator.h>
#include <BGAL/Optimization/ALGLIB/optimization.h>
#include <BGAL/Optimization/LBFGS/LBFGS.h>
#include <BGAL/Optimization/LinearSystem/LinearSystem.h>
#include <BGAL/PointCloudProcessing/PoissonDiskSampler.h>
#include <BGAL/Tessellation3D/Tessellation3D.h>

#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>

#include <Utils/Common.hpp>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

#include "utils_sampling.hpp"
#include <pybind11/eval.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using Path = std::filesystem::path;

void Prepare(BGAL::RemeshParams &params) {
  Eigen::MatrixXf V; // 顶点列表
  Eigen::MatrixXf N; // 顶点法线
  Eigen::MatrixXi F; // 面列表
  std::string file = params.meshName;
  auto outputDir = params.OutputDirPath();
  {
    auto timer = Utils::ScopeTimer("读取Obj");
    igl::readOBJ(params.ObjPath().string(), V, F);
  }

  {
    auto timer = Utils::ScopeTimer("计算法线");
    igl::per_vertex_normals(V, F, N);
  }

  {
    auto timer = Utils::ScopeTimer("泊松采样");
    auto get_name = [&](int i) { return (outputDir / fmt::format("n{}_{}_inputPoints.xyz", i, file)).string(); };
    UtilsSampling::poisson_disk_to_file(0, params.samples, V, N, F, get_name);
  }
  auto tempObj = outputDir / "Temp.obj";
  auto tempOff = outputDir / "Temp.off";
  igl::writeOFF(tempOff.string(), V, F);
  igl::writeOBJ(tempObj.string(), V, F);
}

void DoRemesh(BGAL::RemeshParams params, bool prepare = true) {
  params.Prepare();

  auto objPath = params.ObjPath();
  auto outputDir = params.OutputDirPath();
  // .obj to .off
  if (prepare)
    Prepare(params);
  int &nSamples = params.samples;
  spdlog::info("当前文件 = {}, 采样 = {}", params.meshName, nSamples);
  auto tempObj = outputDir / ("Temp.obj");
  BGAL::_ManifoldModel model(tempObj.string()); // 0.31s

  std::function<double(BGAL::_Point3 & p)> rho = [](BGAL::_Point3 &p) { return 1; };
  BGAL::_LBFGS::_Parameter para;
  para.is_show = true;
  para.epsilon = 0.0001;

  BGAL::_CVT3D cvt(model, rho, para); // 0.36s
  cvt._remesh_params = params;
  auto xyzname = fmt::format("n{}_{}_inputPoints.xyz", nSamples, params.meshName); //
  cvt.loadpoints(xyzname);                                                         // 0.004s
  cvt.calculatet_(nSamples, params.meshName);
}

void RemeshTest() {
  // your can try: Nums = 600      1000   1000          2000
  //               File = mobius1  block  block_smooth  bunny
  BGAL::RemeshParams params;
  params.samples = 5100;
  params.meshName = "bunny";
  params.workDir = "/Volumes/Code/CPP/CWF/data/";
  DoRemesh(params);
}

void DoRemeshEx(BGAL::RemeshParams params, bool prepare) {
  try {
    DoRemesh(params, prepare);
  } catch (const std::exception &e) {
    spdlog::error("{}", e.what());
    params.finishedCallback(params.samples, params.meshName);
  }
}

void Remesh(BGAL::RemeshParams &params) {
  spdlog::set_pattern("%^[%l]%$ %v");
  // DoRemesh(params);
  // 在子线程执行 DoRemesh
  std::thread t(DoRemeshEx, params, true);
  t.detach();
}

// BGAL::RemeshParams 封装成py对应的class
PYBIND11_MODULE(remesh, m) {
  m.def("remesh", &Remesh, "Remesh");

  py::class_<BGAL::RemeshParams>(m, "RemeshParams")
      .def(py::init<>())
      .def_readwrite("batchRemeshCallback", &BGAL::RemeshParams::batchRemeshCallback)
      .def_readwrite("finishedCallback", &BGAL::RemeshParams::finishedCallback)
      .def_readwrite("outputDir", &BGAL::RemeshParams::outputDir)
      .def_readwrite("meshName", &BGAL::RemeshParams::meshName)
      .def_readwrite("workDir", &BGAL::RemeshParams::workDir)
      .def_readwrite("_samples", &BGAL::RemeshParams::samples)
      .def_readwrite("_fnum", &BGAL::RemeshParams::fnum)
      .def_readwrite("_alpha", &BGAL::RemeshParams::alpha)
      .def_readwrite("_eplison", &BGAL::RemeshParams::eplison)
      .def_readwrite("_lambda", &BGAL::RemeshParams::lambda)
      .def_readwrite("_decay", &BGAL::RemeshParams::decay)
      .def_readwrite("bOutputOnlyEnd", &BGAL::RemeshParams::bOutputOnlyEnd)
      .def_readwrite("bOutputXyz", &BGAL::RemeshParams::bOutputXyz)
      .def_readwrite("bOutputRemesh", &BGAL::RemeshParams::bOutputRemesh)
      .def_readwrite("bOutputRVD", &BGAL::RemeshParams::bOutputRVD);
}

int main(int argc, char *argv[]) {
  spdlog::set_pattern("%^[%l]%$ %v");
  auto timer = Utils::ScopeTimer("总计");
  string filepath = "/Users/karrycharon/Desktop/xyz/test.obj";
  // filepath = argv[0];
  spdlog::info("测试函数 RemeshTest");
  // Remesh(filepath, 2356);
  RemeshTest();
  // poisson_test();
  spdlog::info("成功!");
  return 0;
}