#include <fstream>
#include <functional>
#include <iostream>

#ifdef _WIN32
#include <io.h>
#endif
#include <omp.h>
#include <random>

#include <BGAL/Integral/Integral.h>
#include <BGAL/Model/ManifoldModel.h>
#include <BGAL/Model/Model_Iterator.h>
#include <BGAL/PointCloudProcessing/PoissonDiskSampler.h>
#include <BGAL/Optimization/ALGLIB/optimization.h>
#include <BGAL/Optimization/LBFGS/LBFGS.h>
#include <BGAL/Optimization/LinearSystem/LinearSystem.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>

#include <BGAL/Tessellation3D/Tessellation3D.h>

#include <BGAL/CVTLike/CPD.h>
#include <BGAL/CVTLike/CVT.h>

void RemeshTest() {
  // your can try: Nums = 600      1000   1000          2000
  //               File = mobius1  block  block_smooth  bunny

  int Nums = 2000;
  string file = "bunny";
  cout << "Now file: " << file << endl;

  cout << Nums << "   " << file << "   \n";
  string filepath = "../../../../data/";
  string modelname = file;

  // .obj to .off
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ(filepath + modelname + ".obj", V, F);

  igl::writeOFF("Temp.off", V, F);
  igl::writeOBJ("Temp.obj", V, F);

  string xyzpath = "../../../../data/n" + to_string(Nums) + "_" + modelname + "_inputPoints1.xyz";
  auto points = BGAL::poissonDiskSampling(V, Nums);
  // poission disk 采样 Nums 个点 写入 xyz中

  BGAL::_ManifoldModel model("Temp.obj");

  std::function<double(BGAL::_Point3 & p)> rho = [](BGAL::_Point3 &p) { return 1; };

  BGAL::_LBFGS::_Parameter para;
  para.is_show = true;
  para.epsilon = 0.0001;
  BGAL::_CVT3D cvt(model, rho, para);
  int num = Nums;
  cvt.calculate_(num, (char *)modelname.c_str());
}

int alltest() {
  std::cout << "====================CWF3DTest" << std::endl;
  RemeshTest();
  std::cout << "successful!" << std::endl;
  return 0;
}

int main() {
  alltest();
  return 0;
}