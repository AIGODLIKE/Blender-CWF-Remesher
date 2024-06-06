#pragma once
#include <Eigen/Core>
#include <vector>
#include <string>

namespace UtilsSampling {
using Vec3 = std::vector<float>;

/// @param radius : 采样点最小半径. radius <= 0 时用 n_samples 参数来确定近似半径.
/// @param n_samples : 当 radius > 0 忽略, 否则尝试通过 合适的 radius 近似计算
/// @param verts : 顶点列表
/// @param nors : 法线列表(和顶点对应)
/// @param tris : 三角形索引列表
/// @endcode
/// @param [out] samples_pos : resulting samples positions
/// @param [out] samples_nors : resulting samples normals associated to samples_pos[]
/// @warning undefined behavior if (radius <= 0 && n_samples == 0) == true
void poisson_disk(float radius,
                  int n_samples,
                  const std::vector<Vec3> &verts,
                  const std::vector<Vec3> &nors,
                  const std::vector<int> &tris,
                  std::vector<Vec3> &samples_pos,
                  std::vector<Vec3> &samples_nor);

void poisson_disk(float radius,
                  int n_samples,
                  const Eigen::MatrixXf &V,
                  const Eigen::MatrixXf &N,
                  const Eigen::MatrixXi &F,
                  std::vector<Vec3> &samples_pos,
                  std::vector<Vec3> &samples_nor);

void poisson_disk_to_file(float radius, 
                          int &n_samples, 
                          const Eigen::MatrixXf &V, 
                          const Eigen::MatrixXf &N, 
                          const Eigen::MatrixXi &F, 
                          std::function<std::string(int)> get_name);
} // namespace UtilsSampling
