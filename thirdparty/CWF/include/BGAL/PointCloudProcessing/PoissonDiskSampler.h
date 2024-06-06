#include <Eigen/Dense>
#include <vector>

// 三维点
using Point3D = Eigen::Vector3d;

namespace BGAL {
// 生成一个随机点在[min, max]之间
Point3D randomPoint3D(double min, double max);

// 检查点是否在网格内
bool isInsideGrid(const Point3D &p, double min, double max);

// 检查点是否在网格内
bool isInsideGrid(const Point3D &p, const Point3D &min, const Point3D &max);

// 检查点是否与已有点集中的点的距离大于给定最小距离
bool isFarEnough(const Point3D &p, const std::vector<Point3D> &points, double minDist);

// 计算点云的最小距离
double calculateMinDist(const Eigen::MatrixXd& vertices, int expectedPointCount);

// Poisson-disk采样算法
std::vector<Point3D> poissonDiskSampling(const std::vector<Point3D> &vertices, double minDist, int newPointsCount = 30);

// Poisson-disk采样算法，使用 Eigen::MatrixXd 作为输入
std::vector<Point3D> poissonDiskSampling(const Eigen::MatrixXd &vertices, int newPointsCount = 30);
}