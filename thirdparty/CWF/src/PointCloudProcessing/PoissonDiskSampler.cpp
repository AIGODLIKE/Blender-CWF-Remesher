#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <BGAL/PointCloudProcessing/PoissonDiskSampler.h>
#include <vector>

// 三维点
using Point3D = Eigen::Vector3d;

// 生成随机数的引擎
std::default_random_engine rng(std::random_device{}());

// 生成一个随机点在[min, max]之间
Point3D BGAL::randomPoint3D(double min, double max) {
  std::uniform_real_distribution<double> dist(min, max);
  return Point3D(dist(rng), dist(rng), dist(rng));
}

// 检查点是否在网格内
bool BGAL::isInsideGrid(const Point3D &p, double min, double max) { return (p.x() >= min && p.x() <= max && p.y() >= min && p.y() <= max && p.z() >= min && p.z() <= max); }

// 检查点是否在网格内
bool BGAL::isInsideGrid(const Point3D &p, const Point3D &min, const Point3D &max) {
  return (p.x() >= min.x() && p.x() <= max.x() && p.y() >= min.y() && p.y() <= max.y() && p.z() >= min.z() && p.z() <= max.z());
}

// 检查点是否与已有点集中的点的距离大于给定最小距离
bool BGAL::isFarEnough(const Point3D &p, const std::vector<Point3D> &points, double minDist) {
  for (const auto &point : points) {
    if ((p - point).norm() < minDist) {
      return false;
    }
  }
  return true;
}

// 计算点云的最小距离
double BGAL::calculateMinDist(const Eigen::MatrixXd& vertices, int expectedPointCount) {
    Point3D min(vertices.col(0).minCoeff(), vertices.col(1).minCoeff(), vertices.col(2).minCoeff());
    Point3D max(vertices.col(0).maxCoeff(), vertices.col(1).maxCoeff(), vertices.col(2).maxCoeff());
    double diagonalLength = (max - min).norm();
    return diagonalLength / std::sqrt(expectedPointCount);
}

// Poisson-disk采样算法
std::vector<Point3D> BGAL::poissonDiskSampling(const std::vector<Point3D> &vertices, double minDist, int newPointsCount) {
  std::vector<Point3D> result;
  if (vertices.empty())
    return result;

  // 定义3D网格的范围
  double minX = vertices[0].x(), maxX = vertices[0].x();
  double minY = vertices[0].y(), maxY = vertices[0].y();
  double minZ = vertices[0].z(), maxZ = vertices[0].z();

  for (const auto &vertex : vertices) {
    if (vertex.x() < minX)
      minX = vertex.x();
    if (vertex.x() > maxX)
      maxX = vertex.x();
    if (vertex.y() < minY)
      minY = vertex.y();
    if (vertex.y() > maxY)
      maxY = vertex.y();
    if (vertex.z() < minZ)
      minZ = vertex.z();
    if (vertex.z() > maxZ)
      maxZ = vertex.z();
  }

  // 随机选择一个初始点
  std::uniform_int_distribution<> dist(0, vertices.size() - 1);
  Point3D initialPoint = vertices[dist(rng)];
  result.push_back(initialPoint);

  // 活跃点列表
  std::vector<Point3D> activeList = {initialPoint};

  while (!activeList.empty()) {
    std::uniform_int_distribution<> activeDist(0, activeList.size() - 1);
    Point3D point = activeList[activeDist(rng)];
    bool found = false;

    for (int i = 0; i < newPointsCount; ++i) {
      double r = minDist * (1 + rng() / (double)rng.max()); // [minDist, 2 * minDist)
      Point3D newPoint = point + randomPoint3D(-r, r);

      if (isInsideGrid(newPoint, minX, maxX) && isFarEnough(newPoint, result, minDist)) {
        result.push_back(newPoint);
        activeList.push_back(newPoint);
        found = true;
        break;
      }
    }

    if (!found) {
      activeList.erase(activeList.begin() + activeDist(rng));
    }
  }

  return result;
}

// Poisson-disk采样算法，使用 Eigen::MatrixXd 作为输入
std::vector<Point3D> BGAL::poissonDiskSampling(const Eigen::MatrixXd &vertices, int newPointsCount) {
  double minDist = 0;
  std::vector<Point3D> result;
  if (vertices.rows() == 0)
    return result;

  // 定义3D网格的范围
  Point3D min(vertices.col(0).minCoeff(), vertices.col(1).minCoeff(), vertices.col(2).minCoeff());
  Point3D max(vertices.col(0).maxCoeff(), vertices.col(1).maxCoeff(), vertices.col(2).maxCoeff());

  // 随机选择一个初始点
  std::uniform_int_distribution<> dist(0, vertices.rows() - 1);
  Point3D initialPoint = vertices.row(dist(rng));
  result.push_back(initialPoint);

  // 活跃点列表
  std::vector<Point3D> activeList = {initialPoint};
  int findCount = 0;
  while (!activeList.empty()) {
    std::uniform_int_distribution<> activeDist(0, activeList.size() - 1);
    Point3D point = activeList[activeDist(rng)];
    bool found = false;

    for (int i = 0; i < newPointsCount; ++i) {
      double r = minDist * (1 + rng() / (double)rng.max()); // [minDist, 2 * minDist)
      Point3D newPoint = point + randomPoint3D(-r, r);

      if (isInsideGrid(newPoint, min, max) && isFarEnough(newPoint, result, minDist)) {
        result.push_back(newPoint);
        activeList.push_back(newPoint);
        found = true;
        std::cout << "FIND:" << findCount++ << std::endl;
        break;
      }
    }

    if (!found) {
      activeList.erase(activeList.begin() + activeDist(rng));
    }
  }

  return result;
}