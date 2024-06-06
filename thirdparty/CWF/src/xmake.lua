add_rules("mode.debug", "mode.release")

add_packages("eigen", "spdlog", "openmp", "cgal", "libigl", { public = true })

includes("Algorithm")
includes("Optimization")
includes("PQP")

includes("BaseShape")

includes("Draw")
includes("Integral")
includes("PointCloudProcessing")
includes("Tessellation2D")

includes("Model")

includes("Geodesic")
includes("Reconstruction")
includes("Tessellation3D")
includes("CVTLike")