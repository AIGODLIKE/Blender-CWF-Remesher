add_rules("mode.debug", "mode.release")

set_languages("c++14")

add_requires("eigen", "boost", "cgal")

target("test")
    set_kind("binary")
    add_files("*.cpp")
    add_includedirs("../include/BGAL")
    add_packages("boost", "eigen", "cgal")
    add_deps("Algorithm", "BaseShape", "Draw", "Geodesic", "Integral", "Model", "Optimization", "PointCloudProcessing", "PQP", "Reconstruction", "Tessellation2D", "Tessellation3D", "CVTLike")