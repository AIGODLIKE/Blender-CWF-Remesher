add_rules("mode.debug", "mode.release")

add_requires("libigl")

target("MAIN")
    set_kind("binary")
    add_files("main.cpp")
    add_includedirs("../include")
    add_packages("openmp", "libigl")
    add_deps("Algorithm", "BaseShape", "Draw", "Geodesic", "Integral", "Model", "Optimization", "PointCloudProcessing", "PQP", "Reconstruction", "Tessellation2D", "Tessellation3D", "CVTLike")