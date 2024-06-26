add_rules("mode.debug", "mode.release")
add_requires("openmp", "spdlog")

target("CVTLike")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("boost", "openmp", "spdlog")
    add_deps("Algorithm", "BaseShape", "Model", "Tessellation2D", "Tessellation3D", "Optimization")
    add_defines("NOMINMAX")
    add_includedirs(".")
