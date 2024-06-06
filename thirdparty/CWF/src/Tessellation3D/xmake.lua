add_rules("mode.debug", "mode.release")

add_requires("boost", "openmp")

target("Tessellation3D")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("boost", "openmp")
    add_deps("Algorithm", "BaseShape", "Model")
    add_includedirs(".")
