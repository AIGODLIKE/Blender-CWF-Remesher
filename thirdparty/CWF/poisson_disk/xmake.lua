add_rules("mode.debug", "mode.release")

add_requires("vcglib", "eigen", "spdlog", "cgal")

target("PoissonDisk")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("vcglib", "eigen", "spdlog", "cgal", { public = true })
    add_includedirs(".", { public = true })