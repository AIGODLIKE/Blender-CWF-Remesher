add_rules("mode.debug", "mode.release")

add_requires("cgal")

target("Tessellation2D")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("boost", "cgal")
    add_deps("Algorithm", "BaseShape")
    add_includedirs(".")
