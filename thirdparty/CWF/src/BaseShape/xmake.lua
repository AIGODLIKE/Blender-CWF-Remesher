add_rules("mode.debug", "mode.release")

target("BaseShape")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("boost")
    add_deps("Algorithm")
    add_includedirs(".")