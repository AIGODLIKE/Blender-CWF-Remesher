add_rules("mode.debug", "mode.release")

target("Reconstruction")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("boost")
    add_deps("Algorithm", "BaseShape", "Model", "PQP")
    add_includedirs(".")
