add_rules("mode.debug", "mode.release")

target("Model")
    set_kind("static")
    add_files("*.cpp")
    add_packages("boost")
    add_deps("Algorithm", "BaseShape", "PQP")
    add_includedirs(".")
