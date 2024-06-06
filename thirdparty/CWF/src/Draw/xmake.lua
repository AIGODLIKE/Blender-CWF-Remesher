add_rules("mode.debug", "mode.release")

target("Draw")
    set_kind("static")
    add_files("./**.cpp")
    add_packages("boost")
    add_deps("Algorithm", "BaseShape")
    add_includedirs(".")
