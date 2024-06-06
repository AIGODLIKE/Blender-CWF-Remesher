add_rules("mode.debug", "mode.release")

target("Algorithm")
    set_kind("static")
    add_files("./**.cpp")
    add_includedirs(".")