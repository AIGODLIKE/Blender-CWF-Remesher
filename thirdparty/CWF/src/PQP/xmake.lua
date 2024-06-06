add_rules("mode.debug", "mode.release")

target("PQP")
    set_kind("static")
    add_files("./**.cpp")
    add_includedirs(".")
