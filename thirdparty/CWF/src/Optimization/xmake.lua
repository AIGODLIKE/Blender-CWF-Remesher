add_rules("mode.debug", "mode.release")

-- add_requires("spdlog")

target("Optimization")
    set_kind("static")
    add_files("./**.cpp")
    -- add_packages("spdlog")
    add_includedirs(".")
