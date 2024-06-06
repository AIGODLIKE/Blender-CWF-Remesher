add_rules("mode.debug", "mode.release")

add_rules("plugin.compile_commands.autoupdate", {outputdir = ".vscode"})

set_project("BGAL")

set_version("1.0")

set_languages("c++17")
set_encodings("source:utf-8")

if is_plat("windows") and is_mode("debug") then
    set_runtimes("MT")
    -- set_runtimes("MD")
end
-- add_defines("PYTHON_EXECUTABLE=" .. "\"/Users/karrycharon/.pyenv/shims/python\"", { public = true })

add_includedirs("include", { public = true })

-- set_targetdir("build")
-- set_objectdir("build/.objs")
-- set_dependir("build/.deps")

includes("src")

includes("remesh")

includes("poisson_disk")

includes("MAIN")
