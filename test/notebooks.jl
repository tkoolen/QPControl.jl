let
    notebookdir = joinpath(@__DIR__, "..", "notebooks")
    excludedirs = [".ipynb_checkpoints"]
    excludefiles = String[]
    for (root, dir, files) in walkdir(notebookdir)
        basename(root) in excludedirs && continue
        for file in files
            file in excludefiles && continue
            name, ext = splitext(file)
            lowercase(ext) == ".ipynb" || continue
            path = joinpath(root, file)
            Pkg.activate(@__DIR__)
            @eval module $(gensym()) # Each notebook is run in its own module.
                using Test
                using NBInclude
                println("Testing $($(name))")
                # Note: use #NBSKIP in a cell to skip it during tests.
                @nbinclude($path; regex = r"^((?!\#NBSKIP).)*$"s)
            end # module
        end
    end
end
