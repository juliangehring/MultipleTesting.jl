module MultipleTestingDoc

using MultipleTesting
using Lexicon

cd( joinpath(Pkg.dir("MultipleTesting"), "doc") )

config = Config(md_permalink = false,
                md_subheader = :category,
                include_internal = false,
                metadata_order = Symbol[],
                mathjax = true)
index = save("api.md", MultipleTesting, config);
save("index.md", Index([index]));

## notebooks
nb_files = ["pi0-estimation"]
cd( joinpath(Pkg.dir("MultipleTesting"), "doc", "notebooks") )

for nb in nb_files
    run(`jupyter nbconvert --to markdown "$nb.ipynb"`)
end

end
