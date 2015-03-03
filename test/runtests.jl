tests = ["test-pval-adjustment",
         "test-pi0"
         ]

for t in tests
    test_file = "$(t).jl"
    println(" * $(t) ...")
    include(test_file)
end
