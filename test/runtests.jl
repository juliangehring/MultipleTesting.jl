tests = ["test-pval-adjustment",
         "test-pval-pi0-adjustment",
         "test-pi0",
         "test-grenander",
         "test-utils"
         ]

for t in tests
    test_file = "$(t).jl"
    println(" * $(t) ...")
    include(test_file)
end
