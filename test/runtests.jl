tests = ["test-pval-adjustment"
         #"test-qvalue",
         #"test-pi0"
         #"test-utils"
         ]

for t in tests
    test_file = "$(t).jl"
    println(" * $(t) ...")
    include(test_file)
end
