tests = ["test-pval-adjustment",
         "test-pval-pi0-adjustment",
         "test-pi0-estimators",
         "test-grenander",
         "test-utils",
         "test-model",
         "test-types"
         ]

for t in tests
    test_file = "$(t).jl"
    include(test_file)
end
