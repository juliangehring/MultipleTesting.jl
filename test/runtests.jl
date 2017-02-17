tests = [
    "test-types",
    "test-pval-adjustment",
    "test-pval-pi0-adjustment",
    "test-grenander",
    "test-utils",
    "test-model",
    "test-combinations",
    "test-pi0-estimators",
    "test-higher-criticism"
]

for t in tests
    test_file = "$(t).jl"
    include(test_file)
end
