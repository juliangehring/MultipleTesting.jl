#!/usr/bin/env sh

# rebuild all notebooks

# run two times to avoid precompilation in the final output
for i in 1 2
do
    echo "Pass ${i}..."
    # timout: allow long running cells
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=180 *.ipynb
done
