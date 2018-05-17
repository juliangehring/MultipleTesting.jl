.PHONY: test docs

test:
	julia --color=yes test/runtests.jl

docs:
	julia --color=yes docs/make.jl

