COLOR := yes

.PHONY: test docs

all: test docs

test:
	julia --color=$(COLOR) test/runtests.jl

docs:
	julia --color=$(COLOR) docs/make.jl

