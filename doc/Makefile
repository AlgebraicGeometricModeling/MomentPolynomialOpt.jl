.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard ../expl/*.ipynb))

.ipynb.md:
	jupiter nbconvert --to markdown $*.ipynb --stdout > $*.md && cp $*.md ./mrkd/expl && rm $*.md

all: $(MD)
	julia --color=yes make.jl

clean:
	rm ./mrkd/expl/*.md
deploy:
	cp -r html/* /Volumes/Bernard.mourrain/software/MomentTools
