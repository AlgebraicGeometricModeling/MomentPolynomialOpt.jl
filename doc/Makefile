.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard ../expl/*.ipynb))

.ipynb.md:
	jupyter nbconvert --to markdown $*.ipynb
	perl -pi -e "s/\\$$\\$$/\\$$/g" $*.md
	[ ! -d ./mrkd/expl ] || mkdir ./mrkd/expl
	[ ! -f $*.md ] || mv $*.md ./mrkd/expl/
	[ ! -d $*_files ] || cp -rf $*_files ./mrkd/expl/ && rm -rf $*_files

all: $(MD)
	julia --color=yes make.jl

clean:
	rm ./mrkd/expl/*.md
deploy:
	cp -r html/* /Volumes/Bernard.mourrain/software/MomentTools
