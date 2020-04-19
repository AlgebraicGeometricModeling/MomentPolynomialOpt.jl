.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard ../../expl/*.ipynb))

.ipynb.md:
	jupiter nbconvert --to markdown $*.ipynb

all: $(MD)
