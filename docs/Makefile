.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard mrkd/expl/*.ipynb))

.ipynb.md:
	jupyter nbconvert --to markdown $*.ipynb
	perl -pi -e "s/\\$$\\$$/\\$$/g" $*.md
#	[ ! -f $*.md ] || mv $*.md ./mrkd/expl/
#	[ ! -d $*_files ] || cp -rf $*_files ./mrkd/expl/ && rm -rf $*_files

all: $(MD) MomentTools.jl
	julia --color=yes make.jl

MomentTools.jl:
	git clone -b gh-pages https://github.com/AlgebraicGeometricModeling/MomentTools.jl.git

publish:
	cd MomentTools.jl; git add docs; git commit -am "update gh-pages"; git push

clean:
	rm ./mrkd/expl/*.md

