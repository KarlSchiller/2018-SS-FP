all:  build/main.pdf

build:
	mkdir -p build

build/main.pdf: build/Nullmessung1.pdf build/Nullmessung2.pdf build/W1_pos1.pdf \
rohdaten/kein_wuerfel.txt rohdaten/W4_I0.txt rohdaten/W1_pos1.txt \
rohdaten/Würfel_1.txt rohdaten/Würfel_2.txt rohdaten/Würfel_3.txt \
rohdaten/Würfel_4.txt | build

build/Nullmessung1.pdf build/Nullmessung2.pdf build/W1_pos1.pdf: python/berechnungen.py | build
	python python/berechnungen.py

build/main.pdf: FORCE | build
	  TEXINPUTS="$(call translate,build:)" \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
		--interaction=nonstopmode
		--halt-on-error
	main.tex

clean:
	rm -rf build

FORCE:

.PHONY: all clean
