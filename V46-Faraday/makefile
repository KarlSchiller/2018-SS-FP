all:  build/main.pdf

build:
	mkdir -p build

build/main.pdf: build/differenzen.pdf build/feldmessung2.pdf build/winkel.pdf | build

build/differenzen.pdf build/feldmessung2.pdf build/winkel.pdf: auswertung.py \
rohdaten/feldmessung2.txt rohdaten/probe1.txt rohdaten/probe2.txt \
rohdaten/probe3.txt | build
	python auswertung.py

light:
	lualatex -output-directory=build -interaction=nonstopmode main.tex

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
