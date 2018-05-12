all:  build/main.pdf

build:
	mkdir -p build

build/main.pdf: build/alpha.pdf build/cv.pdf | build

build/alpha.pdf build/cv.pdf: auswertung.py \
rohdaten/alpha.txt rohdaten/messung.txt | build
	python auswertung.py

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
