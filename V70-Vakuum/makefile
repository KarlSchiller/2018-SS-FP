all:  build/main.pdf

build:
	mkdir -p build

build/main.pdf: build/volumina.txt build/evak_D.pdf | build

build/volumina.txt: python/volumen.py | build
	python python/volumen.py

build/evak_D.pdf: python/evak.py build/volumina.txt rohdaten/druck_dreh.txt | build
	python python/evak.py

build/evak_T.pdf: python/evak.py build/volumina.txt rohdaten/druck_turbo.txt | build
	python python/evak.py

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
