all:  build/main.pdf

build:
	mkdir -p build
	mkdir -p build/leck build/tab

build/main.pdf: build/volumina.txt build/evak_D.pdf build/evak_T.pdf \
	build/leck/D_0,1.pdf build/leck/T_5e-5.pdf | build

build/volumina.txt: python/volumen.py | build
	python python/volumen.py

build/evak_D.pdf: python/evak.py build/volumina.txt rohdaten/druck_dreh.txt | build
	python python/evak.py

build/evak_T.pdf: python/evak.py build/volumina.txt rohdaten/druck_turbo.txt | build
	python python/evak.py

build/leck/D_0,1.pdf build/leck/T_5e-5.pdf: python/leck.py build/volumina.txt \
rohdaten/leck_dreh_0,1mbar.txt rohdaten/leck_dreh_0,4mbar.txt \
rohdaten/leck_dreh_0,8mbar.txt rohdaten/leck_dreh_1,0mbar.txt \
rohdaten/leck_turbo_5e-5mbar.txt rohdaten/leck_turbo_1e-4mbar.txt \
rohdaten/leck_turbo_1,5e-4mbar.txt rohdaten/leck_turbo_2e-4mbar.txt | build
	python python/leck.py

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
