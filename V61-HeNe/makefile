all:  build/main.pdf

build:
	mkdir -p build

build/main.pdf: build/stabilitaet.pdf build/polarisation.pdf build/T00.pdf build/T01.pdf build/stab-rund.pdf build/stab-flach.pdf| build

build/stabilitaet.pdf: python/stabilitaet.py | build
	python python/stabilitaet.py

build/polarisation.pdf: python/polarisation.py | build
	python python/polarisation.py

build/T00.pdf: python/tem00.py | build
	python python/tem00.py

build/T01.pdf: python/tem10.py | build
	python python/tem10.py

build/stab-rund.pdf: python/stab-rund.py | build
	python python/stab-rund.py

build/stab-flach.pdf: python/stab-rund.py | build
	python python/stab-rund.py

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
