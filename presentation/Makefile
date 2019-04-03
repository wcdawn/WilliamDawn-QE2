# Makefile
PRESENTATION=WilliamDawn-QE2-presentation
INTERMEDIATES=*.aux *.log *.out *.toc *.gls *.glg *.glo *.bbl *.bcf *.blg \
			  *.run.xml *.acn *.acr *.alg *.ist *.nav *.snm

.PHONY : all presentation

all : presentation

presentation : $(PRESENTATION).pdf

$(PRESENTATION).pdf : $(PRESENTATION).tex beamerthemeNCSU.sty\
	../WilliamDawn-QE2.bib
	pdflatex --draftmode --shell-escape $(PRESENTATION)
	biber $(PRESENTATION)
	makeglossaries $(PRESENTATION)
	pdflatex --draftmode --shell-escape $(PRESENTATION)
	pdflatex --shell-escape $(PRESENTATION)
	grep -i "Warn" $(PRESENTATION).log

forcepresentation :
	rm $(PRESENTATION).pdf
	$(MAKE) presentation

clean : 
	rm $(PRESENTATION_INTERMEDIATES)
