# Makefile
DEFENSE=WilliamDawn-QE2-presentation

.PHONY : all defense

all : defense

defense : $(DEFENSE).pdf

$(DEFENSE).pdf : $(DEFENSE).tex beamerthemeNCSU.sty\
	../WilliamDawn-QE2.bib
	pdflatex --draftmode --shell-escape $(DEFENSE)
	biber $(DEFENSE)
	makeglossaries $(DEFENSE)
	pdflatex --draftmode --shell-escape $(DEFENSE)
	pdflatex --shell-escape $(DEFENSE)
	grep -i "Warn" $(DEFENSE).log

forcedefense :
	rm $(DEFENSE).pdf
	$(MAKE) defense

clean : 
	rm $(DEFENSE_AUX) $(DEFENSE_INTERMEDIATES)
