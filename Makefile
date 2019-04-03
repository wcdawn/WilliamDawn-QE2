# Makefile
NAME=WilliamDawn-QE2
INTERMEDIATES=*.aux *.log *.out *.toc *.gls *.glg *.glo *.bbl *.bcf *.blg \
			  *.run.xml *.acn *.acr *.alg *.ist *.nav *.snm

.PHONY : all text forcetext presentation forcepresentation

all : text presentation

text : $(NAME).pdf

$(NAME).pdf : $(NAME).tex $(NAME).bib
	pdflatex --draftmode $(NAME)
	biber $(NAME)
	makeglossaries $(NAME)
	pdflatex --draftmode $(NAME)
	pdflatex $(NAME)
	grep -i "Warn" $(NAME).log

presentation:
	$(MAKE) -C ./presentation/ presentation

forcetext :
	rm $(NAME).pdf
	$(MAKE) text

forcepresentation :
	$(MAKE) -C ./presentation/ forcepresentation


clean :
	rm $(AUX) $(INTERMEDIATES)
	$(MAKE) -C ./presentation/ clean
