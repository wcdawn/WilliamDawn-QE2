# Makefile
NAME=WilliamDawn-QE2

.PHONY : all text

all : text

text : $(NAME).pdf

$(NAME).pdf : $(NAME).tex $(NAME).bib
	pdflatex --draftmode $(NAME)
	biber $(NAME)
	makeglossaries $(NAME)
	pdflatex --draftmode $(NAME)
	pdflatex $(NAME)
	grep -i "Warn" $(NAME).log

forcetext :
	rm $(NAME).pdf
	$(MAKE) text

clean :
	rm $(AUX) $(INTERMEDIATES)
