SOURCE = src/LucasStokey.ipynb
WEBDIR = docs/
SCRIPTDIR = bin/
PDFDIR = pdf/
OUTPUTS= docs/index.html bin/LucasStokey.jl pdf/LucasStokey.pdf
.PHONY: clean

all: $(OUTPUTS)

docs/index.html: $(SOURCE)
	jupyter nbconvert --to HTML $^
	mv src/LucasStokey.html docs/index.html 

bin/LucasStokey.jl: $(SOURCE)
	jupyter nbconvert --to script $^ --output-dir $(SCRIPTDIR)


pdf/LucasStokey.pdf: $(SOURCE)
	jupyter nbconvert --to Markdown $^ --template=src/hidecode.tpl --output-dir $(PDFDIR)
	$(MAKE) -C $(PDFDIR)


clean:
	rm -rf $(OUTPUTS)
	rm -rf $(PDFDIR)Lucas*
