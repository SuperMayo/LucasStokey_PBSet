SOURCE = src/lucasStokey.ipynb
WEBDIR = docs/
SCRIPTDIR = bin/
PDFDIR = pdf/
OUTPUTS= docs/index.html bin/LucasStokey.jl pdf/LucasStokey.pdf
.PHONY: clean

all: $(OUTPUTS)

docs/index.html: $(SOURCE)
	jupyter nbconvert --to HTML $^ --output-dir $(WEBDIR)

bin/LucasStokey.jl: $(SOURCE)
	jupyter nbconvert --to script $^ --output-dir $(SCRIPTDIR)

pdf/LucasStokey.pdf: $(SOURCE)
	jupyter nbconvert --to pdf --template template $^ --output-dir $(PDFDIR)

clean:
	rm -rf $(OUTPUTS)
