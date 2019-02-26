DOKUMENT    = main
BUILDDIR    = build
PDFLATEX    = pdflatex
PDFLATEXOPT = --output-directory=$(BUILDDIR) -include-directory=$(BUILDDIR) -interaction=nonstopmode -extra-mem-top=20000000 --shell-escape
PDFLATEXCMD = $(PDFLATEX) $(PDFLATEXOPT)
BIBTEX      = biber
BIBTEXOPT   = --output-directory=build --input_directory=build
BIBTEXCMD   = $(BIBTEX) $(BIBTEXOPT)
MKINDEX     = makeindex
VIEWER      = evince

all: $(BUILDDIR)/$(DOKUMENT).pdf

$(BUILDDIR)/%.pdf: %.tex
		$(PDFLATEXCMD) $<
		$(BIBTEXCMD)  $(BUILDDIR)/$*
		$(PDFLATEXCMD) $<
		$(PDFLATEXCMD) $<

# Bibliography
bib:
		$(BIBTEXCMD)  $(BUILDDIR)/$*

# Create the basic PDF only. The bibtex and nomenclature procedure is also left.
# This mode is usefull, while this document is edited.
quick:
		$(PDFLATEXCMD) $(DOKUMENT)

# Display the PDF with the user chosen PDF-Viewer.
display: all
		$(VIEWER) $(DOKUMENT).pdf

# Clean the directory. All PDF procedure files and the PDF will be removed!
clean:
		$(RM) $(BUILDDIR)/*
		find . -name "*.aux" -exec $(RM) '{}' +
		find . -name "*.log" -exec $(RM) '{}' +
