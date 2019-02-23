pdflatex -synctex=1 -file-line-error-style -interaction=nonstopmode --output-directory=build -include-directory=build -extra-mem-top=20000000 -file-line-error main.tex
biber --output-directory=build --input_directory=build main
pdflatex -synctex=1 -file-line-error-style -interaction=nonstopmode --output-directory=build -include-directory=build -extra-mem-top=20000000 -file-line-error main.tex
pdflatex -synctex=1 -file-line-error-style -interaction=nonstopmode --output-directory=build -include-directory=build -extra-mem-top=20000000 -file-line-error main.tex
