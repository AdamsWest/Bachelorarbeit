pdflatex -synctex=1 -file-line-error-style -interaction=nonstopmode --output-directory=build -include-directory=build -extra-mem-top=20000000 -file-line-error main.tex
biber --output-directory=build --input_directory=build main
makeindex -s build/main.ist -t build/main.alg -o build/main.acr build/main.acn
makeindex -s build/main.ist -t build/main.glg -o build/main.gls build/main.glo
makeindex -s build/main.ist -t build/main.llg -o build/main.lyi build/main.lyg
makeindex -s build/main.ist -t build/main.glg -o build/main.gyi build/main.gyg
makeindex -s build/main.ist -t build/main.ilg -o build/main.iyi build/main.iyg
pdflatex -synctex=1 -file-line-error-style -interaction=nonstopmode --output-directory=build -include-directory=build -extra-mem-top=20000000 -file-line-error main.tex
pdflatex -synctex=1 -file-line-error-style -interaction=nonstopmode --output-directory=build -include-directory=build -extra-mem-top=20000000 -file-line-error main.tex
