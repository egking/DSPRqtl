#!/bin/bash

RMDFILE=DSPRqtl-intro
Rscript -e "require(knitr); require(markdown); knit('DSPRqtl-intro.Rmd', 'DSPRqtl-intro.md')"
pandoc -N --include-in-header=header.tex -s --from=markdown DSPRqtl-intro.md -o DSPRqtl-intro.pdf
rm DSPRqtl-intro.md
