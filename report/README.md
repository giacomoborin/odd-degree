#### Compilation

`latexmk`  will call `pdflatex`/`bibtex` as often as needed automatically. It is
the default compilation tool used by Overleaf (but is much older than that,
1998) written in an unholy 12k lines of Perl.

Together with the configuration file `latexmkrc` supplied in this repository, it
will not litter auxiliary files everywhere.

Just install `latexmk` (`texlive-binextra` on archlinux and `latexmk` on ubuntu)
and call

    latexmk main.tex

It will do its thing and create `main.pdf` and a directory `.latexmk.aux`
(inside the directory where the tex source file lives, not your current working
directory) with all the auxiliary files.

To `make clean && make` you simply add `-gg` to the command line

```
latexmk -gg main.tex
```

To just `make clean` pass `-c` instead of `-gg`.

If it is having any difficulties with compiling, then just `rm -r .latexmk.aux`.
All files there are ephemeral.

If you additionally install `texfot` (included in `texlive-binextra` on
archlinux, not sure about ubuntu), it will parse out only the interesting
outputs from `(pdf|xe|lua)latex`. The `latexmkrc` shipped will automatically
detect this, and only print "interesting" messages.

# Bibliography

`cryptobib` has been included.

The latex source file has been configured, so that Bibliography entries can go
straight into the latex source (look for "bibliography" at the top of the
document, before the document `begin`s). This means all latex source can live
inside this one file, making searching and general text editing more seamless.
