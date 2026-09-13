#!/bin/bash

mkdir -p cryptobib
curl -#SL https://cryptobib.di.ens.fr/cryptobib/static/files/abbrev3.bib > cryptobib/abbrev3.bib
curl -#SL https://cryptobib.di.ens.fr/cryptobib/static/files/crypto.bib > cryptobib/crypto.bib
