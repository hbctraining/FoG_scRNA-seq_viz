#!/bin/bash

mkdir -p data/
curl -L "https://www.dropbox.com/scl/fo/6n9cfjz0mipjo9o6eqi0x/ALNTS4d-WFAeINFgA5qFlfI?rlkey=blz2flikdexvyk1tsnmd58t1i&st=sahvv3c7&dl=1" -o data/FoGS_scRNAseq.zip
unzip data/FoGS_scRNAseq.zip -d data/