#!/bin/bash

gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=handouts.pdf handouts1.ps handouts2.ps
