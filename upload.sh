#!/bin/bash

# make sure we have latest handouts
./makeHandouts.sh

# make the scripts zip
zip -r scripts.zip Day_1_scripts Day_2_scripts

# make the whole directory structure thing
rm -rf R_course
mkdir R_course
cp -rf Day_1_scripts R_course
cp -rf Day_2_scripts R_course
cp handouts.pdf R_course
cp Slides_day1.pdf R_course
cp Slides_day2.pdf R_course
tar -pczf Rcourse.tar.gz R_course
rm -rf R_course

# upload everything
scp handouts.pdf rs550@logic.sysbiol.cam.ac.uk:/home/rs550/LOGIC_SITE/teaching/Rcourse/
scp Slides_day1.pdf rs550@logic.sysbiol.cam.ac.uk:/home/rs550/LOGIC_SITE/teaching/Rcourse/
scp Slides_day2.pdf rs550@logic.sysbiol.cam.ac.uk:/home/rs550/LOGIC_SITE/teaching/Rcourse/
scp scripts.zip rs550@logic.sysbiol.cam.ac.uk:/home/rs550/LOGIC_SITE/teaching/Rcourse/
scp Rcourse.tar.gz rs550@logic.sysbiol.cam.ac.uk:/home/rs550/LOGIC_SITE/teaching/Rcourse/
