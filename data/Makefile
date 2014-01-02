.PHONY: clean data data2

data: 
	test -s weather.tar.gz || wget http://www.cl.cam.ac.uk/research/dtg/weather/weather.tar.gz
	tar xzf weather.tar.gz
data2:
	test -s weather-raw.csv || wget http://www.cl.cam.ac.uk/research/dtg/weather/weather-raw.csv

clean:
	rm -f *~
	rm -f weather.tar.gz
	rm -rf daily-text
	rm -f weather-raw.csv
