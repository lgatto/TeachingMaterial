.PHONY: clean data

data: 
	test -s weather.tar.gz || wget http://www.cl.cam.ac.uk/research/dtg/weather/weather.tar.gz
	tar xzf weather.tar.gz

clean:
	rm -f *~
	rm -f weather.tar.gz
	rm -rf daily-text
