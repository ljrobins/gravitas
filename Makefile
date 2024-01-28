.PHONY: clean test all install docs 

install:
	source bin/activate && pip install -e .

clean:
	rm -rf build
	rm -rf *.so

test:
	pytest tests/*.py

sphinx:
	cd docs && make html

start-runner:
	actions-runner/run.sh

all: clean install test