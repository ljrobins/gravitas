.PHONY: build

build:
	rm -rf build
	rm -rf setup.py
	poetry version patch
	poetry build


publish:
	poetry publish