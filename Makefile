IN  ?= data/input.zip
OUT ?= work

.PHONY: deps example convert sirius all clean

deps:
	pip install -r requirements.txt

example:
	python scripts/make_example_zip.py examples/example_input.zip

convert:
	python scripts/convert_to_ms.py --zip $(IN) --out $(OUT)

sirius:
	chmod +x scripts/run_sirius.sh && ./scripts/run_sirius.sh $(OUT)/ms_files $(OUT)/project

all: convert sirius

clean:
	rm -rf $(OUT)
