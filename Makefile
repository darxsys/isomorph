default: all

all: release

docs:
	doxygen Doxyfile

release:
	$(MAKE) -C release

clean:
	$(MAKE) -C release clean

.PHONY: default all docs release clean
