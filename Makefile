SUBDIRS=src tests matlab python

all: src tests

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS)

clean:
	@for y in $(SUBDIRS); do $(MAKE) -C $$y/ clean ; done

