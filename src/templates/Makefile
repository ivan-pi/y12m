
FYPP=fypp
FYPPFLAGS=

.PHONY: all
all: legacy/y12mg.f legacy/y12mh.f

legacy/y12mg.f: templates/y12mg.fpp
	$(FYPP) $(FYPPFLAGS) $< $@
legacy/y12mh.f: templates/y12mh.fpp
	$(FYPP) $(FYPPFLAGS) $< $@
