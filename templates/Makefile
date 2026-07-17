# SPDX-License-Identifier: GPL-2.0-only
#
# Regenerates the fixed-form sources in ../src/legacy from the Fypp
# templates in this directory.  Each template expands to a single .f
# file containing both the single-precision (E) and double-precision
# (F) variant of the routine.
#
# The generated files are committed to the repository, so end users do
# not need Fypp; run this Makefile only after editing a template, and
# commit the regenerated sources together with the template change.
# CI verifies that templates and generated sources have not diverged.

FYPP=fypp
FYPPFLAGS=

SRCDIR=../src/legacy
ROUTINES=y12ma y12mb y12mc y12md y12mg y12mh
GENERATED=$(addprefix $(SRCDIR)/,$(addsuffix .f,$(ROUTINES)))

.PHONY: all
all: $(GENERATED)

$(SRCDIR)/%.f: %.fpp
	$(FYPP) $(FYPPFLAGS) $< $@
