# Makefile for handling translation files.
# See README.md for more information.

# Target language for the translation.
LANGUAGE=fr_CA
INSTALL_FILE=$(LANGUAGE)/LC_MESSAGES/fstd2nc.mo

# Rule for compiling the final (binary) translation file.
$(INSTALL_FILE): $(LANGUAGE).po
	mkdir -p $(dir $@)
	msgfmt $< -o $@

# Rule for updating a translation file when new strings are added to the
# program.
$(LANGUAGE).po: fstd2nc.pot
	[ -e $@ ] || msginit -i fstd2nc.pot -l $(LANGUAGE) --no-translator
	msgmerge --update $(LANGUAGE).po fstd2nc.pot
	touch $@

# Rule for generating a template file.
fstd2nc.pot: ../fstd2nc.py
	xgettext $< -d fstd2nc --package-name=fstd2nc --no-location
	mv fstd2nc.po fstd2nc.pot

.PHONY: clean

# Rule for cleaning up the files.
# Removes the auto-generated template and binary files, but leaves the
# translation file (.po).
clean:
	rm -f fstd2nc.pot $(INSTALL_FILE)

