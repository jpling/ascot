#!gmake

USER=ChristopherWilks
REPO=snaptron-experiments

qs: $(REPO)/qs
	ln -s $<

$(REPO)/qs:
	git clone https://github.com/$(USER)/$(REPO).git -- $(REPO)

.PHONY: clean
clean:
	rm -f qs
	rm -rf $(REPO)
