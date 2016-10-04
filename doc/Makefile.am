docdir = $(datadir)/doc/$(PACKAGE)-$(VERSION)
if DOC_INSTALL
doc_DATA = grid_origin.jpg random_sample.jpg column_major.jpg \
	movie_file_layout.jpg sph_coord.jpg imaging_geometry.jpg \
	coll_half_space.jpg pw_inc.jpg flb_inc.jpg sph_sc.jpg \
	sph_pattern.jpg \
	angora.pdf angora.html angora_*.html angora_*.png
endif

EXTRA_DIST = grid_origin.eps grid_origin.jpg \
	random_sample.eps random_sample.jpg \
	column_major.eps column_major.jpg \
	movie_file_layout.eps movie_file_layout.jpg \
	sph_coord.eps sph_coord.jpg \
	imaging_geometry.eps imaging_geometry.jpg \
	coll_half_space.eps coll_half_space.jpg \
	pw_inc.eps pw_inc.jpg \
	flb_inc.eps flb_inc.jpg \
	sph_sc.eps sph_sc.jpg \
	sph_pattern.eps sph_pattern.jpg


info_TEXINFOS = angora.texi

#AM_MAKEINFOHTMLFLAGS = --no-headers --no-split
#AM_MAKEINFOHTMLFLAGS = --no-split

TEXI2HTML = texi2html
TEXI2HTML_FLAGS = --iftex --l2h --split=section

#html_other = angora_*.html angora_*.png

angora.html: $(info_TEXINFOS) $(html_other)
	$(TEXI2HTML) $(TEXI2HTML_FLAGS) $(srcdir)/angora.texi

clean-local:
	-rm -rf angora_*.* *.pm
