SUBDIRS = feHMM_1d feHMM_1d_debug feHMM_2d feHMM_2d_debug feHMM_3d feHMM_3d_debug

# EXTRA_DIST does _not_ work properly with subdirs: automake likes to
# copy them twice which cannot work with read-only sources.
#
# EXTRA_DIST = 0d 1d 2d 3d Common
#dist-hook:
#	@echo Copying 0d, 1d, 2d and 3d to $(distdir)
#	for i in 0 1 2 3 ; do \
#		mkdir -p $(distdir)/$${i}d ; \
#		cp -pR $(srcdir)/$${i}d/*.[ch] $(distdir)/$${i}d/ ; \
#	done

#EXTRA_DIST = common
