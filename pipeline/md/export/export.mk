
export_table:
	$(call _start,$(ANCHOR_EXPORT_DIR))
	$(_R) R/export.r export.table \
		ofn=$(EXPORT_TABLE) \
		$(call _export_variable,$(EXPORT_VARIALBES)) \
		$(EXPORT_VARIALBES_NOEVAL)
	$(_end)

export:
	$(call _start,$(EXPORT_ODIR))
	$(_R) R/export.r export \
		odir=$(EXPORT_ODIR) \
		$(call _export_variable,$(EXPORT_VARIALBES)) \
		$(EXPORT_VARIALBES_NOEVAL)
	tar cvf $(EXPORT_ODIR_TAR) $(EXPORT_ODIR)
	$(_end)
