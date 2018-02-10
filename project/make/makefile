ifeq ($(MESA_DIR),)
ifeq ($($MESA_DIR_INTENTIONALLY_EMPTY),)
  $(error MESA_DIR enviroment variable is not set)
endif
endif




include $(MESA_DIR)/star/work_standard_makefile

# add our MLT module
#Use simply expanded variables '::='
LOAD_STAR_MODS ::= $(LOAD_STAR_MODS) -ltfm_mlt 

STAR = star
