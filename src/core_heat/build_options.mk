PWD=$(shell pwd)
EXE_NAME=heat_model
NAMELIST_SUFFIX=heat
override CPPFLAGS += -DCORE_HEAT

report_builds:
	@echo "CORE=heat"
