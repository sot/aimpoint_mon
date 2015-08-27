# Set the task name
TASK = fid_drift_mon

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

# Set the names of all files that get installed
SHARE = fid_drift_mon.pl plot_drift.py plot_drift_model.py calc_abs_cel_pointing.py \
	plot_starcheck_vs_telem.py
DATA = task_schedule.cfg
WWW = index.html

include /proj/sot/ska/include/Makefile.FLIGHT

WWW_INSTALL = /data/mta4/www/ASPECT/fid_drift/

install:
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_SHARE)
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
	rsync --times $(WWW) $(WWW_INSTALL)
