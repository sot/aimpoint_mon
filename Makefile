# Set the task name
TASK = aimpoint_mon

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

# Set the names of all files that get installed
SHARE = make_web_page.py plot_aimpoint.py update_aimpoint_data.py update_characteristics.py
DATA = index_template.html task_schedule.cfg VERSION

include /proj/sot/ska/include/Makefile.FLIGHT

install:
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_SHARE)
	mkdir -p $(INSTALL_DATA)
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
