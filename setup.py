import os
from setuptools import setup

try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}

entry_points = {'console_scripts': [
    'aimpoint_mon_update_aimpoint_data=aimpoint_mon.update_aimpoint_data:main',
    'aimpoint_mon_update_observed_aimpoints=aimpoint_mon.update_observed_aimpoints:main',
    'aimpoint_mon_plot_aimpoint=aimpoint_mon.plot_aimpoint:main',
    'aimpoint_mon_make_web_page=aimpoint_mon.make_web_page:main']}

data_files = [(os.path.join('share', 'aimpoint_mon'), ['task_schedule.cfg'])]

setup(name='aimpoint_mon',
      author='Tom Aldcroft',
      description='Chandra aimpoint monitor',
      author_email='taldcroft@cfa.harvard.edu',
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      zip_safe=False,
      entry_points=entry_points,
      packages=['aimpoint_mon'],
      package_data={'aimpoint_mon': ['data/index_template.html']},
      include_package_data=True,
      data_files=data_files,
      tests_require=['pytest'],
      cmdclass=cmdclass,
      )
