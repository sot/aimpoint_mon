from setuptools import setup

entry_points = {'console_scripts': [
    'aimpoint_mon_update_aimpoint_data=aimpoint_mon.update_aimpoint_data:main',
    'aimpoint_mon_update_observed_aimpoints=aimpoint_mon.update_observed_aimpoints:main',
    'aimpoint_mon_plot_aimpoint=aimpoint_mon.plot_aimpoint:main',
    'aimpoint_mon_make_web_page=aimpoint_mon.make_web_page:main']}

setup(name='aimpoint_mon',
      author='Tom Aldcroft',
      description='Chandra aimpoint monitor',
      author_email='taldcroft@cfa.harvard.edu',
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      zip_safe=False,
      entry_points=entry_points,
      packages=['aimpoint_mon'],
      package_data={'aimpoint_mon': ['data/index_template.html', 'task_schedule.cfg']},
      include_package_data=True,
      tests_require=['pytest'],
      )
