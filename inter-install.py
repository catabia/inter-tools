#inter-install.py
'''
# -----------------------------------------------------------------------
#
# inter-install.py
#
# by Hannah Catabia
# Last Modified:  WHEN!?!!?
#
# This script installs the Python packages needed to run
# Inter-Tools.
#
# Inter-Tools: A Python toolkit for interactome research
#
# by Hannah Catabia, Caren Smith, and Jose Ordovas
#  
# 
# -----------------------------------------------------------------------
# 
# This program will use pip to install the eight Python 
# packages that are needed in inter-build.py and
# inter-map.py.  It only needs to be run once from
# the command line:
#
# Usage example:
# > python inter-install.py
#
'''

import sys
import pip
from subprocess import check_call, CalledProcessError

def system_check():
	if int(sys.version.split()[0].split('.')[2]) < 9:
		print '\nThis program only runs on Python 2 >= 2.7.9.  Please install and try again.'

	print '\n Checking that the appropriate Python modules are installed...'
	print
	pip.main(['install', 'mygene'])
	print
	pip.main(['install', 'pandas'])
	print
	pip.main(['install', 'rpy2'])
	print
	pip.main(['install', 'matplotlib'])
	print
	pip.main(['install', 'matplotlib_venn'])
	print
	pip.main(['install', 'networkx'])
	print
	pip.main(['install', 'numpy'])
	# checking if R is installed
	try:
		check_call(['which', 'R'])
	except CalledProcessError:
		print 'This program requires R  3.3.3 to run.  Please install R 3.3.3!'
		print
		sys.exit(0)

if __name__ == '__main__':
	system_check()