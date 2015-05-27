SHELL:=/bin/bash

test:
	cram -v tests/*.t

#
# Deploy to the GNU modules system on the PacBio cluster.
# Doesn't update dependencies.
#
DEPLOY_VE:=/mnt/software/d/dumbview/current/VE
deploy:
	(. $(DEPLOY_VE)/bin/activate && python setup.py install)


.PHONY: test deploy
