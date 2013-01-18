#!/bin/sh
# Discover all the tests and create a coverage report.
nosetests  --with-coverage --cover-erase\
           --cover-package=neoart --cover-html
