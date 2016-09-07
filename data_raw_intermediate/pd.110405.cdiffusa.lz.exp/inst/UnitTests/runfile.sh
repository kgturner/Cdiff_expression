#!/bin/bash

TEST_FILE=$1

R --slave <<-EOF
	library('RUnit')
	library('pd.110405.cdiffusa.lz.exp')
	res <- runTestFile('${TEST_FILE}',
        	rngKind='default', rngNormalKind='default')
	printTextProtocol(res, showDetails=FALSE)
EOF
