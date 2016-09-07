library("pd.110405.cdiffusa.lz.exp")
library("RUnit")

options(warn=1)

## RUnit Test Suites

dirs <- '.'
testFilePat <- ".*_test\\.R$"

allSuite <- defineTestSuite(name="Test Suite for pd.110405.cdiffusa.lz.exp",
                            dirs=dirs,
                            testFileRegexp=testFilePat,
                            rngKind="default",
                            rngNormalKind="default")

testData <- runTestSuite(allSuite)
printTextProtocol(testData, showDetails=FALSE)

