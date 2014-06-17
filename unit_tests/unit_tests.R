kRepoLocation <- Sys.getenv("GIT_REPO_LOC")
library(RUnit)

test.suite <- defineTestSuite("fast_point_mixtures",
                              dirs=file.path(kRepoLocation, "fast_point_mixtures/unit_tests"),
                              testFileRegexp = "^.+runit\\.R",
                              testFuncRegexp = "^test.+",
                              rngKind = "Marsaglia-Multicarry",
                              rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(test.suite)
printTextProtocol(testResult)
