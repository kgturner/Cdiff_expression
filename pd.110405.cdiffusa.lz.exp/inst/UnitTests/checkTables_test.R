pdinfo <- pd.110405.cdiffusa.lz.exp

## if SNP 5.0 or 6.0, genomewide == TRUE
genomewide <- length(grep("CNV", dbListTables(db(pdinfo)))) > 0

test_table_set <- function() {
    ##
    ## this test assumes we know what tables have been defined.  if we
    ## change the schema, this test needs to be updated
    ##

  if (!genomewide){
    alltabs.sort <- sort(c("featureSet", "mmfeature", "pm_mm", "pmfeature",
                           "qcmmfeature", "qcpm_qcmm", "qcpmfeature",
                           "sequence", "sqlite_stat1", "table_info",
                           "fragmentLength"))
  }else{
    alltabs.sort <- sort(c("featureSet", "pmfeature", "sequence",
                            "table_info", "featureSetCNV",
                            "pmfeatureCNV", "sequenceCNV", "sqlite_stat1",
                            "fragmentLength", "fragmentLengthCNV"))
  }
    
    dd  <- pdinfo@getdb()
    ll  <- sort(dbListTables(dd))
    checkEquals(alltabs.sort, ll)
}

test_table_info <- function() {
    ##
    ## this test verifies that the table_info table knows the correct
    ## number of rows for featureSet, and that we can get the featureSet
    ## table into R to count them for 50K chip this is reasonable --
    ## actually a little slow but maybe worth it
    ##
    ## replacing "select *" by "select fsetid" to speed up
  
    dd  <- pdinfo@getdb()
    tinfo  <- dbGetQuery(dd, "select * from table_info")
    nrfs  <- as.numeric(tinfo[ tinfo$tbl == "featureSet", "row_count" ])
    fstab  <- dbGetQuery(dd, "select fsetid from featureSet")
    checkEquals(nrfs, nrow(fstab))

    if (genomewide){
      nrfs  <- as.numeric(tinfo[ tinfo$tbl == "featureSetCNV", "row_count" ])
      fstab  <- dbGetQuery(dd, "select fsetid from featureSetCNV")
      checkEquals(nrfs, nrow(fstab))
    }

  }

test_table_info <- function() {
    ##
    ## this test verifies that the table_info table knows the correct
    ## number of rows for featureSet -- commented out a table extraction
    ##
    dd  <- pdinfo@getdb()
    tinfo  <- dbGetQuery(dd, "select * from table_info")
    nrfs  <- as.numeric(tinfo[ tinfo$tbl == "featureSet", "row_count" ])
#    fstab  <- dbGetQuery(dd, "select * from featureSet")
    nrc  <- dbGetQuery(dd, "select count(*) from featureSet")[[1]]
    checkEquals(nrfs, as.numeric(nrc))
    if (genomewide){
      nrfs  <- as.numeric(tinfo[ tinfo$tbl == "featureSetCNV", "row_count" ])
      fstab  <- dbGetQuery(dd, "select fsetid from featureSetCNV")
      checkEquals(nrfs, nrow(fstab))
    }

}

test_geom <- function() {
    dd <- pdinfo@getdb()
    xylim = pdinfo@geometry
    xs = dbGetQuery(dd, "select x from pmfeature")[[1]]
    ys = dbGetQuery(dd, "select y from pmfeature")[[1]]
    checkTrue(max(xs) <= xylim[2]) # zero offset?
    checkTrue(max(ys) <= xylim[1]) # zero offset?

    if(genomewide){
      xs = dbGetQuery(dd, "select x from pmfeatureCNV")[[1]]
      ys = dbGetQuery(dd, "select y from pmfeatureCNV")[[1]]
      checkTrue(max(xs) <= xylim[2]) # zero offset?
      checkTrue(max(ys) <= xylim[1]) # zero offset?
    }
}

