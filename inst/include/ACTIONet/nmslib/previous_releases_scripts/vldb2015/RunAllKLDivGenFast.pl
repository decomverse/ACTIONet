#!/usr/bin/perl -w
require "./RunCommon.pm";

my $DataDir=$ENV{"DATA_DIR"};
defined($DataDir) or die("Define the environemnt variable DATA_DIR, which contains data files.");
my $CacheDir=$ENV{"CACHE_DIR"};
defined($CacheDir) or die("Define the environemnt variable CACHE_DIR, which contains gold standard cache files.");
-d $DataDir or die("DATA_DIR='$DataDir' is not a directory");
my $BinDir=$ENV{"BIN_DIR"};
defined($BinDir) or die("Define the environemnt variable BIN_DIR, which contains binary files.");
-d $BinDir or die("BIN_DIR='$BinDir' is not a directory");
-f "$BinDir/experiment" or die("Missing file '$BinDir/experiment'");


my $SpaceType = "kldivgenfast"; 

my @DataSetName       = ("wikipedia_lda8",      "wikipedia_lda128",
                                                "wikipedia_lda128",
                                                "wikipedia_lda128",
                                                "wikipedia_lda128",
                                                "wikipedia_lda128",
                                                "wikipedia_lda128",
                                                "wikipedia_lda128"
                        );
my @DataSetFile       = ("wikipedia_lda8.txt",  "wikipedia_lda128.txt",
                                                "wikipedia_lda128.txt",
                                                "wikipedia_lda128.txt",
                                                "wikipedia_lda128.txt",
                                                "wikipedia_lda128.txt",
                                                "wikipedia_lda128.txt",
                                                "wikipedia_lda128.txt"
                        );
my @CachePrefix       = ("wikipedia_lda_8_kldiv","wikipedia_lda_128_kldiv",
                                                "wikipedia_lda1M_128_kldiv",
                                                "wikipedia_lda500K_128_kldiv",
                                                "wikipedia_lda250K_128_kldiv",
                                                "wikipedia_lda125K_128_kldiv",
                                                "wikipedia_lda62K_128_kldiv",
                                                "wikipedia_lda31K_128_kldiv"
                        );

my @TuneParamFilePrefix=("tunning/ResultsKLFast/OutFile.wikipedia_lda8", "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                                                "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                                                "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                                                "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                                                "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                                                "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                                                "tunning/ResultsKLFast/OutFile.wikipedia_lda128",
                        );
my @BucketSize        = (50,                     50,
                                                 50,
                                                 50,
                                                 50,
                                                 50,
                                                 50,
                                                 50
                        );
my @MaxNumData        = (2134995,                2134995,
                                                 1000000,
                                                 500000,
                                                 250000,
                                                 125000,
                                                 62000,
                                                 31000
                        );
my @NN                = (7,                   11,
                                              11,
                                              11,
                                              11,
                                              11,
                                              11,
                                              11);

my @IndexAttempts     = (3,                   3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3);
my @SearchAttempts    = ([1..10],             [1..7,10,15,20,25,30],
                                              [1..10],
                                              [1..10],
                                              [1..10],
                                              [1..10],
                                              [1..10],
                                              [1..10]
                        );


my @NumPivotIndex  =(32, 32,
                          32,
                          32,
                          32,
                          32,
                          32,
                          32
);

my @NumPivot       =(2000, 2000,
                            2000,
                            2000,
                            2000,
                            2000,
                            2000,
                            2000
);

my @NumPivotSearch  =([15,20,25,27,29,30,31], [5,8,11,14,17,20,25,30],
                                        [5,8,11,14,17,20],
                                        [5,8,11,14,17,20],
                                        [5,8,11,14,17,20],
                                        [5,8,11,14,17,20],
                                        [5,8,11,14,17,20],
                                        [5,8,11,14,17,20]
);

my @DbScanFrac =([], [],
                     [],
                     [],
                     [],
                     [],
                     [],
                     []
);


my @MaxNumQuery   = (200,                 200,
                                          200,
                                          200,
                                          200,
                                          200,
                                          200,
                                          200);
my @TestSetQty    = (5,                   5,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1);
my @Use           = (1,                   1,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0,
                                          0);


my $DataSetQty = scalar(@DataSetName);

scalar(@DataSetFile)  == $DataSetQty or die("DataSetFile has wrong number of elems");
scalar(@NumPivot)  == $DataSetQty or die("NumPivot has wrong number of elems");
scalar(@NumPivotIndex)  == $DataSetQty or die("NumPivotIndex has wrong number of elems");
scalar(@NumPivotSearch)  == $DataSetQty or die("NumPivotSearch has wrong number of elems");
scalar(@DbScanFrac)  == $DataSetQty or die("DbScanFrac has wrong number of elems");
scalar(@MaxNumData)   == $DataSetQty or die("MaxNumData has wrong number of elems");
scalar(@MaxNumQuery)  == $DataSetQty or die("MaxNumQuery has wrong number of elems");
scalar(@CachePrefix)  == $DataSetQty or die("CachePrefix has wrong number of elems");
scalar(@Use)          == $DataSetQty or die("Use has wrong number of elems");
scalar(@TestSetQty)   == $DataSetQty or die("TestSetQty has wrong number of elems");
scalar(@NN)           == $DataSetQty or die("NN has wrong number of elems");
scalar(@IndexAttempts)== $DataSetQty or die("IndexAttempts has wrong number of elems");
scalar(@SearchAttempts)== $DataSetQty or die("SearchAttempts has wrong number of elems");

#RunTest(1);
RunTest(10);

sub RunTest {
  my ($K) = @_;

  for (my $dn = 0; $dn < $DataSetQty; ++$dn) {
    my $Name = $DataSetName[$dn];
    my $DataQty = $MaxNumData[$dn];

    if (!$Use[$dn]) {
      print "Skipping $Name K=$K (qty=$DataQty)\n";
      next; 
    }

    my $OutFileDir="ResultsKLDivGen/$Name";
    !system("mkdir -p $OutFileDir") or die("Cannot create $OutFileDir");
    my $OutFilePrefix="$OutFileDir/res";

    if (1) {
      RunNAPP($BinDir, $DataDir, $CacheDir,
             $SpaceType,
             $DataSetFile[$dn], 
             $DataQty, 
             $TestSetQty[$dn], 
             $MaxNumQuery[$dn], 
             $CachePrefix[$dn],
             $K, $OutFilePrefix,
             $NumPivot[$dn],
             $NumPivotIndex[$dn],
             $NumPivotSearch[$dn],
             $DbScanFrac[$dn]
            );
    }

    if (1) {
      RunSmallWorld($BinDir, $DataDir, $CacheDir,
             $SpaceType,
             $DataSetFile[$dn], 
             $DataQty, 
             $TestSetQty[$dn], 
             $MaxNumQuery[$dn], 
             $CachePrefix[$dn],
             $K, $OutFilePrefix,
             $NN[$dn],
             $IndexAttempts[$dn],
             $SearchAttempts[$dn]);
    }

    if (1) {
      RunVPTree($BinDir, $DataDir, $CacheDir,
             $SpaceType,
             $DataSetFile[$dn], 
             $DataQty, 
             $TestSetQty[$dn], 
             $MaxNumQuery[$dn], 
             $CachePrefix[$dn],
             $K, $OutFilePrefix,
             $TuneParamFilePrefix[$dn],
             $BucketSize[$dn]);
    }

  }
}

