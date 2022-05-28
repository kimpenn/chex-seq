#!/usr/bin/env perl
## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
use strict;
use warnings;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;


my ( $filePrefix, $inFile, $inFh );
my ( $outFile, $outFh );
my @classes;
my $endType = "PE";
my $debug = "info";
my $version = 0;
my $help = 0;

my $headers = "length\tcount\texpect\tmax.err\terror counts";

our $VERSION = "0.22";
our $LOGGER = get_logger(__PACKAGE__);

sub parseArgs {
    GetOptions(
        "filePrefix=s"          => \$filePrefix,
        "endType=s"             => \$endType,
        "outFile=s"             => \$outFile,
        "version|v"             => \$version,
        "debug=s"               => \$debug,
        "help|h"                => \$help,
    ) or usage();
    if ($version) {
        print "$0 v$VERSION\n";
        exit(0);
    }
    if ($help) {
        &usage();
        exit(0);
    }
}

sub usage {
    print STDERR <<DOC
perl $0 --filePrefix data/E.chex/analyzed/Sample_scCLTdegenNuc333/chextrim/Sample_scCLTdegenNuc333.chextrim 
        [--endType PE]
        [--outFile data/E.chex/analyzed/Sample_scCLTdegenNuc333/chextrim/Sample_scCLTdegenNuc333.chextrim.stats.txt]
        [--version|-v]
        [--debug info]
DOC
}

sub setLogger {
    if ( $debug eq "fatal" ) {
        $LOGGER->level($FATAL);
    }
    elsif ( $debug eq "error" ) {
        $LOGGER->level($ERROR);
    }
    elsif ( $debug eq "warn" ) {
        $LOGGER->level($WARN);
    }
    elsif ( $debug eq "info" ) {
        $LOGGER->level($INFO);
    }
    elsif ( $debug eq "debug" ) {
        $LOGGER->level($DEBUG);
    }
    elsif ( $debug eq "trace" ) {
        $LOGGER->level($TRACE);
    }
    my $appender =
      Log::Log4perl::Appender->new("Log::Log4perl::Appender::Screen");
    my $layout = Log::Log4perl::Layout::PatternLayout->new(
        "[%d{yyyy-MM-dd HH:mm:ss.SSS Z} %p] %m");
    $appender->layout($layout);
    $LOGGER->add_appender($appender);
}

sub getClasses {
    my $endType = pop;
    my @classes;
    if ( $endType eq 'SE' ) { 
        @classes = qw( A B BB C CC D E F FF G );
    } elsif ( $endType eq 'PE' ) { 
        @classes = qw( A AA AAX AAA AAAX AAAA 
                       B BB BBX BBB BBBX BBBX 
                       C CX CC CCX CCC 
                       D DX DD DDX DDD 
                       E EE EEX 
                       F FF FFX 
                       G GX 
                       H HX 
                       );
    } else {
        $LOGGER->fatal("endType can be only 'SE' or 'PE'!\n") && die();
    }
    return @classes;

}


&parseArgs(); 
&setLogger();

if (!defined($filePrefix)) {
    $LOGGER->fatal("You didn't provide --filePrefix!\n") && die();
}
    
if (!defined($outFile)) { 
    $outFile = $filePrefix . "stats.txt";
}


$LOGGER->info("{ version = $version, filePrefix = $filePrefix, endType = $endType, outFile = $outFile, debug = $debug, VERSION = $VERSION }\n");

@classes = &getClasses($endType);
open($outFh, ">", $outFile) or $LOGGER->fatal("Cannot open $outFile for output!\n") && die();

print $outFh ( join("\t", qw( Category Length Count Expect MaxError ErrorCounts )), "\n" );

for my $class ( @classes ) {
    if ( $endType eq "SE") {
        $inFile = $filePrefix . "." . $class . ".log";
    } else {
        $inFile = $filePrefix . ".PE" . $class . ".log";
    }

    unless (-f $inFile) {
        $LOGGER->warn("Cannot find info file $inFile for read!\n"); 
    }; 
    open($inFh, "<", $inFile) or $LOGGER->fatal("Cannot open info file $inFile!") && die();
    my $flag = 0;
    my $line;
    
    while( <$inFh> ) {
        chomp;
        if (/$headers/) {
            $flag = 1;
        }
        next unless $flag;
        if ($flag && ($_ ne $headers)) {
            $line = $_;
            $line =~ s/\s+$//;
            $line =~ s/^\s+//;
            my @fields = split("\t", $line);
            if ($#fields == 4) {
                print $outFh (join("\t", ( $class, @fields ) ), "\n"); 
            } else {
                $LOGGER->warn("Why does this line not have 5 columns? $class : $line\n") unless $line eq '';
            }
        }
    }
    close $inFh;
}
$LOGGER->info("All done.\n");
