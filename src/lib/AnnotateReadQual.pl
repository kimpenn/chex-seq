#!/usr/bin/env perl
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
## Categories:
## A: for PE, per insert, both primers are proper (applicable to paired-end only)
## - A1: fragments long enough to avoid read-through
## - A2: short fragments that have read-through
## B, C:
## - for PE, per insert, only the "focus primer" is proper, and the other primer is not (B: 2p being the focus, C: pC the focus)
## -- B1/C1, long fragments
## -- B2/C2, short fragments
## - for SE, per read, the "focus primer" is proper (B: 2p the focus, C: pC the focus)
## -- B1/C1, long fragments
## -- B2/C2, short fragments
## D:
## - for PE, per insert, neither focus nor non-focus primer is proper
## - for SE, per read, neigther focus nor non-focus primer is proper
##
## Barcode orientation:
## F1R2 or F2R1
## F: the 'focus' barcode, the barcode we are interested in given a class; it is 2p in class A/B, pC in class C;
## F1R2: the focus barcode is found in R1, the non-focus barcode in R2;
## F2R1: the focus barcode is found in R2, the non-focus barcode in R1;

use strict;
use warnings;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;

my ( $inFile,  $inFh );
my ( $outFile, $outFh );
my $endType             = "PE";
my $flags               = "";
my $PrimerLengthMapFile = "src/lib/ChexPrimerTable.conf";
my %PrimerLengthMap;
my ( $primerIdx2p, $primerIdxpC ) = ( "505", "302" );
my $PrimerPairMap = { "2p" => "pC", "pC" => "2p" };

# 505 is the old Chex2ndPrimer (20bp), corresponding to probe 15NTGTG/T;
# 304 is the old Chex1stPrimer (22bp+13Cs);
# 505b is still used in C_V2b;
# old 304 is replaced by 302 (called Chex-App).
my $bcPrimer = "2p";
my $nbPrimer = $PrimerPairMap->{$bcPrimer};
my ( $minL2p, $minLpC ) = ( 6, 6 );
my $version = 0;
my $debug   = "info";
my $help    = 0;
our $VERSION = "0.52";
our $LOGGER  = get_logger(__PACKAGE__);

{

    package SamReader;

    my ( $endType, $rangePrimerL, $bcPrimer, $nbPrimer );

    sub initClass {
        my $class = shift;
        ( $endType, $rangePrimerL, $bcPrimer, $nbPrimer ) = @_;
    }

    sub new {
        my $class = shift;
        my ( $inFile, $outFile, $flags ) = @_;

        my $self = bless {
            inFile  => $inFile,
            outFile => $outFile,
            flags   => $flags,
            inFh    => undef,
            outFh   => undef,
        }, $class;
        return $self;
    }

    sub getTrimStr {
        my $self = shift;
        my ($readName) = @_;
        $readName = ( split( " ", $readName ) )[0];
        my ( $readID, $trimStr ) = split( "::", $readName );
        if ( $endType eq "PE" )
        { # D:0-R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~51,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0
            if ( defined($trimStr) ) {
                $trimStr = ( split( "-", $trimStr ) )[1];
            }
            else {
                ( $readID, undef, $trimStr ) = ( split( "-", $readID ) );
            }
        }
        else {    # D:0,2p_fw~20,pC_rc~17
            my @a = split( ",", $trimStr );
            $trimStr = join( ",", @a[ 1 .. $#a ] );
        }
        return $trimStr;
    }

    sub getPrimerLocs {
        my $self = shift;
        my ($trimStr) = @_;
        my (
            $R1_2p_fw, $R1_pC_fw, $R1_2p_rc, $R1_pC_rc,
            $R2_2p_fw, $R2_pC_fw, $R2_2p_rc, $R2_pC_rc
        ) = ( 0, 0, 0, 0, 0, 0, 0, 0 );
        if ( $endType eq "PE" ) {
            if ( $trimStr =~
/R1_2p_fw~(\d+),R1_pC_fw~(\d+),R1_2p_rc~(\d+),R1_pC_rc~(\d+),R2_2p_fw~(\d+),R2_pC_fw~(\d+),R2_2p_rc~(\d+),R2_pC_rc~(\d+)/
              )
            {
                (
                    $R1_2p_fw, $R1_pC_fw, $R1_2p_rc, $R1_pC_rc,
                    $R2_2p_fw, $R2_pC_fw, $R2_2p_rc, $R2_pC_rc
                ) = ( $1, $2, $3, $4, $5, $6, $7, $8 );
            }
            else {
                $LOGGER->fatal(
"Paired-end read ID containing unexpected CHEX primer locations!\n"
                ) && die();
            }
        }
        else {
            if ( $trimStr =~ /2p_fw~(\d+)/ ) { $R1_2p_fw = $1; }
            if ( $trimStr =~ /2p_rc~(\d+)/ ) { $R1_2p_rc = $1; }
            if ( $trimStr =~ /pC_fw~(\d+)/ ) { $R1_pC_fw = $1; }
            if ( $trimStr =~ /pC_rc~(\d+)/ ) { $R1_pC_rc = $1; }
        }
        my $R1 = {
            "2p" => { "fw" => $R1_2p_fw, "rc" => $R1_2p_rc },
            "pC" => { "fw" => $R1_pC_fw, "rc" => $R1_pC_rc },
        };
        my $R2 = {
            "2p" => { "fw" => $R2_2p_fw, "rc" => $R2_2p_rc },
            "pC" => { "fw" => $R2_pC_fw, "rc" => $R2_pC_rc },
        };
        return { R1 => $R1, R2 => $R2 };
    }

    sub getPrimerScore
    { ## The score is with respect to fwPrimer: fwPrimer has it found in forward strand at the 5'end of the read, while rcPrimer may (read-through) or may not be found at the 3'end of the read.
        my $self = shift;
        my ( $PrimerLocs, $fwPrimer, $rcPrimer ) = @_;
        my ( $minPrimerL, $maxPrimerL ) =
          ( $rangePrimerL->{"min"}, $rangePrimerL->{"max"} );
        my $primerScore;
        if ( $PrimerLocs->{$fwPrimer}->{"fw"} < $minPrimerL->{$fwPrimer} )
        { ## No forward primer found or forward primer found too close to the 5' end;
            $primerScore = 0;
        }
        elsif ( $PrimerLocs->{$fwPrimer}->{"fw"} > $maxPrimerL->{$fwPrimer} )
        {    ## foward primer found but at a position too far away from 5' end;
            $primerScore = 1;
        }
        else {    ## forward primer found at a proper distance
            if ( $PrimerLocs->{$rcPrimer}->{"rc"} > $maxPrimerL->{$rcPrimer} )
            { ## reverse primer's RC found but at a position too far away from the 3' end;
                $primerScore = 2;
            }
            elsif ( $PrimerLocs->{$rcPrimer}->{"rc"} > 0 )
            {    ## reverse primer's RC is found at a proper distance
                $primerScore = 3;
            }
            else
            {    ## reverse primer's RC is not found; this is the best scenario;
                $primerScore = 4;
            }
        }
        return $primerScore;
    }

    sub getPrimerScores {
        my $self = shift;
        my ($PrimerLocs) = @_;
        my $PrimerScores;
        if ( $endType eq "PE" ) {
            my ( $R1, $R2 ) = ( $PrimerLocs->{R1}, $PrimerLocs->{R2} );
            my $R1_bc = $self->getPrimerScore( $R1, $bcPrimer, $nbPrimer );
            my $R1_nb = $self->getPrimerScore( $R1, $nbPrimer, $bcPrimer );
            my $R2_bc = $self->getPrimerScore( $R2, $bcPrimer, $nbPrimer );
            my $R2_nb = $self->getPrimerScore( $R2, $nbPrimer, $bcPrimer );
            $PrimerScores = {
                R1 => { $bcPrimer => $R1_bc, $nbPrimer => $R1_nb },
                R2 => { $bcPrimer => $R2_bc, $nbPrimer => $R2_nb },
            };
        }
        else {
            my $R1    = $PrimerLocs->{R1};
            my $R1_bc = $self->getPrimerScore( $R1, $bcPrimer, $nbPrimer );
            my $R1_nb = $self->getPrimerScore( $R1, $nbPrimer, $bcPrimer );
            $PrimerScores =
              { R1 => { $bcPrimer => $R1_bc, $nbPrimer => $R1_nb }, };
        }
        return $PrimerScores;
    }

    sub getPrimerQual {
        my $self = shift;
        my ($PrimerScores) = @_;
        my ( $qual, $focusTag );
        if ( $endType eq "PE" ) {
            my ( $readInPair, $mateInPair ) = ( "R1", "R2" );
            my $bcPrimerScoreRead = $PrimerScores->{$readInPair}->{$bcPrimer};
            my $nbPrimerScoreMate = $PrimerScores->{$mateInPair}->{$nbPrimer};
            my $bcPrimerScoreMate = $PrimerScores->{$mateInPair}->{$bcPrimer};
            my $nbPrimerScoreRead = $PrimerScores->{$readInPair}->{$nbPrimer};
            my $isA1Read = $bcPrimerScoreRead + $nbPrimerScoreMate == 8;
            my $isA1Mate = $bcPrimerScoreMate + $nbPrimerScoreRead == 8;
            my $isA2Read = ( $bcPrimerScoreRead + $nbPrimerScoreMate == 7 )
              || ( $bcPrimerScoreRead == 3 && $nbPrimerScoreMate == 3 );
            my $isA2Mate = ( $bcPrimerScoreMate + $nbPrimerScoreRead == 7 )
              || ( $bcPrimerScoreMate == 3 && $nbPrimerScoreRead == 3 );
            my $isB1Read = $bcPrimerScoreRead == 4 && $nbPrimerScoreMate <= 2;
            my $isB1Mate = $bcPrimerScoreMate == 4 && $nbPrimerScoreRead <= 2;
            my $isB2Read = $bcPrimerScoreRead == 3 && $nbPrimerScoreMate <= 2;
            my $isB2Mate = $bcPrimerScoreMate == 3 && $nbPrimerScoreRead <= 2;
            my $isC1Read = $bcPrimerScoreMate <= 2 && $nbPrimerScoreRead == 4;
            my $isC1Mate = $bcPrimerScoreRead <= 2 && $nbPrimerScoreMate == 4;
            my $isC2Read = $bcPrimerScoreMate <= 2 && $nbPrimerScoreRead == 3;
            my $isC2Mate = $bcPrimerScoreRead <= 2 && $nbPrimerScoreMate == 3;
            my $isDRead  = $bcPrimerScoreRead <= 2 && $nbPrimerScoreMate <= 2;
            my $isDMate  = $bcPrimerScoreMate <= 2 && $nbPrimerScoreRead <= 2;

            if ( $isA1Read xor $isA1Mate )
            { # It is impossible for both to be true, because CHEXTRIM only searches for bcPrimer in R2 when it is not found in R1. That is, CHEXTRIM won't find bcPrimer in both reads.
                $qual     = "A1";
                $focusTag = $isA1Read ? "F1R2" : "F2R1";
            }
            elsif ( $isA2Read xor $isA2Mate ) {
                $qual     = "A2";
                $focusTag = $isA2Read ? "F1R2" : "F2R1";
            }
            elsif ( $isB1Read xor $isB1Mate ) {
                $qual     = "B1";
                $focusTag = $isB1Read ? "F1R2" : "F2R1";
            }
            elsif ( $isB2Read xor $isB2Mate ) {
                $qual     = "B2";
                $focusTag = $isB2Read ? "F1R2" : "F2R1";
            }
            elsif ( $isC1Read xor $isC1Mate ) {
                $qual     = "C1";
                $focusTag = $isC1Read ? "F1R2" : "F2R1";
            }
            elsif ( $isC2Read xor $isC2Mate ) {
                $qual     = "C2";
                $focusTag = $isC2Read ? "F1R2" : "F2R1";
            }
            elsif ( $isDRead && $isDMate ) {
                $qual     = "D";
                $focusTag = "";
            }
            else {
                $LOGGER->error("Unknown case happened!\n") && die();
            }
        }
        else {
            my $readInPair        = "R1";
            my $bcPrimerScoreRead = $PrimerScores->{$readInPair}->{$bcPrimer};
            my $nbPrimerScoreRead = $PrimerScores->{$readInPair}->{$nbPrimer};
            my $isB1Read = $bcPrimerScoreRead == 4 && $nbPrimerScoreRead <= 2;
            my $isB2Read = $bcPrimerScoreRead == 3 && $nbPrimerScoreRead <= 2;
            my $isC1Read = $bcPrimerScoreRead <= 2 && $nbPrimerScoreRead == 4;
            my $isC2Read = $bcPrimerScoreRead <= 2 && $nbPrimerScoreRead == 3;
            my $isDRead  = $bcPrimerScoreRead <= 2 && $nbPrimerScoreRead <= 2;

            if ($isB1Read) {
                $qual     = "B1";
                $focusTag = "F1";
            }
            elsif ($isB2Read) {
                $qual     = "B2";
                $focusTag = "F1";
            }
            elsif ($isC1Read) {
                $qual     = "C1";
                $focusTag = "F1";
            }
            elsif ($isC2Read) {
                $qual     = "C2";
                $focusTag = "F1";
            }
            elsif ($isDRead) {
                $qual     = "D";
                $focusTag = "";
            }
            else {
                $LOGGER->error("Unknown case happened!\n") && die();
            }
        }
        return ( $qual, $focusTag );
    }

    sub init {
        my $self = shift;
        my ( $inFile,  $inFh )  = ( $self->{inFile},  $self->{inFh} );
        my ( $outFile, $outFh ) = ( $self->{outFile}, $self->{outFh} );
        my $flags = $self->{flags};
        open( $inFh, "samtools view -h $flags $inFile |" )
          or $LOGGER->fatal("Cannot open $inFile for read!\n") && die();
        if ( $outFile =~ /\.sam$/i ) {
            open( $outFh, ">", $outFile )
              or $LOGGER->fatal("Cannot open $outFile for write!\n") && die();
        }
        elsif ( $outFile =~ /\.bam$/i ) {
            open( $outFh, "| samtools view -b > $outFile" )
              or $LOGGER->fatal("Cannot open $outFile for write!\n") && die();
        }
        else {
            $LOGGER->fatal("Unknown output format!\n") && die();
        }
        $self->{inFh}  = $inFh;
        $self->{outFh} = $outFh;
        return $self;
    }

    sub iter {
        my $self = shift;
        my ( $inFh, $outFh ) = ( $self->{inFh}, $self->{outFh} );
        while (<$inFh>) {
            if (/^@/) {
                print $outFh $_;
                next;
            }
            my @f            = split( "\t", $_ );
            my $readName     = $f[0];
            my $rest         = join( "\t", @f[ 1 .. $#f ] );
            my $trimStr      = $self->getTrimStr( $readName, $endType );
            my $PrimerLocs   = $self->getPrimerLocs( $trimStr, $endType );
            my $PrimerScores = $self->getPrimerScores($PrimerLocs);
            my ( $qual, $focusTag ) =
              $self->getPrimerQual( $PrimerScores, $endType );
            print $outFh join( "\t", ( "$readName:$qual:$focusTag", $rest ) );
        }
    }

    sub fin {
        my $self = shift;
        $LOGGER->info("Closing filehandles of input and output...\n");
        close $self->{inFh}
          or $LOGGER->warn("Cannot close filehandle for read!\n");
        close $self->{outFh}
          or $LOGGER->warn("Cannot close filehandle for write!\n");
    }

}

sub loadPrimerLengthMap {
    $PrimerLengthMapFile = $_[0];
    if ( !-e $PrimerLengthMapFile or !-f $PrimerLengthMapFile ) {
        $LOGGER->fatal(
"PrimerLengthMapFile $PrimerLengthMapFile does not exist or is not a file!\n"
        ) && die();
    }
    my %PrimerLengthMap;
    open my $fh, "<", $PrimerLengthMapFile
      or
      $LOGGER->fatal("Cannot open PrimerLengthMapFile $PrimerLengthMapFile!\n")
      && die
      ();
    while (<$fh>) {
        next if /^#/;
        chomp;
        my @f = split( /\s+/, $_ );
        $PrimerLengthMap{ $f[0] } = $f[1];
    }
    close $fh
      or $LOGGER->warn(
"Failed closing file handle for PrimerLengthMapFile $PrimerLengthMapFile!\n"
      );
    %PrimerLengthMap;
}

sub usage {
    print <<DOC;
Summary:
    Read a SAM/BAM file with CHEXTRIM annotated read names, parse the trim length stats and output a corresponding SAM/BAM which has barcode/primer quality information appended to the read names. 

    For example, a record from the input file (PE) is
    NB501328:184:HWHYTBGX5:2:12210:9896:10350-1:N:0:CTTGTA::D:0-R1_2p_fw~20,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~0,R2_2p_fw~0,R2_pC_fw~22,R2_2p_rc~0,R2_pC_rc~0     163     chr1    629222  3     ...

    The output will be
    NB501328:184:HWHYTBGX5:2:12210:9896:10350-1:N:0:CTTGTA::D:0-R1_2p_fw~20,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~0,R2_2p_fw~0,R2_pC_fw~22,R2_2p_rc~0,R2_pC_rc~0:A1:F1R2     163     chr1    629222  3     ...

    'A1' denotes the quality category of this read, and 'F1R2' means the focus primer is in found in R1 while the non-focus primer is found in R2. By "focus primer" we mean the primer on which we focus and evaluate the quality: it is the barcode primer '2p' for class A and B; it is the non-barcode primer 'pC' for class C. It does not apply for D, hence D will always have an empty field.

Usage:
    perl $0 --inFile Sample_XXX.star.posSorted.bam --outFile Sample_XXX.star.posSorted.annotated.bam [--endType PE] [--flags ''] [--PrimerLengthMapFile src/lib/ChexPrimerTable.conf] [--primerIdx2p 505] [--primerIdxpC 302] [--minL2p 6] [--minLpC 6] [--bcPrimer 2p] [--debug info] [--version] [--help]

Options:
    --inFile, -i    input file, BAM/SAM format;
    --outFile, -o   output file BAM/SAM format;
    --flags         additional flags to be passed to samtools, e.g. '-f 0x100' (default: '');
    --endType       choose from 'PE' or 'SE' (default: PE);
    --primerIdx2p   index of the 2p primer (default: 505)
    --primerIdxpC   index of the pC primer (default: 307)
    --minL2p        minimum length (bp) of the 2p primer (inclusive, default: 6)
    --minLpC        minimum length (bp) of the pC primer (inclusive, default: 6)
    --PrimerLengthMapFile   Configurtion file registering each primer's length (default: lib/ChexPrimerTable.conf)
    --bcPrimer      barcode primer, choose from '2p' and 'pC';
    --help, -h      print usage and exit;
    --debug         debug level, choose from fatal, error, warn, info, debug, trace;
    --version, -v   print version and exit;
DOC
}

GetOptions(
    "inFile|i=s"            => \$inFile,
    "outFile|o=s"           => \$outFile,
    "flags=s"               => \$flags,
    "endType=s"             => \$endType,
    "primerIdx2p=s"         => \$primerIdx2p,
    "primerIdxpC=s"         => \$primerIdxpC,
    "minL2p=i"              => \$minL2p,
    "minLpC=i"              => \$minLpC,
    "PrimerLengthMapFile=s" => \$PrimerLengthMapFile,
    "bcPrimer=s"            => \$bcPrimer,
    "help|h"                => \$help,
    "debug=s"               => \$debug,
    "version|v"             => \$version,
) or &usage() && exit(-1);

( &usage() && exit(0) ) if $help;
( ( print "$0 v$VERSION\n" ) && exit(0) ) if $version;
die("--endType must be either 'PE' or 'SE'!\n")
  if $endType ne "PE" && $endType ne "SE";

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
my $appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::Screen");
my $layout   = Log::Log4perl::Layout::PatternLayout->new(
    "[%d{yyyy-MM-dd HH:mm:ss.SSS Z} %p] %m");
$appender->layout($layout);
$LOGGER->add_appender($appender);

## Load primer length database
$LOGGER->fatal("--PrimerLengthMapFile not found!\n") && exit(-1)
  if !defined($PrimerLengthMapFile)
  or !-e $PrimerLengthMapFile;
%PrimerLengthMap = &loadPrimerLengthMap($PrimerLengthMapFile);
for my $k ( sort keys %PrimerLengthMap ) {
    $LOGGER->debug("$k => $PrimerLengthMap{$k}\n");
}

my ( $maxL2p, $maxLpC ) =
  ( $PrimerLengthMap{$primerIdx2p}, $PrimerLengthMap{$primerIdxpC} );
my $minPrimerL   = { "2p"  => $minL2p,     "pC"  => $minLpC };
my $maxPrimerL   = { "2p"  => $maxL2p,     "pC"  => $maxLpC };
my $rangePrimerL = { "min" => $minPrimerL, "max" => $maxPrimerL };

$LOGGER->info(
"{ inFile = $inFile, outFile = $outFile, endType = $endType, flags = $flags, primerIdx2p = $primerIdx2p, primerIdxpC = $primerIdxpC, PrimerLengthMapFile = $PrimerLengthMapFile, minPrimerL = [2p = $minL2p, pC = $minLpC], maxPrimerL = [2p = $maxL2p, pC = $maxLpC], bcPrimer = $bcPrimer, debug = $debug, VERSION = $VERSION }\n"
);

## Initialize class variables
SamReader->initClass( $endType, $rangePrimerL, $bcPrimer, $nbPrimer );
my $samReader = SamReader->new( $inFile, $outFile, $flags );
$LOGGER->info("Opening file handles...\n");
$samReader->init();
$LOGGER->info("Iterating records...\n");
$samReader->iter();
$LOGGER->info("Closing file handles...\n");
$samReader->fin();
$LOGGER->info("All done.\n");
