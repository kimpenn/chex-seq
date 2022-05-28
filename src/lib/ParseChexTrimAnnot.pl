#!/usr/bin/env perl
## Youtao Lu <luyoutao@sas.upenn.edu>
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
## - for SE, per read, the "focus primer" is proper (B: 2p being the focus, C: pC the focus)
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

my ( $inFile1, $inFile2 );
my $outFile;
my $endType             = "PE";
my $PrimerLengthMapFile = "src/lib/PrimerLengthTable.conf";
my %PrimerLengthMap;
my ( $primerIdx2p, $primerIdxpC ) = ( "505", "302" );

## 505 is the old Chex2ndPrimer (20bp), corresponding to probe 15NTGTG/T;
## 304 is the old Chex1stPrimer (22bp+13Cs)
## 505b is still used in C_V2b; old 304 is replaced by 302 (called Chex-App).
my $PrimerPairMap = { "2p" => "pC", "pC" => "2p" };
my $bcPrimer      = "2p";
my $nbPrimer      = $PrimerPairMap->{$bcPrimer};
my ( $minL2p, $minLpC ) = ( 6, 6 );
my $version = 0;
my $help    = 0;
my $debug   = "info";
our $VERSION = "0.43";
our $LOGGER  = get_logger(__PACKAGE__);

{

    package Read;
    my ( $endType, $rangePrimerL );

    sub initClass {
        my $class = shift;
        ( $endType, $rangePrimerL ) = @_;
    }

    sub new {
        my $class = shift;
        my ( $readName, $seq, $sep, $phred, $readInPair ) = @_;
        my $self = bless {
            readName     => $readName,
            seq          => $seq,
            sep          => $sep,
            phred        => $phred,
            readInPair   => $readInPair,
            readID       => undef,
            trimStr      => undef,
            yI           => undef,
            yD           => undef,
            PrimerLocs   => undef,
            PrimerScores => undef,
        }, $class;
        return $self;
    }

    sub parseReadName {
        my $self       = shift;
        my $readInPair = $self->{readInPair};
        my $readName   = $self->{readName};
        my ( $readID, $trimStr, $yI, $yD ) = (undef, undef, "", "");
        if ( $readInPair eq "R1" ) {

# @NB501328:212:HFW2MBGX9:2:22305:8919:10250-1:N:0:ACTGAT::D:0-R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~51,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0 L:53
            ( $readName, undef ) = split( " ", $readName );

# @NB501328:212:HFW2MBGX9:2:22305:8919:10250-1:N:0:ACTGAT::D:0-R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~51,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0
            ( $readID, $trimStr ) = split( "::", $readName );

# @NB501328:212:HFW2MBGX9:2:22305:8919:10250-1:N:0:ACTGAT
# D:0-R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~51,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0
            if (defined($trimStr)) {
                ( $readID, $yI ) = split( "-", $readID );
            # @NB501328:212:HFW2MBGX9:2:22305:8919:10250
                ( $yD, $trimStr ) = split( "-", $trimStr );
                $yD = $1 if $yD =~ /D:([0-9]+)/;
            } else {
                ( $readID, $yI, $trimStr ) = split( "-", $readID );
            }
            $yI =~ /1:N:[0-9]+:([ACGTN]+)/;
            $yI = $1;
            $LOGGER->trace("\$readID = $readID, \$yI = $yI, \$trimStr = $trimStr\n");
# R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~51,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0
        }
        elsif ( $readInPair eq "R2" ) {

# @NB501328:212:HFW2MBGX9:2:22305:8919:10250 2:N:0:ACTGAT::D:0 R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0 L:53
            ( $readID, undef, undef, undef ) = split( " ", $readName );
        }
        else {
            $LOGGER->fatal("Something other than R1 or R2 found!\n") && die();
        }
        $self->{readID}  = substr( $readID, 1 );
        $self->{trimStr} = $trimStr;
        $self->{yI}      = $yI;
        $self->{yD}      = $yD;
    }

    sub setPrimerScore
    { ## The score is with respect to fwPrimer: fwPrimer has it found in forward strand at the 5'end of the read, while rcPrimer may (read-through) or may not be found at the 3'end of the read.
        my $self = shift;
        my ( $fwPrimer, $rcPrimer ) = @_;
        my $PrimerLocs = $self->{PrimerLocs};
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
        $self->{PrimerScores}->{$fwPrimer} = $primerScore;
        $LOGGER->trace(
"\$self->{PrimerScores}->{$fwPrimer} = $self->{PrimerScores}->{$fwPrimer}\n"
        );
    }

    sub setPrimerScores {
        my $self = shift;
        $self->setPrimerScore( "2p", "pC" );
        $self->setPrimerScore( "pC", "2p" );
        $LOGGER->trace(
            "\$self->{PrimerScores}->{'2p'} = $self->{PrimerScores}->{'2p'}\n");
        $LOGGER->trace(
            "\$self->{PrimerScores}->{pC} = $self->{PrimerScores}->{pC}\n");
    }
}

{

    package ReadPair;
    use IO::Zlib;
    use File::Path qw( make_path );

## we treat barcode/primer information as class variables so that we don't set them up upon every object creation
    my ( $endType, $rangePrimerL, $bcPrimer, $nbPrimer );

    sub initClass {
        my $class = shift;
        ( $endType, $rangePrimerL, $bcPrimer, $nbPrimer ) = @_;
    }

    sub new {
        my $class = shift;
        my ( $inFile1, $inFile2, $outFile ) = @_;
        Read->initClass( $endType, $rangePrimerL );

        my $self = bless {
            inFile1  => $inFile1,
            inFile2  => $inFile2,
            outFile  => $outFile,
            inFh1    => undef,
            inFh2    => undef,
            outFh    => undef,
            R1       => undef,
            R2       => undef,
            qual     => undef,
            focusTag => undef,
        }, $class;
        return $self;
    }

    sub initFH {
        my $self = shift;
        my ( $inFile1, $inFh1 ) = ( $self->{inFile1}, undef );
        my ( $inFile2, $inFh2 ) = ( $self->{inFile2}, undef );
        my ( $outFile, $outFh ) = ( $self->{outFile}, undef );

        $LOGGER->info("Initializing filehandles for $inFile1...\n");
        if ( $inFile1 =~ /(gz|gzip)$/i ) {
            $inFh1 = IO::Zlib->new( $inFile1, "r" )
              or $LOGGER->fatal("Cannot open $inFile1 for read!\n") && die();
        }
        else {
            open( $inFh1, '<', $inFile1 )
              or $LOGGER->fatal("Cannot open $inFile1 for read!\n") && die();
        }

        $LOGGER->info( "Initializing filehandles for "
              . ( defined($outFile) ? $outFile : "STDOUT" )
              . "...\n" );
        if ( !defined($outFile) ) {
            $outFh = *STDOUT;
        }
        elsif ( $outFile =~ /(gz|gzip)$/i ) {
            $outFh = IO::Zlib->new( $outFile, "w" )
              or $LOGGER->fatal("Cannot open $outFile for write!\n") && die();
        }
        else {
            open $outFh, '>', $outFile
              or $LOGGER->fatal("Cannot open $outFile for write!\n") && die();
        }

        if ( $endType eq "PE" ) {
            $LOGGER->info("Initializing filehandles for $inFile2...\n");
            if ( $inFile2 =~ /(gz|gzip)$/i ) {
                $inFh2 = IO::Zlib->new( $inFile2, "r" )
                  or $LOGGER->fatal("Cannot open $inFile2 for read!\n")
                  && die();
            }
            else {
                open( $inFh2, '<', $inFile2 )
                  or $LOGGER->fatal("Cannot open $inFile2 for read!\n")
                  && die();
            }
        }
        else {
            if ( defined($inFile2) && -e $inFile2 ) {
                $LOGGER->warn(
"You provided $inFile2 as R2 but also specified endType = $endType, ignoring it!\n"
                );
            }
        }
        print $outFh (
            join(
                "\t", qw(
                  ReadID Index Rmdup
                  R1_2p_fw R1_2p_rc R1_pC_fw R1_pC_rc
                  R2_2p_fw R2_2p_rc R2_pC_fw R2_pC_rc
                  Qual Strand )
            ),
            "\n"
        );
        $self->{inFh1} = $inFh1;
        $self->{inFh2} = $inFh2;
        $self->{outFh} = $outFh;
        return $self;
    }

    sub nextReadPair {
        my $hasNext = 1;
        my $self    = shift;
        my $inFh1   = $self->{inFh1};
        my $inFh2   = $self->{inFh2};
        my $rname1  = readline($inFh1);
        if ( !defined($rname1) ) { $hasNext = 0; }
        if ($hasNext) {
            readline($inFh1);
            readline($inFh1);
            readline($inFh1);
            $self->{R1} = Read->new( $rname1, undef, undef, undef, "R1" );
            if ( $endType eq "PE" ) {
                my $rname2 = readline($inFh2);
                readline($inFh2);
                readline($inFh2);
                readline($inFh2);
                $self->{R2} = Read->new( $rname2, undef, undef, undef, "R2" );
            }
        }
        return $hasNext;
    }

    sub parseReadNames {
        my $self = shift;
        $self->{R1}->parseReadName();
        if ( $endType eq "PE" ) {
            $self->{R2}->parseReadName();
            if ( $self->{R1}->{readID} ne $self->{R2}->{readID} ) {
                $LOGGER->fatal("Your R1 and R2 have different read IDs!\n")
                  && die();
            }
        }
    }

    sub setPrimerLocs {
        my $self = shift;
        my ( $R1, $R2 ) = ( $self->{R1}, $self->{R2} );
        my $trimStr1 =
          $R1->{trimStr};    # ChexTrim stats is only fully preserved in R1
        my (
            $R1_2p_fw, $R1_pC_fw, $R1_2p_rc, $R1_pC_rc,
            $R2_2p_fw, $R2_pC_fw, $R2_2p_rc, $R2_pC_rc
        );
        my ( $PrimerLocs1, $PrimerLocs2 );
        if ( $endType eq "PE" ) {
            if ( $trimStr1 =~
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
            if ( $trimStr1 =~ /2p_fw~(\d+)/ ) { $R1_2p_fw = $1; }
            if ( $trimStr1 =~ /2p_rc~(\d+)/ ) { $R1_2p_rc = $1; }
            if ( $trimStr1 =~ /pC_fw~(\d+)/ ) { $R1_pC_fw = $1; }
            if ( $trimStr1 =~ /pC_rc~(\d+)/ ) { $R1_pC_rc = $1; }
            ( $R2_2p_fw, $R2_2p_rc, $R2_pC_fw, $R2_pC_rc ) = ( "", "", "", "" );
        }
        $PrimerLocs1 = {
            "2p" => { "fw" => $R1_2p_fw, "rc" => $R1_2p_rc },
            "pC" => { "fw" => $R1_pC_fw, "rc" => $R1_pC_rc },
        };
        $PrimerLocs2 = {
            "2p" => { "fw" => $R2_2p_fw, "rc" => $R2_2p_rc },
            "pC" => { "fw" => $R2_pC_fw, "rc" => $R2_pC_rc },
        };
        $self->{R1}->{PrimerLocs} = $PrimerLocs1;
        $self->{R2}->{PrimerLocs} = $PrimerLocs2;
    }

    sub setPrimerScores {
        my $self = shift;
        $self->{R1}->setPrimerScores();
        if ( $endType eq "PE" ) {
            $self->{R2}->setPrimerScores();
        }
    }

    sub setPairedEndQual {
        my $self = shift;
        my ( $readInPair, $mateInPair ) = ( "R1", "R2" );
        my $bcPrimerScoreRead =
          $self->{$readInPair}->{PrimerScores}->{$bcPrimer};
        my $nbPrimerScoreMate =
          $self->{$mateInPair}->{PrimerScores}->{$nbPrimer};
        my $bcPrimerScoreMate =
          $self->{$mateInPair}->{PrimerScores}->{$bcPrimer};
        my $nbPrimerScoreRead =
          $self->{$readInPair}->{PrimerScores}->{$nbPrimer};
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
            $self->{qual}     = "A1";
            $self->{focusTag} = $isA1Read ? "F1R2" : "F2R1";
        }
        elsif ( $isA2Read xor $isA2Mate ) {
            $self->{qual}     = "A2";
            $self->{focusTag} = $isA2Read ? "F1R2" : "F2R1";
        }
        elsif ( $isB1Read xor $isB1Mate ) {
            $self->{qual}     = "B1";
            $self->{focusTag} = $isB1Read ? "F1R2" : "F2R1";
        }
        elsif ( $isB2Read xor $isB2Mate ) {
            $self->{qual}     = "B2";
            $self->{focusTag} = $isB2Read ? "F1R2" : "F2R1";
        }
        elsif ( $isC1Read xor $isC1Mate ) {
            $self->{qual}     = "C1";
            $self->{focusTag} = $isC1Read ? "F1R2" : "F2R1";
        }
        elsif ( $isC2Read xor $isC2Mate ) {
            $self->{qual}     = "C2";
            $self->{focusTag} = $isC2Read ? "F1R2" : "F2R1";
        }
        elsif ( $isDRead && $isDMate ) {
            $self->{qual}     = "D";
            $self->{focusTag} = "";
        }
        else {
            $LOGGER->error("We don't really expect this case!\n") && die();
        }
    }

    sub setSingleEndQual {
        my $self              = shift;
        my $bcPrimerScoreRead = $self->{R1}->{PrimerScores}->{$bcPrimer};
        my $nbPrimerScoreRead = $self->{R1}->{PrimerScores}->{$nbPrimer};
        $LOGGER->trace(
"\$self->{R1}->{PrimerScores}->{$bcPrimer} = $self->{R1}->{PrimerScores}->{$bcPrimer}\n"
        );
        $LOGGER->trace(
"\$self->{R1}->{PrimerScores}->{$nbPrimer} = $self->{R1}->{PrimerScores}->{$nbPrimer}\n"
        );

        my $isB1Read = $bcPrimerScoreRead == 4 && $nbPrimerScoreRead <= 2;
        my $isB2Read = $bcPrimerScoreRead == 3 && $nbPrimerScoreRead <= 2;
        my $isC1Read = $bcPrimerScoreRead <= 2 && $nbPrimerScoreRead == 4;
        my $isC2Read = $bcPrimerScoreRead <= 2 && $nbPrimerScoreRead == 3;
        my $isDRead  = $bcPrimerScoreRead <= 2 && $nbPrimerScoreRead <= 2;

        if ($isB1Read) {
            $self->{qual}     = "B1";
            $self->{focusTag} = "F1";
        }
        elsif ($isB2Read) {
            $self->{qual}     = "B2";
            $self->{focusTag} = "F1";
        }
        elsif ($isC1Read) {
            $self->{qual}     = "C1";
            $self->{focusTag} = "F1";
        }
        elsif ($isC2Read) {
            $self->{qual}     = "C2";
            $self->{focusTag} = "F1";
        }
        elsif ($isDRead) {
            $self->{qual}     = "D";
            $self->{focusTag} = "";
        }
        else {
            $LOGGER->error("This case should never happen!\n") && die();
        }
    }

    sub setPrimerQuals {
        my $self = shift;
        if ( $endType eq "PE" ) {
            $self->setPairedEndQual();
        }
        else {
            $self->setSingleEndQual();
        }
    }

    sub output {
        my $self   = shift;
        my $outFh  = $self->{outFh};
        my $readID = $self->{R1}->{readID};
        my ( $yI, $yD ) = ( $self->{R1}->{yI}, $self->{R1}->{yD} );
        my $yP = $self->{R1}->{trimStr};
        my ( $qual, $focusTag ) = ( $self->{qual}, $self->{focusTag} );

        my $PrimerLocs1 = $self->{R1}->{PrimerLocs};
        my $PrimerLocs2 = $self->{R2}->{PrimerLocs};
        my @PrimerLocs  = (
            $PrimerLocs1->{'2p'}->{fw}, $PrimerLocs1->{'2p'}->{rc},
            $PrimerLocs1->{'pC'}->{fw}, $PrimerLocs1->{'pC'}->{rc},
            $PrimerLocs2->{'2p'}->{fw}, $PrimerLocs2->{'2p'}->{rc},
            $PrimerLocs2->{'pC'}->{fw}, $PrimerLocs2->{'pC'}->{rc},
        );
        $LOGGER->trace(
"\$yP = $yP: [$PrimerLocs1->{'2p'}->{fw}, $PrimerLocs1->{'pC'}->{fw}, $PrimerLocs1->{'2p'}->{rc}, $PrimerLocs1->{'pC'}->{rc}, $PrimerLocs2->{'2p'}->{fw}, $PrimerLocs2->{'pC'}->{fw}, $PrimerLocs2->{'2p'}->{rc}, $PrimerLocs2->{'pC'}->{rc}]\n"
        );
        print $outFh (
            join( "\t", ( $readID, $yI, $yD, @PrimerLocs, $qual, $focusTag ) ),
            "\n"
        );
    }

    sub finFH {
        my $self = shift;
        my ( $inFile1, $inFile2, $outFile ) =
          ( $self->{inFile1}, $self->{inFile2}, $self->{outFile} );
        my ( $inFh1, $inFh2, $outFh ) =
          ( $self->{inFh1}, $self->{inFh2}, $self->{outFh} );
        $LOGGER->info("Closing filehandles of $inFile1...\n");
        close $self->{inFh1}
          or $LOGGER->warn("Cannot close filehandle of $inFile1 for read!\n");
        $LOGGER->info( "Closing filehandles of "
              . ( defined($outFile) ? $outFile : "STDOUT" )
              . "...\n" );
        close $outFh
          or $LOGGER->warn( "Cannot close filehandle of "
              . ( defined($outFile) ? $outFile : "STDOUT" )
              . " for write!\n" );
        if ( $endType eq "PE" ) {
            $LOGGER->info("Closing filehandles of $inFile2...\n");
            close $self->{inFh2}
              or
              $LOGGER->warn("Cannot close filehandle of $inFile2 for read!\n");
        }
    }
}

sub parseArgs {
    GetOptions(
        "inFile1=s"             => \$inFile1,
        "inFile2=s"             => \$inFile2,
        "outFile|o=s"           => \$outFile,
        "endType=s"             => \$endType,
        "PrimerLengthMapFile=s" => \$PrimerLengthMapFile,
        "primerIdx2p=s"         => \$primerIdx2p,
        "primerIdxpC=s"         => \$primerIdxpC,
        "minL2p=i"              => \$minL2p,
        "minLpC=i"              => \$minLpC,
        "bcPrimer=s"            => \$bcPrimer,
        "version|v"             => \$version,
        "help|h"                => \$help,
        "debug=s"               => \$debug,
    ) or usage() && exit(-1);
}

sub usage {
    print <<DOC;
Summary:
    Parse the read names from FASTQ files annotated by CHEXTRIM module.

Usage:
    perl $0 --inFile1 unaligned_1.fq.gz [--inFile2 unaligned_2.fq.gz] --outFile out.tsv.gz [--endType PE] [--PrimerLengthMapFile src/lib/ChexPrimerTable.conf] [--primerIdx2p 505] [--primerIdxpC 307] [--minL2p 6] [--minLpC 6] [-bcPrimer 2p] [--debug info]

Options:
    --inFile1       input FASTQ file for R1, can be plain or GZIP (.gz);
    --inFile2       input FASTQ file for R2, can be plain or GZIP (.gz); 
                    must be present if --endType is 'PE';
    -o, --outFile   output file, tab separated, can be plain or GZIP (.gz);
    --endType       choose from 'PE' or 'SE' (default: 'PE');
    --primerIdx2p   index of the 2p primer (default: 505);
    --primerIdxpC   index of the polyC primer (default: 302);
    --PrimerLengthMapFile   Config file for primers' length (default: src/lib/ChexPrimerTable.conf);
    --minL2p        minimum length (bp) of the 2p primer to be found (inclusive, default: 6);
    --minLpC        minimum length (bp) of the polyC primer to be found (inclusive, default: 6);
    --debug         log level, choose from 'fatal', 'error', 'warn', 'info', 'debug', 'trace' (default: info);

Output:
    A typical read name from CHEXTRIM looks like:
    \@NB501328:229:HGMMYBGXB:3:12601:3543:11957-1:N:0:AGTCAA::D:0-R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~53,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0

    A typical read name from TRIM looks like:
    \@NB501328:229:HGMMYBGXB:3:12601:3543:11957-1:N:0:AGTCAA::D:0-R1_2p_fw~21,R1_pC_fw~0,R1_2p_rc~0,R1_pC_rc~53,R2_2p_fw~0,R2_pC_fw~23,R2_2p_rc~0,R2_pC_rc~0 L:53

    This program should tolerate the above variants and output 13 columms:
    1) read ID (ReadID), NB501328:229:HGMMYBGXB:3:12601:3543:11957;
    2) Illumina sequencing index (Index), AGTCAA;
    3) RMDUP duplicate copy number (Rmdup), 0;
    4) CHEXTRIM primer trimmed length (R1_2p_fw), 21; 
    5) CHEXTRIM primer trimmed length (R1_2p_rc), 0; 
    6) CHEXTRIM primer trimmed length (R1_pC_fw), 0; 
    7) CHEXTRIM primer trimmed length (R1_pC_rc), 53; 
    8) CHEXTRIM primer trimmed length (R2_2p_fw), 0; 
    9) CHEXTRIM primer trimmed length (R2_2p_rc), 0; 
    10) CHEXTRIM primer trimmed length (R2_pC_fw), 23; 
    11) CHEXTRIM primer trimmed length (R2_pC_rc), 0; 
    12) CHEXTRIM barcode/primer primer quality category (Qual), C1;
    13) CHEXTRIM barcode/primer strand orientation (Strand), F2R1; 

Example:
    perl $0 --inFile1 data/E.chex/analyzed/Sample_scCLTdegenNuc333/chextrim/unaligned_1.fq.gz --inFile2 data/E.chex/analyzed/Sample_scCLTdegenNuc333/chextrim/unaligned_2.fq.gz 2>/dev/null
ReadID  Index   Rmdup   R1_2p_fw        R1_2p_rc        R1_pC_fw        R1_pC_rc        R2_2p_fw        R2_2p_rc        R2_pC_fw        R2_pC_rc        Qual    Strand
NB501328:229:HGMMYBGXB:3:12601:3543:11957       AGTCAA  0       21      0       0       53      0       0       23      0       C1      F2R1
NB501328:229:HGMMYBGXB:1:23206:25951:4173       AGTCAA  0       20      0       0       6       0       0       24      0       A2      F1R2
NB501328:229:HGMMYBGXB:1:12112:15752:4296       AGTCAA  0       21      0       0       6       0       0       42      0       D
NB501328:229:HGMMYBGXB:1:22107:20741:16647      AGTCAA  0       21      0       0       49      0       0       23      0       C1      F2R1
NB501328:229:HGMMYBGXB:3:12410:15209:16644      AGTCAA  0       18      0       0       43      0       0       23      0       C1      F2R1
...
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

&parseArgs();
( print "$0 v$VERSION\n" ) && exit(0) if $version;
&usage() && exit(0) if $help;
&setLogger();
$LOGGER->fatal("R1 input is not provided!\n") && die(1) if !defined($inFile1);
$LOGGER->fatal("Cannot find $inFile1!\n") && die(1) if !-e $inFile1;
$LOGGER->fatal("R2 input is not provided but endType is set to 'PE'!\n")
  && die(2)
  if $endType eq "PE"
  && !defined($inFile2);
$LOGGER->fatal("Cannot find $inFile2 but endType is set to 'PE'!\n") && die(2)
  if $endType eq "PE" && !-e $inFile2;

## Load primer length database
%PrimerLengthMap = &loadPrimerLengthMap($PrimerLengthMapFile);

for my $k ( sort keys %PrimerLengthMap ) {
    $LOGGER->debug("$k => $PrimerLengthMap{$k}\n");
}

my ( $maxL2p, $maxLpC ) =
  ( $PrimerLengthMap{$primerIdx2p}, $PrimerLengthMap{$primerIdxpC} );
my $minPrimerL   = { "2p"  => $minL2p,     "pC"  => $minLpC };
my $maxPrimerL   = { "2p"  => $maxL2p,     "pC"  => $maxLpC };
my $rangePrimerL = { "min" => $minPrimerL, "max" => $maxPrimerL };

$LOGGER->info( "{ inFile1 = $inFile1, inFile2 = "
      . ( defined($inFile2) ? $inFile2 : "" )
      . ", outFile = "
      . ( defined($outFile) ? $outFile : "" )
      . ", endType = $endType, primerIdx2p = $primerIdx2p, primerIdxpC = $primerIdxpC, PrimerLengthMapFile = $PrimerLengthMapFile, minPrimerL = (2p = $minL2p, pC = $minLpC), maxPrimerL = (2p = $maxL2p, pC = $maxLpC), bcPrimer = $bcPrimer, nbPrimer = $nbPrimer, debug = $debug, VERSION = $VERSION }\n"
);

## Initialize class variables
$LOGGER->info("Initializing class ...\n");
ReadPair->initClass( $endType, $rangePrimerL, $bcPrimer, $nbPrimer );
my $ReadPair = ReadPair->new( $inFile1, $inFile2, $outFile );
$LOGGER->info("Opening file handels ...\n");
$ReadPair->initFH();
$LOGGER->info(
    $endType eq "PE"
    ? "Iterating read pairs ...\n"
    : "Iterating reads...\n"
);
my $i = 0;
while ( $ReadPair->nextReadPair() ) {
    $i++;
    $ReadPair->parseReadNames();
    $ReadPair->setPrimerLocs();
    $ReadPair->setPrimerScores();
    $ReadPair->setPrimerQuals();
    $ReadPair->output();
    $LOGGER->info("Processed $i\n") if $i % 10000 == 0;
}
$LOGGER->info("Closing file handels ...\n");
$ReadPair->finFH();
$LOGGER->info("Total $i\n");
$LOGGER->info("All done.\n");
