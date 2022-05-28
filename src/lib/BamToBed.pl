#!/usr/bin/env perl
## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
## Turn a adjacent duplicates removed bam into a bed file
## The point here is to keep the strand of the fragment to same as the focus read
## in a read pair.
## 
## Options:
## 1. whether to output read pairs that have no focus ('F') read;
## 2. output types: 1) 5End, 2) Read 3) Frag, 4) ReadAndFrag
## 
## Example:
## Given PE data
## Input:read pair F1------>     <----------R2                    F1------>     <----------R2
##       singleton        F------>                                       F------>
##       singleton          <-------R                                      <-------R
## Output:
##                    doesOutNoFocus=Y                                doesOutNoFocus=N
##
## 5End           [F]------>     <---------[R]                   [F]------->
##                       [F]------>                                     [F]------>
##                          <------[R]
##
## Read           [F------->]   [<----------R]                   [F-------->]
##                       [F------->]                                    [F------->]
##                          [<-------R]
##
## Frag           [F----------------------->R] (paired only)     [F---------------------->R] (paired only)
##
## ReadAndFrag    [F----------------------->R] (paired or        [F---------------------->R] (paired or
##                       [F------->] singleton)                         [F-------->]  singleton)
##                          [<-------R]
##
## Given SE data
## Input:          F------->                                      F------->
##                               <-----------                                  <-----------
## Output:
##                    doesOutNoFocus=Y                                doesOutNoFocus=N
## 5End           [F]------>     <---------[-]                   [F]------->
## Read           [F------->]   [<-----------]                   [F-------->]
## Frag:          Not applicable
## ReadAndFrag:   Not applicable
############################################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use POSIX;
use File::Basename;
use File::Path;
use File::Copy;

my ( $inFile, $outFile );
my $isNameSorted   = 0;
my $endType        = "PE";
my $outType        = "5End";
my $samFlags       = "";
my $focusTag       = "F";
my $doesOutNoFocus = 0;
my $addSoftClipped = 0;
my $samSortNCores  = 1;
our $VERSION = "0.71";
my $verbose = 0;
my $debug   = 0;

sub usage {
    print STDERR
"perl $0 --inFile inFile.bam --outFile outFile.bed [--isNameSorted] [--endType PE] [--outType '5End'|'Read'|'Frag'|'ReadAndFrag'] [--addSoftClipped] [--focusTag F] [--doesOutNoFocus] [--samFlags ''] [--samSortNCores 1] [--verbose] [--debug]
    Options:
    
    --inFile        path of the input bam file 
    --outFile       path of the output bed file
    --isNameSorted  [optional] whether the input bam is sorted by name (default: false)
    --endType       [optional] input file being PE or SE sequencing (choose from 'PE', 'SE'; default: 'PE')
    --outType       [optional] which range to output to the bed file (choose from '5End', 'Read', 'Frag', 'ReadAndFrag'; default: '5End')
    --samFlags      [optional] additional arguments for `samtools view ` (default: '')
    --focusTag      [optional] if the bam file is tagged with focusTag, what is it? (default 'F')
    --doesOutNoFocus [optional] whether or not include the location of NoFocus reads? (default false)
    --addSoftClipped [optional] if read has focusTag found, whether to accommodate soft-clipped length to infer focusTag location (default: false)
    --samSortNCores [optional] how many threads should `samtools sort` use? (default: 1)
    --verbose       [optional] whether to output running information
    --debug         [optional] whether to output more detailed information
";
    exit(-1);
}

GetOptions(
    "inFile=s"        => \$inFile,
    "outFile=s"       => \$outFile,
    "isNameSorted"    => \$isNameSorted,
    "endType=s"       => \$endType,
    "outType=s"       => \$outType,
    "addSoftClipped"  => \$addSoftClipped,
    "focusTag=s"      => \$focusTag,
    "doesOutNoFocus"  => \$doesOutNoFocus,
    "samFlags=s"      => \$samFlags,
    "samSortNCores=i" => \$samSortNCores,
    "verbose"         => \$verbose,
    "debug"           => \$debug,
) or &usage();

if ($verbose) {
    print STDERR "[Parameter settings:]
    inFile = $inFile
    outFile = $outFile
    isNameSorted = $isNameSorted
    endType = $endType
    outType = $outType
    addSoftClipped = $addSoftClipped
    doesOutNoFocus = $doesOutNoFocus
    focusTag = $focusTag
    samFlags = $samFlags
    samSortNCores = $samSortNCores
    VERSION = $VERSION
    debug = $debug
";
}

{

    package Read;
    use Carp qw /croak carp/;

    my ( $endType, $focusTag, $addSoftClipped ) =
      ( "PE", "F", 0 );    ## Class static varaibles

    sub import {
        shift;
        my %args = @_;
        ( $endType, $focusTag, $addSoftClipped ) =
          map { $args{$_} } qw( endType focusTag addSoftClipped );
    }

    sub new {
        my $class = shift;
        my %args  = @_;
        my $rec   = $args{rec};
        $rec = defined($rec) ? $rec : "";
        my $self = bless {
            rec          => $rec,
            readID       => "",
            flag         => 0,
            seqname      => "",
            pos          => -1,
            cigar        => "",
            readInPair   => "",
            mateInPair   => "",
            tag          => "",
            hasFocusTag  => 0,
            isMateMapped => 0,
            strand       => '*',
            coord5       => -1,     # coorespond to the 5'end start of the read
            coord3       => -1,
            score        => 0,
        }, $class;
        return $self;
    }

    sub strand {
        my $self = shift;
        return $self->{flag} & 0x10 ? '-' : '+';
    }

    sub isMateMapped {
        my $self         = shift;
        my $isMateMapped = 0;
        if ( $endType eq "PE" ) {
            $isMateMapped = $self->{flag} & 0x08 ? 0 : 1;
        }
        return $isMateMapped;
    }

    sub whichReadInPair {
        my $self = shift;
        my $readInPair;
        if ( $endType eq "SE" ) {
            $readInPair = "R1";
        }
        else {
            $readInPair = ( $self->{flag} & 0x40 ) ? "R1" : "R2";
            if ( !defined($readInPair) ) {
                croak "Paired-end data but cannot tell R1 or R2!";
            }
        }
        return $readInPair;
    }

    sub getMateInPair {
        my $self       = shift;
        my $mateInPair = "";
        if ( $endType eq "PE" ) {
            my $readInPair = $self->{readInPair};
            $mateInPair = { R1 => "R2", R2 => "R1" }->{$readInPair};
        }
        return $mateInPair;
    }

    sub hasFocusTag {
        my $self        = shift;
        my $tag         = $self->{tag};
        my $hasFocusTag = 0;
        my $readInPair  = $self->{readInPair};
        if ( $tag =~ m/$focusTag([1-2])/ ) {
            $hasFocusTag = $readInPair eq "R$1" ? 1 : 0;
        }
        return $hasFocusTag;
    }

    sub coord5 {
        my $self        = shift;
        my $coord5      = $self->{pos};
        my $hasFocusTag = $self->{hasFocusTag};
        my $cigar       = $self->{cigar};
        my $strand      = $self->{strand};
        if ( $strand eq '-' ) {
            if ( $addSoftClipped && $hasFocusTag )
            {    # Note, addSoftClipped only applies when hasFocusTag
                if ( $cigar =~ m/(\d+)S$/ ) {
                    $coord5 += $1;
                }
            }
            while ( $cigar =~ s/(\d+)[MND=X]// ) {
                $coord5 += $1;
            }
        }
        else {
            if ( $addSoftClipped && $hasFocusTag ) {
                if ( $cigar =~ m/^(\d+)S/ ) {
                    $coord5 -= $1;
                }
            }
        }
        if ( $coord5 < 1 ) {
            carp
"Warning: 5'End coord5: $coord5 is negative after adjusting soft clipping!";
            $coord5 = 1;
        }
        return $coord5;
    }

    sub coord3
    { ## for 3End we should not tolarate degenerate sequence mismatches, just ignore soft-clipped part;
        my $self   = shift;
        my $coord3 = $self->{pos};
        my $cigar  = $self->{cigar};
        my $strand = $self->{strand};
        if ( $strand eq '+' ) {
            while ( $cigar =~ s/(\d+)[MND=X]// ) {
                $coord3 += $1;
            }
        }
        if ( $coord3 < 1 ) {
            carp
"Warning: 3'End coord3: $coord3 is negative after adjusting insertion!";
            $coord3 = 1;
        }
        return $coord3;
    }

    sub parseFields {
        my $self   = shift;
        my @fields = split "\t", $self->{rec};
        $self->{readID}       = $fields[0];
        $self->{flag}         = $fields[1];
        $self->{seqname}      = $fields[2];
        $self->{pos}          = $fields[3];
        $self->{score}        = $fields[4];
        $self->{cigar}        = $fields[5];
        $self->{strand}       = $self->strand();
        $self->{isMateMapped} = $self->isMateMapped();
        $self->{readInPair}   = $self->whichReadInPair();
        $self->{mateInPair}   = $self->getMateInPair();
        my @readIDFields = split ":", $self->{readID};
        $self->{tag}         = pop @readIDFields;
        $self->{hasFocusTag} = $self->hasFocusTag();
        $self->{coord5}      = $self->coord5();
        $self->{coord3}      = $self->coord3();
    }
}

{

    package ReadPair;    # For SE data, store R1 only;
    use Carp qw ( croak carp );
    use POSIX qw ( round );

    ## Class static varaibles
    my ( $endType, $focusTag, $outType, $inFh, $outFh, $doesOutNoFocus ) =
      ( "PE", "F", "5End", undef, undef, 0 );

    sub import {
        shift;
        my %args = @_;
        ( $endType, $focusTag, $outType, $inFh, $outFh, $doesOutNoFocus ) =
          map { $args{$_} }
          qw ( endType  focusTag  outType inFh   outFh   doesOutNoFocus );

        if ($debug) {
            print STDERR
"[ReadPair::import] inFh = $inFh, outFh = $outFh, doesOutNoFocus = $doesOutNoFocus\n";
        }
    }

    sub new {
        my $class = shift;
        my %args  = @_;
        my ( $RA, $RB ) = ( $args{RA}, $args{RB} )
          ; ## Because R1 and R2 may be out of order, we denote them as A and B.

        my $self = bless {
            RA        => $RA,
            RB        => $RB,
            R1        => undef,
            R2        => undef,
            hasFocus  => 0,
            focusRead => "",
        }, $class;
        return $self;
    }

    sub orderMates {
        my $self = shift;
        if ( $endType eq "SE" ) { my $self->{R1} = $self->{RA}; }
        else {
            if ( defined( $self->{RA} ) && !defined( $self->{RB} ) )
            {    # a singleton
                if ($debug) {
                    print STDERR
"[ReadPair::orderMates, RA defined but not RB] $self->{RA}->{readInPair}\t$self->{RA}->{seqname}\t$self->{RA}->{strand}\t$self->{RA}->{coord5}\t$self->{RA}->{coord3}\n";
                }
                my $Read       = $self->{RA};
                my $readInPair = $self->{RA}->{readInPair};
                $self->{$readInPair} = $self->{RA};
            }
            elsif ( defined( $self->{RA} ) && defined( $self->{RB} ) )
            {    # paired-end and both aligned
                if ($debug) {
                    print STDERR
"[ReadPair::orderMates, RA RB both defined] $self->{RA}->{readInPair}\t$self->{RA}->{seqname}\t$self->{RA}->{strand}\t$self->{RA}->{coord5}\t$self->{RA}->{coord3}\n";
                    print STDERR
"[ReadPair::orderMates, RA RB both defined] $self->{RB}->{readInPair}\t$self->{RB}->{seqname}\t$self->{RB}->{strand}\t$self->{RB}->{coord5}\t$self->{RB}->{coord3}\n";
                }
                my ( $readInPairA, $readInPairB ) =
                  ( $self->{RA}->{readInPair}, $self->{RB}->{readInPair} );
                $self->{$readInPairA} = $self->{RA};
                $self->{$readInPairB} = $self->{RB};
                if ($debug) {
                    print STDERR
"[ReadPair::orderMates, RA RB both defined] $self->{R1}->{readInPair}\t$self->{R1}->{seqname}\t$self->{R1}->{strand}\t$self->{R1}->{coord5}\t$self->{R1}->{coord3}\n";
                    print STDERR
"[ReadPair::orderMates, RA RB both defined] $self->{R2}->{readInPair}\t$self->{R2}->{seqname}\t$self->{R2}->{strand}\t$self->{R2}->{coord5}\t$self->{R2}->{coord3}\n";
                }
            }
            else {
                croak "Neither of RA and RB is present!";
            }
        }
    }

    sub setFocusRead {
        my $self = shift;
        my @whichHasFocusTag;
        my ( $hasFocus, $focusRead ) = ( 0, "" );
        if ( $endType eq "SE" ) {
            if ( $self->{"R1"}->{hasFocusTag} ) {
                $hasFocus  = 1;
                $focusRead = "R1";
            }
        }
        else {
            @whichHasFocusTag =
              grep { defined( $self->{$_} ) && $self->{$_}->{"hasFocusTag"} }
              qw( R1 R2 )
              ; ## Note, we need to test if $self->{R1/2} is defined not or, otherwise Perl's autovivification will make $self->{R1} and $self->{R2} defined af this line.
            if ( $#whichHasFocusTag > 0 ) {
                croak "Both R1 and R2 have focusTag! This is not possible...";
            }
            elsif ( $#whichHasFocusTag <= -1 ) {
                if ( !defined( $self->{R1} ) ) {
                    carp
"R1 missing and no focusTag found in R2! Labelling this singleton as without focusTag";
                }
                elsif ( !defined( $self->{R2} ) ) {
                    carp
"R2 missing and no focusTag found in R1! Labelling this singleton as without focusTag";
                }
                else {
                    carp
"No focusTag found in either R1 or R2! Labelling this read pair as without focusTag";
                }
            }
            else {
                $hasFocus  = 1;
                $focusRead = $whichHasFocusTag[0];
            }
        }
        $self->{hasFocus}  = $hasFocus;
        $self->{focusRead} = $focusRead;
    }

    sub getReadLocs {
        my $self = shift;
        my ($Read) = @_;
        my ( $readID, $seqname, $strand, $score ) =
          ( $Read->{readID}, $Read->{seqname}, $Read->{strand},
            $Read->{score} );
        my ( $coord5, $coord3 ) = ( $Read->{coord5}, $Read->{coord3} );
        my ( $left, $right ) = ( -1, -1 );
        if ( $outType eq "5End" ) {
            if ( $strand eq '-' ) {
                $right = $coord5;
                $left  = $right - 1;
                $left  = $left < 1 ? 1 : $left;
            }
            else {
                $left  = $coord5;
                $right = $left + 1;
            }
        }
        else {
            $left  = $strand eq '-' ? $coord3 : $coord5;
            $right = $strand eq '-' ? $coord5 : $coord3;
        }
        return {
            seqname => $seqname,
            left    => $left,
            right   => $right,
            readID  => $readID,
            strand  => $strand,
            score   => $score
        };
    }

    sub getFragLocs {
        my $self = shift;
        my ( $R1, $R2 ) = @_;
        my $ReadPair = { R1 => $R1, R2 => $R2 };
        my $ReadPairLocs =
          { R1 => $self->getReadLocs($R1), R2 => $self->getReadLocs($R2) };
        my ( $seqname1, $seqname2 ) =
          ( $ReadPairLocs->{R1}->{seqname}, $ReadPairLocs->{R2}->{seqname} );
        if ( $seqname1 ne $seqname2 ) {
            croak "seqname1: $seqname1 does not match seqname2: $seqname2!";
        }
        my ( $readID1, $readID2 ) =
          ( $ReadPairLocs->{R1}->{readID}, $ReadPairLocs->{R2}->{readID} );
        if ( $readID1 ne $readID2 ) {
            croak "readID1: $readID1 does not match readID2: $readID2!";
        }
        my ( $hasFocus, $focusRead ) =
          ( $self->{hasFocus}, $self->{focusRead} );
        my ( $strand, $left, $right, $score ) = ( "*", -1, -1, 0 );
        if ($hasFocus) {
            $strand = $ReadPairLocs->{$focusRead}->{strand};
            $score  = $ReadPairLocs->{$focusRead}->{score};
            my $noFocusRead = $ReadPair->{$focusRead}->{mateInPair};
            ( $left, $right ) = (
                $ReadPair->{$focusRead}->{coord5},
                $ReadPair->{$noFocusRead}->{coord5}
            );
            if ( $strand eq '-' ) { ( $left, $right ) = ( $right, $left ); }
            if ($debug) {
                print STDERR
"[ReadPair::getFragLocs] hasFocus = $hasFocus, focusRead = $focusRead, noFocusRead = $noFocusRead\n";
                print STDERR
                  "[ReadPair::getFragLocs] left = $left, right = $right\n";
            }
        }
        else {
            my $strand1 = $ReadPairLocs->{R1}->{strand};
            $score = round(
                0.5 * (
                    $ReadPairLocs->{R1}->{score} + $ReadPairLocs->{R2}->{score}
                )
            );
            ( $left, $right ) =
              ( $ReadPair->{R1}->{coord5}, $ReadPair->{R2}->{coord5} );
            if ( $strand1 eq '-' ) { ( $left, $right ) = ( $right, $left ); }

        }
        return {
            seqname => $seqname1,
            left    => $left,
            right   => $right,
            readID  => $readID1,
            strand  => $strand,
            score   => $score
        };
    }

    sub output {
        my $self     = shift;
        my $ReadPair = { R1 => $self->{R1}, R2 => $self->{R2} };
        if ($debug) {
            if ( defined( $ReadPair->{R1} ) ) {
                print STDERR "[ReadPair::output] ReadPair->{R1} is defined.\n";
            }
            else {
                print STDERR
                  "[ReadPair::output] ReadPair->{R1} is not defined.\n";
            }
            if ( defined( $ReadPair->{R2} ) ) {
                print STDERR "[ReadPair::output] ReadPair->{R2} is defined.\n";
            }
            else {
                print STDERR
                  "[ReadPair::output] ReadPair->{R2} is not defined.\n";
            }
        }
        my ( $hasFocus, $focusRead ) =
          ( $self->{hasFocus}, $self->{focusRead} );
        if ( $endType eq "SE" ) {
            if ( $hasFocus && ( $focusRead ne "R1" ) ) {
                croak "SE read hasFocus but focusRead: $focusRead is not R1!";
            }
            my $R           = $ReadPair->{R1};
            my $ReadLocs    = $self->getReadLocs($R);
            my $hasFocusTag = $R->{hasFocusTag};
            my ( $seqname, $left, $right, $readID, $strand, $score ) =
              map { $ReadLocs->{$_} }
              qw( seqname left right readID strand score );
            if ( $hasFocusTag || $doesOutNoFocus ) {
                print $outFh join(
                    "\t",
                    (
                        $seqname, $left - 1, $right - 1, "$readID/1",
                        $score,   $strand
                    )
                  ),
                  "\n";
            }
        }
        else {
            if ( $outType eq "Read" || $outType eq "5End" ) {
                for my $Read ( "R1", "R2" ) {
                    if ( defined( $ReadPair->{$Read} ) ) {
                        my $R           = $ReadPair->{$Read};
                        my $ReadLocs    = $self->getReadLocs($R);
                        my $hasFocusTag = $R->{hasFocusTag};
                        my ( $seqname, $left, $right, $readID, $strand, $score )
                          = map { $ReadLocs->{$_} }
                          qw( seqname left right readID strand score );
                        if ( $hasFocusTag || $doesOutNoFocus ) {
                            my $i = $Read;
                            $i =~ s/^R//;
                            print $outFh join(
                                "\t",
                                (
                                    $seqname,   $left - 1,
                                    $right - 1, "$readID/$i",
                                    $score,     $strand
                                )
                              ),
                              "\n";
                        }
                    }
                }
            }
            elsif ( $outType eq "Frag" ) {
                if ( defined( $ReadPair->{R1} ) && defined( $ReadPair->{R2} ) )
                {
                    my $FragLocs =
                      $self->getFragLocs( $ReadPair->{R1}, $ReadPair->{R2} );
                    my ( $seqname, $left, $right, $readID, $strand, $score ) =
                      map { $FragLocs->{$_} }
                      qw( seqname left right readID strand score );
                    if ( $hasFocus || $doesOutNoFocus ) {
                        print $outFh join(
                            "\t",
                            (
                                $seqname,   $left - 1,
                                $right - 1, "$readID/3",
                                $score,     $strand
                            )
                          ),
                          "\n";
                    }
                }
            }
            elsif ( $outType eq "ReadAndFrag" ) {
                if ( defined( $ReadPair->{R1} ) && defined( $ReadPair->{R2} ) )
                {
                    my $FragLocs =
                      $self->getFragLocs( $ReadPair->{R1}, $ReadPair->{R2} );
                    my ( $seqname, $left, $right, $readID, $strand, $score ) =
                      map { $FragLocs->{$_} }
                      qw( seqname left right readID strand score );
                    if ( $hasFocus || $doesOutNoFocus ) {
                        print $outFh join(
                            "\t",
                            (
                                $seqname,   $left - 1,
                                $right - 1, "$readID/3",
                                $score,     $strand
                            )
                          ),
                          "\n";
                    }
                }
                if ( defined( $ReadPair->{R1} ) && !defined( $ReadPair->{R2} ) )
                {
                    my $R           = $ReadPair->{R1};
                    my $ReadLocs    = $self->getReadLocs($R);
                    my $hasFocusTag = $R->{hasFocusTag};
                    my ( $seqname, $left, $right, $readID, $strand, $score ) =
                      map { $ReadLocs->{$_} }
                      qw( seqname left right readID strand score );
                    if ( $hasFocusTag || $doesOutNoFocus ) {
                        print $outFh join(
                            "\t",
                            (
                                $seqname,   $left - 1,
                                $right - 1, "$readID/1",
                                $score,     $strand
                            )
                          ),
                          "\n";
                    }
                }
                elsif ( !defined( $ReadPair->{R1} )
                    && defined( $ReadPair->{R2} ) )
                {
                    my $R           = $ReadPair->{R2};
                    my $ReadLocs    = $self->getReadLocs($R);
                    my $hasFocusTag = $R->{hasFocusTag};
                    my ( $seqname, $left, $right, $readID, $strand, $score ) =
                      map { $ReadLocs->{$_} }
                      qw( seqname left right readID strand score );
                    if ( $hasFocusTag || $doesOutNoFocus ) {
                        print $outFh join(
                            "\t",
                            (
                                $seqname,   $left - 1,
                                $right - 1, "$readID/2",
                                $score,     $strand
                            )
                          ),
                          "\n";
                    }
                }
            }
        }
    }
}

{

    package Alignments;
    use Carp qw( croak carp );
    use File::Copy qw( move );

    my ( $inFile, $inFh, $nameSortedFile );
    my ( $outFile, $outFh );
    my ($doesOutNoFocus);
    my ( $isNameSorted, $endType, $samFlags, $samSortNCores );
    my ($addSoftClipped);
    my ( $verbose, $debug );

    sub checkArgs {
        my $self = shift;

        if ( !-e $inFile ) {
            croak "$inFile does not exist!";
        }

        if ( defined($outFile) && ( $outFile ne "" ) ) {
            my $outDir = File::Basename::dirname($outFile);
            if ( !-d $outDir ) {
                carp
"Directory for outFile $outFile does not exist! Going to create one...";
                make_path($outDir)
                  or Carp::croak "Cannot create directory $outDir!";
            }
        }
        else { croak "outFile $outFile cannot be empty!"; }

        if ( $endType ne "SE" && $endType ne "PE" ) {
            croak "--endType: $endType should be either PE or SE!";
        }

        if ( $focusTag ne "F" && $focusTag ne "R" ) {
            croak "--focusTag: $focusTag should be either F or R!";
        }
    }

    sub new {
        my $class = shift;
        my %args  = @_;
        (
            $inFile,         $outFile,        $isNameSorted,
            $endType,        $samFlags,       $samSortNCores,
            $addSoftClipped, $doesOutNoFocus, $verbose,
            $debug
          )
          = map { $args{$_} }
          qw (
          inFile outFile isNameSorted endType samFlags samSortNCores
          addSoftClipped doesOutNoFocus verbose debug );

        my $self = bless {}, $class;

        Read->import(
            endType        => $endType,
            focusTag       => $focusTag,
            addSoftClipped => $addSoftClipped
        );
        $self;
    }

    sub sortBam {
        my $self = shift;
        my ( $inFile, $outFile, $by, $samSortNCores ) = @_;
        my $exit;
        my $outTmpPrefix = $outFile;
        $outTmpPrefix =~ s/\.bam$//;
        if ( $by eq "name" ) {
            $exit = system(
"samtools sort -n -o $outFile -T $outTmpPrefix -\@ $samSortNCores $inFile"
            );
        }
        elsif ( $by eq "pos" ) {
            $exit = system(
"samtools sort -o $outFile -T $outTmpPrefix -\@ $samSortNCores $inFile"
            );
        }
        if ( $exit != 0 ) {
            croak "Cannot sort $inFile by $by and save it to $outFile!";
        }
    }

    sub sortBed {
        my $self = shift;
        my ( $inFile, $outFile ) = @_;
        my $exit;
        $exit = system("bedtools sort -i $inFile > $outFile");
        if ( $exit != 0 ) {
            croak "Cannot sort $inFile and save it to $outFile!";
        }
    }

    sub initFh {
        my $self = shift;
        if ( !$isNameSorted ) {
            if ($debug) {
                print STDERR
"Now that $inFile is not sorted by name, we need to sort it first\n";
            }
            $nameSortedFile = $inFile;
            if ( $nameSortedFile !~ s/posSorted\.bam$/nameSorted\.bam/ ) {
                $nameSortedFile =~ s/\.bam$/\.nameSorted\.bam/;
            }
            if ( $inFile eq $nameSortedFile ) {
                Carp::croak "$nameSortedFile is identical to $inFile!";
            }
            if ($debug) {
                print STDERR "[Alignments::initFh] Sorting $inFile\n";
            }
            $self->sortBam( $inFile, $nameSortedFile, "name", $samSortNCores );
            if ($debug) {
                print STDERR
"[Alignments::initFh] Done sorting $inFile by name and save it to $nameSortedFile\n";
            }
            $inFile = $nameSortedFile;
            if ($debug) {
                print STDERR
"[Alignments::initFh] Switching to open $nameSortedFile as the name sorted input\n";
            }
        }

        open( $inFh, "samtools view $samFlags $inFile |" )
          or croak "Cannot open inFile $inFile!";
        open( $outFh, ">", $outFile ) or croak "Cannot open outFile: $outFile!";
    }

    sub initReadPair {
        my $self = shift;
        if ($debug) {
            print STDERR
              "[Alignments::initReadPair] Initializing Class ReadPair...\n";
        }
        ReadPair->import(
            endType        => $endType,
            focusTag       => $focusTag,
            outType        => $outType,
            inFh           => $inFh,
            outFh          => $outFh,
            doesOutNoFocus => $doesOutNoFocus
        );
    }

    sub closeFh {
        my $self = shift;
        if ( defined($inFh) ) {
            if ($debug) {
                if ( !$isNameSorted ) {
                    print STDERR
                      "[Alignments::closeFh] Closing $nameSortedFile...\n";
                }
                else {
                    print STDERR "[Alignments::closeFh] Closing $inFile...\n";
                }
            }
            close $inFh;
        }
        if ( defined($outFh) ) {
            if ($debug) {
                print STDERR "[Alignments::closeFh] Closing $outFile...\n";
            }
            close $outFh;
        }
        if ( !$isNameSorted ) {
            if ($debug) {
                print STDERR
                  "[Alignments::closeFh] Removing $nameSortedFile...\n";
            }
            unlink $nameSortedFile
              or croak
              "Cannot delete the intermediate name sorted $nameSortedFile!";
        }
    }

    sub iterate {
        my $self = shift;
        if ( $endType eq "SE" ) {
            while (<$inFh>) {
                my $rec = $_;
                if ($debug) {
                    print STDERR "[Alignments::iterate::while] rec = $rec";
                }
                my $RA = Read->new( rec => $rec );
                $RA->parseFields();
                my $ReadPair = ReadPair->new( RA => $RA, RB => undef );
                $ReadPair->orderMates();
                $ReadPair->setFocusRead();
                $ReadPair->output();
            }
        }
        else {
            my $currRec = <$inFh>;
            my ( $currRead,   $nextRead );
            my ( $currReadID, $nextReadID );
            my @currReads;
            if ( defined($currRec) ) {
                if ($debug) {
                    print STDERR "[Alignments::iterate] currRec = $currRec";
                }
                $currRead = Read->new( rec => $currRec );
                $currRead->parseFields();
                if ($debug) {
                    print STDERR
"[Alignments::iterate] currRead->readID = $currRead->{readID}\n";
                }
                push @currReads, $currRead;
                $currReadID = $currRead->{readID};
                while (<$inFh>) {
                    my $nextRec = $_;
                    if ($debug) {
                        print STDERR
                          "[Alignments::iterate::while] nextRec = $nextRec";
                    }
                    $nextRead = Read->new( rec => $nextRec );
                    $nextRead->parseFields();
                    if ($debug) {
                        print STDERR
"[Alignments::iterate::while] nextRead->readID = $nextRead->{readID}\n";
                    }
                    $nextReadID = $nextRead->{readID};
                    if ( !defined($currReadID) )
                    { # Case1. following Case3, nextRead is from a new reqd pair, store it in @currReads, and replace currRead with it; No output
                        $currRead   = $nextRead;
                        $currReadID = $currRead->{readID};
                        push @currReads, $currRead;
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case1] currReads has ",
                              $#currReads + 1, " elements\n";
                        }
                    }
                    elsif ( defined($currReadID)
                        && ( $currReadID ne $nextReadID ) )
                    { # Case2. currRead is a singleton, output it anyway; then assign next to current one, and store it in @currReads;
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case2] currReads has ",
                              $#currReads + 1, " elements\n";
                        }
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case2] currReads[0]->readID = $currReads[0]->{readID}\n";
                        }

                        my $ReadPair =
                          ReadPair->new( RA => $currReads[0], RB => undef );
                        $ReadPair->orderMates();
                        $ReadPair->setFocusRead();
                        if ($debug) {
                            if ( defined( $ReadPair->{R1} ) ) {
                                print STDERR
"[Alignments::iterate::while::Case2] ReadPair->R1->readID = $ReadPair->{R1}->{readID}\n";
                            }
                            else {
                                print STDERR
"[Alignments::iterate::while::Case2] ReadPair->{R1} is not defined\n";
                            }
                            if ( defined( $ReadPair->{R2} ) ) {
                                print STDERR
"[Alignments::iterate::while::Case2] ReadPair->R2->readID = $ReadPair->{R2}->{readID}\n";
                            }
                            else {
                                print STDERR
"[Alignments::iterate::while::Case2] ReadPair->{R2} is not defined\n";
                            }
                            print STDERR
"[Alignments::iterate::while::Case2] ReadPair->{hasFocus} = $ReadPair->{hasFocus}\n";
                        }
                        $ReadPair->output();
                        $currRead   = $nextRead;
                        $currReadID = $currRead->{readID};
                        @currReads  = ();
                        push @currReads, $currRead;
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case2] currReads has ",
                              $#currReads + 1, " elements\n";
                        }
                    }
                    elsif ( $currReadID eq $nextReadID )
                    { # Case3. currRead and nextRead come from a pair, so output them and clean up the storage @currReads
                        $currRead = $nextRead;
                        push @currReads, $currRead;
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case3] currReads has ",
                              $#currReads + 1, " elements\n";
                        }
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case3] currReads[0]->readID = $currReads[0]->{readID}\n";
                        }
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case3] currReads[1]->readID = $currReads[1]->{readID}\n";
                        }
                        my $ReadPair = ReadPair->new(
                            RA => $currReads[0],
                            RB => $currReads[1]
                        );
                        $ReadPair->orderMates();
                        $ReadPair->setFocusRead();
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case3] ReadPair->R1->readID = $ReadPair->{R1}->{readID}\n";
                        }
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case3] ReadPair->R2->readID = $ReadPair->{R2}->{readID}\n";
                        }
                        if ($debug) {
                            print STDERR
"[Alignments::iterate::while::Case3] ReadPair->{hasFocus} = $ReadPair->{hasFocus}\n";
                        }
                        $ReadPair->output();
                        ( $currRead,   $nextRead )   = ();
                        ( $currReadID, $nextReadID ) = ();
                        @currReads = ();
                    }
                }
                if (@currReads)
                { ## We need to test if there is anything left in @currReads. If any, we need to process it (it must be Case2: a singleton), otherwise we should just stop.
                    my $ReadPair =
                      ReadPair->new( RA => $currReads[0], RB => undef );
                    $ReadPair->orderMates();
                    $ReadPair->setFocusRead();
                    if ($debug) {
                        print STDERR
"[Alignments::iterate::finishup] ReadPair->R1->readID = $ReadPair->{R1}->{readID}\n";
                    }
                    $ReadPair->output();
                }
            }
        }
    }

    sub postSort {
        my $self          = shift;
        my $posSortedFile = $outFile . ".posSortedTmp.bed";
        $self->sortBed( $outFile, $posSortedFile );
        move( $posSortedFile, $outFile )
          or croak "Cannot rename $posSortedFile to $outFile!";
    }
}

{

    package main;
    my $Alignments = Alignments->new(
        inFile         => $inFile,
        outFile        => $outFile,
        isNameSorted   => $isNameSorted,
        endType        => $endType,
        outType        => $outType,
        addSoftClipped => $addSoftClipped,
        focusTag       => $focusTag,
        samFlags       => $samFlags,
        samSortNCores  => $samSortNCores,
        doesOutNoFocus => $doesOutNoFocus,
        verbose        => $verbose,
        debug          => $debug,
    );

    if ($verbose) { print STDERR "Checking arguments...\n"; }
    $Alignments->checkArgs();
    if ($verbose) { print STDERR "Initializing file handles...\n"; }
    $Alignments->initFh();
    if ($verbose) { print STDERR "Initializing Class ReadPair...\n"; }
    $Alignments->initReadPair();
    if ($verbose) { print STDERR "Preparing seq info...\n"; }
    $Alignments->iterate();
    if ($verbose) { print STDERR "Done with iteration\n" }
    if ($verbose) { print STDERR "Closing file handles...\n" }
    $Alignments->closeFh();
    if ($verbose) { print STDERR "Sorting output bed files by position...\n"; }
    $Alignments->postSort();
    if ($verbose) { print STDERR "All done!\n"; }
}
