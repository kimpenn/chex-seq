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
use POSIX;
use File::Basename;
use File::Path;
use File::Copy;
use Storable;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;

my ( $inFile, $outFile, $outDupFile, $outBedFile, $outStatsFile );
my ( $doesOutDup, $doesOutBed ) = ( 0, 0 );
my $priority            = 'A1,A2,B1,B2,C1,C2';
my $endType             = "PE";
my $ncells              = 1;
my $mindist             = 20;
my $ignoreStrand        = 0;
my $nameSorted          = 0;
my $scheme              = "AS";
my $as_per_mate         = 0;
my $nm_per_mate         = 0;
my $penalty             = 2;
my $flags               = "";
my $focusTag            = "A1:F,A2:F,B1:F,B2:F,C1:R,C2:R";
my $addSoftClipped      = 0;
my $doesOutNoFocus      = 0;
my $doesOutBedZeroCount = 0;
my $ncores              = 1;
our $VERSION = "2.01";
my $debug   = "info";
my $help    = 0;
my $version = 0;
our $LOGGER = get_logger(__PACKAGE__);

sub usage {
    print STDERR <<DOC;
Summary: 
    Remove duplicate priming sites within a specific genomic interval (--mindist) with priority on better barcode/primer quality (--priority). The location of a priming site is determined by the read that carry the focus tag, which is specified by the policy (--focusTag). The threshold of duplicate priming is determined by the number of cells (--ncells): for each cell, more than 2 priming sites is deemed as duplicates.

Usage:
    perl $0 --inFile inFile.bam --outFile NoDup.bam [--outDupFile AdjDup.bam] [--outStatsFile AdjCnts.bed][--ncells 50] [--priority 'A1,A2,B1,B2,C1,C2'] [--endType PE] [--scheme AS] [--penalty 2] [--mindist 20] [--ignoreStrand] [--addSoftClipped] [--focusTag 'A1:F,A2:F,B1:F,B2:F,C1:R,C2:R'] [--flags ''] [--ncores 1] [--doesOutNoFocus] [--doesOutBedZeroCount] [--help] [--version] [--debug debug]

Options:
    --inFile, -i        path of the input BAM;
    --outFile, -o       path of the output BAM that contains reads after adjacent duplicates removed;
    --outDupFile        path of the output BAM that contains filtered reads (default: empty, i.e. no output);
    --outBedFile        path of the output BED that contains nonzero barcode counts per bin per chromosome (default: empty, i.e. no output);
    --outStatsFile      path of the output TSV that contains the final stats (default: empty, i.e. to STDOUT);
    --nameSorted        whether the input is sorted by query name (default: no);
    --ncores            number of threads for `samtools sort` (default: 1);
    --scheme            score scheme, choose from 'AS', 'heuristic_raw' and 'heuristic_nonoverlap' (default: 'AS');
    --AS-per-mate       boolean switch, if on, a read pair's AS equals to the sum of two mates' AS. This should be 
                        turned off if STAR, otherwise if Bowtie2 or BWA inputs (default: off); 
    --NM-per-mate       boolean switch, if on, instead of 'nM', a read pair's #. of mismatches equals to the sum of 
                        two mates' 'NM' (default: off)--this requires presence of tag 'NM'; note that STAR by default does NOT output 'NM' (default: off); 
    --penalty           linear penalty of # of mismatches in heuristic scores (default: 2);
    --priority          order of the barcode/primer quality (default: 'A1,A2,B1,B2,C1,C2' which means A1>A2>B1>B2>C1>C2);
    --ncells            number of cells (default: 1 if unset, otherwise 50 if 'NA');
    --endType           'SE' or 'PE' (default: 'PE');
    --mindist           minimum distance (bp) to be considered duplicate priming (default: 20); 
    --ignoreStrand      whether or not ignore the strands when counting duplicate barcodes? (default: false); 
    --addSoftClipped    when a barcode's 5' end softclipped by the aligner, whether or not recover the softclipped length of the barcode? (default: false);
    --focusTag       a string of key-value pairs denoting which read (for PE) has the barcode for focus location tracking (default: 'A1:F,A2:F,B1:F,B2:F,C1:R,C2:R');
    --flags             additional flags for `samtools view` (default: '');
    --doesOutNoFocus    whether output read (pair) that contains no focusTag (default: false);
    --doesOutBedZeroCount   whether to output bins with no focusTag count (default: false);
    --debug             choose from 'fatal', 'error', 'warn', 'info', 'debug' and 'trace' (default: 'info');
    --help, -h          print usage and exit;
    --version, -v       print version and exit;
DOC
}

GetOptions(
    "inFile|i=s"          => \$inFile,
    "outFile|o=s"         => \$outFile,
    "outDupFile=s"        => \$outDupFile,
    "outBedFile=s"        => \$outBedFile,
    "outStatsFile=s"      => \$outStatsFile,
    "nameSorted"          => \$nameSorted,
    "priority=s"          => \$priority,
    "endType=s"           => \$endType,
    "mindist=i"           => \$mindist,
    "ignoreStrand"        => \$ignoreStrand,
    "addSoftClipped"      => \$addSoftClipped,
    "focusTag=s"          => \$focusTag,
    "flags=s"             => \$flags,
    "ncells=s"            => \$ncells,
    "scheme=s"            => \$scheme,
    "penalty=s"           => \$penalty,
    "ncores=i"            => \$ncores,
    "doesOutNoFocus"      => \$doesOutNoFocus,
    "doesOutBedZeroCount" => \$doesOutBedZeroCount,
    "help|h"              => \$help,
    "version|v"           => \$version,
    "debug=s"             => \$debug,
) or ( &usage() && exit(-1) );

( &usage() && exit(0) ) if $help;
( ( print "$0 v$VERSION\n" ) && exit(0) ) if $version;

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

$LOGGER->info( "{ inFile = $inFile, outFile = $outFile, outDupFile = "
      . ( defined($outDupFile) ? $outDupFile : "" )
      . ", outBedFile = "
      . ( defined($outBedFile) ? $outBedFile : "" )
      . ", outStatsFile = "
      . ( defined($outStatsFile) ? $outStatsFile : "" )
      . ", nameSorted = $nameSorted, ncells = $ncells, priority = $priority, endType = $endType, mindist = $mindist, ignoreStrand = $ignoreStrand, addSoftClipped = $addSoftClipped, focusTag = $focusTag, flags = $flags, scheme = $scheme, penalty = $penalty, ncores = $ncores, doesOutNoFocus = $doesOutNoFocus, VERSION = $VERSION, debug = $debug }\n"
);

{

    package Read;
    my ( $endType, $addSoftClipped, $nm_per_mate, $as_per_mate, $qranks,
        $focusTagMap );

    sub init_class {
        my $class = shift;
        $LOGGER->info("Initializing Read...\n");
        my %args = @_;
        (
            $endType, $qranks, $focusTagMap, $addSoftClipped, $nm_per_mate,
            $as_per_mate
          )
          = map { $args{$_} }
          qw( endType qranks focusTagMap addSoftClipped nm_per_mate as_per_mate );
    }

    sub new {
        my $class = shift;
        my ($rec) = @_;
        my $self  = bless {
            rec          => $rec,
            readID       => undef,
            flag         => undef,
            chr          => undef,
            pos          => undef,
            cigar        => undef,
            tags         => undef,
            clip_l       => 0,
            clip_r       => 0,
            mappedLength => undef,
            whichInPair  => "R1",
            hasFocus     => 0,
            isRev        => undef,
            coord5       => undef,
            strandedness => undef,
            qual         => undef,
            NM           => undef,
            nM           => undef,
            AS           => undef,
        }, $class;
        $self;
    }

    sub parse_flag {
        my $self = shift;
        ( $self->{whichInPair} = $self->{flag} & 0x40 ? "R1" : "R2" )
          if $endType eq "PE";
        $self->{isRev} = $self->{flag} & 0x10 ? 1 : 0;
    }

    sub get_strand {
        my $self = shift;
        return $self->{isRev} ? '-' : '+';
    }

    sub parse_cigar {
        my $self   = shift;
        my $cigar  = $self->{cigar};
        my $maplen = 0;
        $self->{clip_r} = $1 if $cigar =~ m/(\d+)S$/;
        $self->{clip_l} = $1 if $cigar =~ m/^(\d+)S/;
        while ( $cigar =~ s/(\d+)[MX=]// ) { $maplen += $1 }
        $self->{mappedLength} = $maplen;
    }

    sub set_coord5 {
        my $self   = shift;
        my $coord5 = $self->{pos};
        if ( $self->{isRev} ) {
            $coord5 += $self->{clip_r} if $addSoftClipped;
            $coord5 += $self->{mappedLength};
        }
        else {
            $coord5 -= $self->{clip_l} if $addSoftClipped;
            if ( $coord5 < 0 ) {
                $LOGGER->warn(
"5' end coordinate: $coord5 is negative after adjusting softclipping!\n"
                );
                $coord5 = 0;
            }
        }
        $self->{coord5} = $coord5;
    }

    sub set_hasFocus {
        my $self = shift;
        my $t    = ${$focusTagMap}{ $self->{qual} };
        if ( defined($t) ) {
            if ( $self->{strandedness} =~ m/$t([12])/ ) {
                $self->{hasFocus} = $self->{whichInPair} eq "R$1" ? 1 : 0;
            }
        }
    }

    sub parse_NM {
        my $self = shift;
        my $NM;
        my @t = grep { /^NM/ } @{ $self->{tags} };
        $NM = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{NM} = $NM;
    }

    sub parse_nM {
        my $self = shift;
        my $nM;
        my @t = grep { /^nM/ } @{ $self->{tags} };
        $nM = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{nM} = $nM;
    }

    sub parse_AS {
        my $self = shift;
        my $AS;
        my @t = grep { /^AS/ } @{ $self->{tags} };
        $AS = substr( $t[0], 5 ) if ( $#t == 0 );
        $self->{AS} = $AS;
    }

    sub parse_fields {
        my $self     = shift;
        my @F        = split( "\t", $self->{rec} );
        my $readName = $F[0];
        $self->{flag}  = $F[1];
        $self->{chr}   = $F[2];
        $self->{pos}   = $F[3];
        $self->{cigar} = $F[5];
        @F             = @F[ 11 .. $#F ];
        $self->{tags}  = \@F;
        $self->parse_flag();     # set $isRev, $whichInPair
        $self->parse_cigar();    # set $mappedLength, $clip_l, $clip_r
        $self->set_coord5();     # set the position of the 5' end
        my @t = split( ":", $readName, -1 );
        $self->{strandedness} = pop @t;    # get F1R2 or F2R1
        $self->{qual}         = pop @t;    # get A1, B2, ...
        my $readID = shift @t;
        $self->{readID} = ( split( "-", $readName ) )[0];
        $self->set_hasFocus();
        $self->parse_AS();
        $nm_per_mate ? $self->parse_NM() : $self->parse_nM();
    }
}

{

    package ReadPair;                      # For SE data, store R1 only;

    ## Class static varaibles
    my (
        $endType,        $focusTagMap, $addSoftClipped, $inFh,
        $outFh,          $doesOutDup,  $outDupFh,       $outStatsFh,
        $doesOutNoFocus, $scheme,      $penalty,        $nm_per_mate,
        $as_per_mate,    %stats,       $qranks,
    );

    sub init_class {
        shift;
        my %args = @_;
        (
            $endType,        $focusTagMap, $addSoftClipped, $inFh,
            $outFh,          $doesOutDup,  $outDupFh,       $outStatsFh,
            $doesOutNoFocus, $scheme,      $penalty,        $nm_per_mate,
            $as_per_mate,    $qranks
          )
          = map { $args{$_} }
          qw ( endType   focusTagMap   addSoftClipped inFh
          outFh   doesOutDup outDupFh outStatsFh
          doesOutNoFocus scheme penalty nm_per_mate
          as_per_mate qranks );

        Read->init_class(
            endType        => $endType,
            qranks         => $qranks,
            focusTagMap    => $focusTagMap,
            addSoftClipped => $addSoftClipped,
            nm_per_mate    => $nm_per_mate,
            as_per_mate    => $as_per_mate,
        );

        $LOGGER->debug(
"inFh = $inFh, outFh = $outFh, outDupFh = $outDupFh, outStatsFh = $outStatsFh, doesOutDup = $doesOutDup, doesOutNoFocus = $doesOutNoFocus, scheme = $scheme, penalty = $penalty, nm_per_mate = $nm_per_mate\n"
        );
    }

    sub new {
        my $class = shift;
        my %args  = @_;
        my ( $R1, $R2 ) = map { $args{$_} } qw( R1 R2 );
        my $self = bless {
            R1               => $R1,
            R2               => $R2,
            rawLength        => -1,
            nonoverlapLength => -1,
            mismatches       => -1,
            score            => -1,
            hasFocus         => 0,
            focusRead        => undef,
            qual             => '',
            qrank            => 0,
            chr              => undef,
            strand           => undef,
            coord5           => undef,
            isDup            => 0,
        }, $class;
        $self;
    }

    sub is_discordant {
        my $self = shift;
        return 1 if ( $self->{R1}->{chr} ne $self->{R2}->{chr} );
        return 0;
    }

    sub parse_mappedLength {
        my $self = shift;
        if ( $endType eq "PE" ) {
            if ( defined( $self->{R1} ) && defined( $self->{R2} ) ) {
                $self->{rawLength} =
                  $self->{R1}->{mappedLength} + $self->{R2}->{mappedLength};
                if ( $self->is_discordant() ) {
                    $self->{nonoverlapLength} = -1;
                    return;
                }
                my $overlap =
                    $self->{R1}->{isRev}
                  ? $self->{R2}->{pos} + $self->{R2}->{mappedLength} -
                  $self->{R1}->{pos}
                  : $self->{R1}->{pos} +
                  $self->{R1}->{mappedLength} -
                  $self->{R2}->{pos};
                $overlap = $overlap > 0 ? $overlap : 0;
                $self->{nonoverlapLength} = $self->{rawLength} - $overlap;
                return;
            }
            elsif ( defined( $self->{R1} ) ) {
                $self->{rawLength} = $self->{nonoverlapLength} =
                  $self->{R1}->{mappedLength};
                return;
            }
            else {
                $self->{rawLength} = $self->{nonoverlapLength} =
                  $self->{R2}->{mappedLength};
                return;
            }
        }
        else {
            $self->{rawLength} = $self->{nonoverlapLength} =
              $self->{R1}->{mappedLength};
            return;
        }
    }

    sub parse_mismatches {
        my $self = shift;
        if ( $endType eq "PE" ) {
            if ($nm_per_mate) {
                my $s1 = defined( $self->{R1} ) ? $self->{R1}->{NM} : 0;
                my $s2 = defined( $self->{R2} ) ? $self->{R2}->{NM} : 0;
                $self->{mismatches} = $s1 + $s2;
            }
            else {
                $self->{mismatches} =
                  defined( $self->{R1} )
                  ? $self->{R1}->{nM}
                  : $self->{R2}->{nM};
            }
        }
        else {
            if ($nm_per_mate) {
                $self->{mismatches} = $self->{R1}->{NM};
            }
            else {
                $self->{mismatches} = $self->{R1}->{nM};
            }
        }
    }

    sub parse_AS {
        my $self = shift;
        my @t;
        if ( $endType eq "PE" ) {
            if ($as_per_mate) {
                my $s1 = defined( $self->{R1} ) ? $self->{R1}->{AS} : 0;
                my $s2 = defined( $self->{R2} ) ? $self->{R2}->{AS} : 0;
                $self->{score} = $s1 + $s2;
            }
            else {
                $self->{score} =
                  defined( $self->{R1} )
                  ? $self->{R1}->{AS}
                  : $self->{R2}->{AS};
            }
        }
        else {
            $self->{score} = $self->{R1}->{AS};
        }
    }

    sub set_primerqual {
        my $self = shift;
        if ( $endType eq "PE" ) {
            if ( defined( $self->{R1} ) ) {
                $self->{qual} = $self->{R1}->{qual}
                  if defined( $self->{R1}->{qual} );
            }
            else {
                $self->{qual} = $self->{R2}->{qual}
                  if defined( $self->{R2}->{qual} );
            }
        }
        else {
            $self->{qual} = $self->{R1}->{qual}
              if defined( $self->{R1}->{qual} );
        }
        $self->{qrank} = ${$qranks}{ $self->{qual} };
    }

    sub set_mapqual {
        my $self = shift;
        $self->parse_mappedLength();
        $self->parse_mismatches();
        if ( $scheme eq "AS" ) {
            $self->parse_AS();
        }
        else {
            if ( $scheme eq "heuristic_raw" ) {
                $self->{score} =
                  $self->{rawLength} - $penalty * $self->{mismatches};
            }
            else {
                $self->{score} =
                  $self->{nonoverlapLength} - $penalty * $self->{mismatches};
            }
        }
    }

    sub set_focusLoc {
        my $self = shift;
        if ( $endType eq "SE" ) {
            my $read = $self->{R1};
            if ( $read->{hasFocus} ) {
                $self->{hasFocus}  = 1;
                $self->{focusRead} = "R1";
                $self->{chr}       = $read->{chr};
                $self->{strand}    = $read->get_strand();
                $self->{coord5}    = $read->{coord5};
                $LOGGER->trace(
"focusRead = $self->{focusRead}, chr = $read->{chr}, strand = $read->{strand}, coord5 = $read->{coord5}\n"
                );
            }
        }
        else {
            my @whichHasFocus =
              grep { defined( $self->{$_} ) && $self->{$_}->{hasFocus} }
              qw( R1 R2 );
            $LOGGER->fatal("Both R1 and R2 have focus! This is not possible...")
              && die()
              if $#whichHasFocus == 1;
            $LOGGER->warn(
"No focus found in either R1 or R2 or a singleton with no focus!\n"
              )
              && return
              if $#whichHasFocus == -1;

            $self->{hasFocus}  = 1;
            $self->{focusRead} = $whichHasFocus[0];
            my $read = $self->{ $self->{focusRead} };
            $self->{chr}    = $read->{chr};
            $self->{strand} = $read->get_strand();
            $self->{coord5} = $read->{coord5};
            $LOGGER->trace(
"focusRead = $self->{focusRead}, chr = $self->{chr}, strand = $self->{strand}, coord5 = $self->{coord5}\n"
            );
        }
    }

    sub output_reads {
        my $self = shift;
        if ( $self->{hasFocus} || $doesOutNoFocus ) {
            if ( $endType eq "SE" ) {
                if ( !$self->{isDup} ) {
                    print $outFh $self->{R1}->{rec}, "\n";
                }
                elsif ($doesOutDup) {
                    print $outDupFh $self->{R1}->{rec}, "\n";
                }
            }
            else {
                if ( !$self->{isDup} ) {
                    if ( defined( $self->{R1} ) ) {
                        print $outFh $self->{R1}->{rec}, "\n";
                    }
                    if ( defined( $self->{R2} ) ) {
                        print $outFh $self->{R2}->{rec}, "\n";
                    }
                }
                elsif ($doesOutDup) {
                    if ( defined( $self->{R1} ) ) {
                        print $outDupFh $self->{R1}->{rec}, "\n";
                    }
                    if ( defined( $self->{R2} ) ) {
                        print $outDupFh $self->{R2}->{rec}, "\n";
                    }
                }
            }
        }
    }

    sub do_stats {
        my $self = shift;
        $stats{nReadPairs}++;
        if ( $self->{hasFocus} ) {
            $stats{nHasFocus}{both}++;
            $stats{nHasFocus}{ $self->{focusRead} }++;
            if ( $self->{isDup} ) {
                $stats{nHasFocus}{isDup}++;
            }
            else {
                $stats{nHasFocus}{isNotDup}++;
            }
        }
        else {
            $stats{nHasNoFocus}{both}++;
        }
    }

    sub output_stats {
        my $self = shift;
        $LOGGER->info("Writing final stats report...\n");
        my $nReadPairs = defined( $stats{nReadPairs} ) ? $stats{nReadPairs} : 0;
        my $nHasFocus =
          defined( $stats{nHasFocus}{both} ) ? $stats{nHasFocus}{both} : 0;
        my $nHasFocusInR1 =
          defined( $stats{nHasFocus}{R1} ) ? $stats{nHasFocus}{R1} : 0;
        my $nHasFocusInR2 =
          defined( $stats{nHasFocus}{R2} ) ? $stats{nHasFocus}{R2} : 0;
        my $nHasFocusIsDup =
          defined( $stats{nHasFocus}{isDup} ) ? $stats{nHasFocus}{isDup} : 0;
        my $nHasFocusIsNotDup =
          defined( $stats{nHasFocus}{isNotDup} )
          ? $stats{nHasFocus}{isNotDup}
          : 0;
        my $nHasNoFocus =
          defined( $stats{nHasNoFocus}{both} ) ? $stats{nHasNoFocus}{both} : 0;
        print $outStatsFh
"nReadPairs\tnHasFocus\tnHasFocusInR1\tnHasFocusInR2\tnHasFocusIsDup\tnHasFocusIsNoDup\tnHasNoFocus\n";
        print $outStatsFh
"$nReadPairs\t$nHasFocus\t$nHasFocusInR1\t$nHasFocusInR2\t$nHasFocusIsDup\t$nHasFocusIsNotDup\t$nHasNoFocus\n";
    }
}

{

    package SamReader;
    my ( $inFile, $inFh, $nameSortedFile );
    my ( $outFile,      $outFh );
    my ( $outDupFile,   $outDupFh, $doesOutDup, $doesOutNoFocus );
    my ( $outBedFile,   $outBedFh, $doesOutBed, $doesOutBedZeroCount );
    my ( $outStatsFile, $outStatsFh );
    my ( @SeqLevels,  @SeqLengths, %cnts );
    my ( $nameSorted, $endType,    $flags, $ncores );
    my ( $ncells,     $mindist,    $ignoreStrand );
    my ( $priority, %qranks );
    my ( $focusTag, %focusTagMap );

    sub new {
        my $class = shift;
        my %args  = @_;
        (
            $inFile,     $outFile,        $outDupFile,
            $outBedFile, $outStatsFile,   $doesOutBedZeroCount,
            $nameSorted, $priority,       $endType,
            $flags,      $ncores,         $ncells,
            $mindist,    $ignoreStrand,   $addSoftClipped,
            $focusTag,   $doesOutNoFocus, $scheme,
            $penalty,    $nm_per_mate,    $as_per_mate
          )
          = map { $args{$_} }
          qw (
          inFile outFile outDupFile
          outBedFile outStatsFile doesOutBedZeroCount
          nameSorted priority endType
          flags ncores ncells
          mindist ignoreStrand addSoftClipped
          focusTag doesOutNoFocus scheme
          penalty nm_per_mate as_per_mate
        );
        my $self = bless { readPairs => [], }, $class;
        $self;
    }

    sub check_args {
        my $self = shift;
        $LOGGER->info("Checking arguments...\n");
        if ( !-e $inFile ) {
            $LOGGER->fatal("$inFile does not exist!\n") && die();
        }

        if ( defined($outFile) && ( $outFile ne "" ) ) {
            my $outDir = File::Basename::dirname($outFile);
            if ( !-d $outDir ) {
                $LOGGER->fatal(
"Directory for outFile $outFile does not exist! Going to create one...\n"
                ) && die();
                make_path($outDir)
                  or $LOGGER->fatal("Cannot create directory $outDir!\n")
                  && die();
            }
        }
        else { $LOGGER->fatal("outFile $outFile cannot be empty!\n") && die(); }

        if ( $endType ne "SE" && $endType ne "PE" ) {
            $LOGGER->fatal("--endType: $endType should be either PE or SE!\n")
              && die();
        }

        if ( int($mindist) != $mindist || $mindist < 1 ) {
            $LOGGER->fatal("mindist should be integer and at least 1 (bp)\n")
              && die();
        }

        if ( $ncells eq "NA" ) {
            $LOGGER->warn("Cell count is NA! Setting it to 50...\n");
            $ncells = 50;
        }
        else {
            if ( $ncells =~ /[A-Za-z]/ ) {
                $LOGGER->fatal("Cell count is not a number! Exiting...\n")
                  && die();
            }
            $ncells = int($ncells);
            if ( $ncells < 0 ) {
                $LOGGER->fatal("Cell count must be nonnegative! Exiting...\n")
                  && die();
            }
            if ( $ncells == 0 ) {
                $LOGGER->warn(
"Cell count is 0! We don't do duplicate removal in this case...\n"
                );
                $ncells = -1;
            }
        }
    }

    sub parse_qualfocus {
        my @quals = split( ",", $priority );
        $LOGGER->debug("quals = @quals\n");
        $qranks{''} = 0;
        for my $i ( 0 .. $#quals ) {
            $qranks{ $quals[$i] } = $i + 1;
        }
        for my $i (@quals) {
            $LOGGER->debug("qranks: $i => $qranks{$i}\n");
        }
        my @qualfocus = split( ",", $focusTag );
        for my $qf (@qualfocus) {
            my @t = split( ":", $qf );
            if ( $#t != 1 ) {
                $LOGGER->fatal("--focusTag: $qf is not key:value pair!\n")
                  && die();
            }
            if ( $t[1] ne "F" && $t[1] ne "R" ) {
                $LOGGER->fatal("--focusTag: $t[1] should be either F or R!\n")
                  && die();
            }
            $focusTagMap{ $t[0] } = $t[1];
        }
        for my $i (@quals) {
            $LOGGER->debug("focusTagMap: $i => $focusTagMap{$i}\n");
        }
    }

    sub init_fh {
        my $self = shift;
        $LOGGER->info("Initializing file handles...\n");
        if ( !$nameSorted ) {
            $LOGGER->debug(
"Now that $inFile is not sorted by name, we need to sort it first\n"
            );
            $nameSortedFile = $inFile;
            if ( $nameSortedFile !~ s/posSorted\.bam$/nameSorted\.bam/i ) {
                $nameSortedFile =~ s/\.bam$/\.nameSorted\.bam/i;
            }
            if ( $inFile eq $nameSortedFile ) {
                $LOGGER->fatal("$nameSortedFile is identical to $inFile!")
                  && die();
            }
            $LOGGER->debug("Sorting $inFile\n");
            $self->sort_bam( $inFile, $nameSortedFile, "name", $ncores );
            $LOGGER->debug(
                "Done sorting $inFile by name and save it to $nameSortedFile\n"
            );
            $inFile = $nameSortedFile;
            $LOGGER->debug(
                "Switching to open $nameSortedFile as the name sorted input\n");
        }

        open $inFh, "samtools view $flags $inFile |"
          or $LOGGER->fatal("Cannot open inFile $inFile!") && die();
        open $outFh, "|-", "samtools view -Sb - > $outFile"
          or $LOGGER->fatal("Cannot open outFile: $outFile!") && die();

        if ( defined($outDupFile) && ( $outDupFile ne "" ) ) {
            my $outDupDir = File::Basename::dirname($outDupFile);
            if ( !-d $outDupDir ) {
                $LOGGER->fatal(
"Directory for outDupFile $outDupFile does not exist! Going to create one...\n"
                ) && die();
                make_path($outDupDir)
                  or $LOGGER->fatal("Cannot create directory $outDupDir!\n")
                  && die();
            }
            $doesOutDup = 1;
            $LOGGER->debug(
                "outDupFile is defined, so we set doesOutDup = 1.\n");
            open $outDupFh, "|-", "samtools view -Sb - > $outDupFile"
              or $LOGGER->fatal("Cannot open outDupFile: $outDupFile!")
              && die();
        }

        if ( defined($outBedFile) && ( $outBedFile ne "" ) ) {
            my $outBedDir = File::Basename::dirname($outBedFile);
            if ( !-d $outBedDir ) {
                $LOGGER->fatal(
"Directory for outBedFile $outBedFile does not exist! Going to create one...\n"
                ) && die();
                make_path($outBedDir)
                  or $LOGGER->fatal("Cannot create directory $outBedDir!\n")
                  && die();
            }
            $doesOutBed = 1;
            $LOGGER->debug(
                "outBedFile is defined, so we set doesOutBed = 1.\n");
            open $outBedFh, ">", $outBedFile
              or $LOGGER->fatal("Cannot open outBedFile: $outBedFile!")
              && die();
        }

        if ( defined($outStatsFile) && ( $outStatsFile ne "" ) ) {
            my $outStatsDir = File::Basename::dirname($outStatsFile);
            if ( !-d $outStatsDir ) {
                $LOGGER->fatal(
"Directory for outStatsFile $outStatsFile does not exist! Going to create one...\n"
                ) && die();
                make_path($outStatsDir)
                  or $LOGGER->fatal("Cannot create directory $outStatsDir!\n")
                  && die();
            }
            open $outStatsFh, ">", $outStatsFile
              or $LOGGER->fatal("Cannot open outStatsFile: $outStatsFile!")
              && die();
        }
        else {
            $outStatsFh = *STDOUT;
        }
    }

    sub init_class {
        my $self = shift;
        $LOGGER->debug("Initializing ReadPair...\n");

        ReadPair->init_class(
            endType        => $endType,
            focusTag       => $focusTag,
            addSoftClipped => $addSoftClipped,
            qranks         => \%qranks,
            focusTagMap    => \%focusTagMap,
            inFh           => $inFh,
            outFh          => $outFh,
            outDupFh       => $outDupFh,
            outStatsFh     => $outStatsFh,
            doesOutDup     => $doesOutDup,
            doesOutNoFocus => $doesOutNoFocus,
            scheme         => $scheme,
            nm_per_mate    => $nm_per_mate,
            as_per_mate    => $as_per_mate,
            penalty        => $penalty,
        );
    }

    sub init_seqinfo {
        my $self = shift;
        $LOGGER->info("Preparing seq info...\n");
        $LOGGER->debug("Preparing seq info...\n");
        open my $inFh0, "samtools view -H $inFile |"
          or $LOGGER->fatal("Cannot open inFile $inFile!") && die();
        while (<$inFh0>) {
            print $outFh $_;
            print $outDupFh $_ if $doesOutDup;
            if (/^\@SQ\s+SN:(chr.+)\s+LN:(\d+)/) {
                push @SeqLevels,  $1;
                push @SeqLengths, $2;
            }
        }
        close $inFh0;
    }

    sub fin_fh {
        my $self = shift;
        $LOGGER->info("Closing file handles...\n");
        if ( defined($inFh) ) {
            if ( !$nameSorted ) {
                $LOGGER->debug("Closing $nameSortedFile...\n");
            }
            else {
                $LOGGER->debug("Closing $inFile...\n");
            }
            close $inFh;
        }
        if ( defined($outFh) ) {
            $LOGGER->debug("Closing $outFile...\n");
            close $outFh;
        }
        if ( $doesOutDup && defined($outDupFh) ) {
            $LOGGER->debug("Closing $outDupFile...\n");
            close $outDupFh;
        }
        if ( $doesOutBed && defined($outBedFh) ) {
            $LOGGER->debug("Closing $outBedFile...\n");
            close $outBedFh;
        }
        if ( !$nameSorted ) {
            $LOGGER->debug("Removing $nameSortedFile...\n");
            unlink $nameSortedFile
              or $LOGGER->fatal(
                "Cannot delete the intermediate name sorted $nameSortedFile!")
              && die
              ();
        }
    }

    sub init_cnts {
        my $self = shift;
        $LOGGER->info("Initializing counts...\n");
        my %x = map { $_ => {} } @SeqLevels;
        if ($ignoreStrand) {
            $cnts{"*"} = Storable::dclone( \%x );
        }
        else {
            ## Make a deep copy here; otherwise '-' will overwrite '+' strand!
            $cnts{"+"} = Storable::dclone( \%x );
            $cnts{"-"} = Storable::dclone( \%x );
        }
        $LOGGER->info("Done initializing counts.\n");
    }

    sub test_dup {
        my $self = shift;
        my ( $readPair, $cutoff ) = @_;
        if ( $readPair->{hasFocus} ) {
            my ( $chr, $strand, $coord5 ) =
              ( $readPair->{chr}, $readPair->{strand}, $readPair->{coord5} );
            my $binIdx = POSIX::round( $coord5 / $mindist );
            my $c =
                $ignoreStrand
              ? $cnts{"*"}{$chr}{$binIdx}++
              : $cnts{$strand}{$chr}{$binIdx}++;
            $c = defined($c) ? $c : 0;
            $readPair->{isDup} = 1 if $cutoff != -1 && $c >= $cutoff;
            $LOGGER->debug(
"chr = $chr, strand = $strand, coord5 = $coord5, binIdx = $binIdx, c = $c, qual = $readPair->{qual}, qrank = $readPair->{qrank}, score = $readPair->{score}, rawLength = $readPair->{rawLength}, nonoverlapLength = $readPair->{nonoverlapLength}, mismatches = $readPair->{mismatches}, hasFocus = $readPair->{hasFocus}, isDup = $readPair->{isDup}\n"
            );
        }
    }

    sub slurp_reads {
        my $self = shift;
        $LOGGER->info("Slurping all reads...\n");
        if ( $endType eq "SE" ) {
            while (<$inFh>) {
                chomp;
                my $rec = $_;
                $LOGGER->trace("rec = $rec");
                my $read = Read->new($rec);
                $read->parse_fields();
                my $readPair = ReadPair->new( R1 => $read );
                $readPair->set_focusLoc();
                $readPair->set_primerqual();
                $readPair->set_mapqual();
                push @{ $self->{readPairs} }, $readPair;
            }
        }
        else {
            my $currRec = <$inFh>;
            chomp($currRec);
            my ( $currRead, $nextRead );
            my @currReads;
            if ( defined($currRec) ) {
                $LOGGER->trace("currRec = $currRec");
                $currRead = Read->new($currRec);
                $currRead->parse_fields();
                $LOGGER->trace("currRead->readID = $currRead->{readID}\n");
                while (<$inFh>) {
                    chomp;
                    my $nextRec = $_;
                    $LOGGER->trace("nextRec = $nextRec");
                    $nextRead = Read->new($nextRec);
                    $nextRead->parse_fields();
                    $LOGGER->trace("nextRead->readID = $nextRead->{readID}\n");

                    if ( !defined($currRead) )
                    { # Case1, following Case3, nextRead is from a new read pair, store it in @currReads, and replace currRead with it; No output
                        $currRead = $nextRead;
                        next;
                    }
                    if ( ( $currRead->{readID} ne $nextRead->{readID} ) )
                    { # Case2, currRead is a singleton, output it anyway; then assign next to current one, and store it in @currReads;
                        my $readPair = ReadPair->new(
                            $currRead->{whichInPair} => $currRead );
                        $readPair->set_focusLoc();
                        $readPair->set_primerqual();
                        $readPair->set_mapqual();
                        push @{ $self->{readPairs} }, $readPair;
                        $currRead = $nextRead;
                        next;
                    }
                    if ( $currRead->{readID} eq $nextRead->{readID} )
                    { # Case3, currRead and nextRead come from a pair, so output them and clean up the storage @currReads
                        my $readPair = ReadPair->new(
                            $currRead->{whichInPair} => $currRead,
                            $nextRead->{whichInPair} => $nextRead
                        );
                        $readPair->set_focusLoc();
                        $readPair->set_primerqual();
                        $readPair->set_mapqual();
                        push @{ $self->{readPairs} }, $readPair;
                        ( $currRead, $nextRead ) = ( undef, undef );
                        next;
                    }
                }
                if ( defined($currRead) )
                { ## Case2, a singleton leftover, otherwise we should just stop.
                    my $readPair =
                      ReadPair->new( $currRead->{whichInPair} => $currRead );
                    $readPair->set_focusLoc();
                    $readPair->set_primerqual();
                    $readPair->set_mapqual();
                    push @{ $self->{readPairs} }, $readPair;
                }
            }
        }
        $LOGGER->info("Done slurping all reads.\n");
    }

    sub sort_quals {
        my $self = shift;
        $LOGGER->info(
"Sorting reads (read pairs) by barcode/primer quality and mapping quality...\n"
        );
        @{ $self->{readPairs} } = sort {
                 $a->{qrank} <=> $b->{qrank}
              or $b->{score} <=> $a->{score}
              or $b->{nonoverlapLength} <=> $a->{nonoverlapLength}
        } @{ $self->{readPairs} };
        $LOGGER->info(
            "Done sorting barcode/primer quality and mapping quality.\n");
    }

    sub detect_dup {
        my $self = shift;
        $LOGGER->info("Detecting duplicates...\n");
        my $cutoff = $ncells == -1 ? -1 : 2 * $ncells;
        $LOGGER->debug("ncells = $ncells, cutoff = $cutoff\n");
        for my $readPair ( @{ $self->{readPairs} } ) {
            $self->test_dup( $readPair, $cutoff );
            $readPair->output_reads();
            $readPair->do_stats();
        }
        $LOGGER->info("Done with duplicates detection.\n");
    }

    sub sort_bam {
        my $self = shift;
        my ( $inFile, $outFile, $by, $ncores ) = @_;
        my $exit;
        my $outTmpPrefix = $outFile;
        $outTmpPrefix =~ s/\.bam$//i;
        if ( $by eq "name" ) {
            $exit = system(
"samtools sort -n -o $outFile -T $outTmpPrefix -\@ $ncores $inFile"
            );
        }
        elsif ( $by eq "pos" ) {
            $exit = system(
                "samtools sort -o $outFile -T $outTmpPrefix -\@ $ncores $inFile"
            );
        }
        if ( $exit != 0 ) {
            $LOGGER->fatal(
                "Cannot sort $inFile by $by and save it to $outFile!")
              && die();
        }
    }

    sub postproc {
        my $self = shift;
        $LOGGER->info("Sorting the output BAM by position...\n");
        my $posSortedFile = $outFile . ".posSortedTmp.bam";
        $self->sort_bam( $outFile, $posSortedFile, "pos", $ncores );
        File::Copy::move( $posSortedFile, $outFile );
        if ($doesOutDup) {
            my $posSortedDupFile = $outDupFile . ".posSortedTmp.bam";
            $self->sort_bam( $outDupFile, $posSortedDupFile, "pos", $ncores );
            File::Copy::move( $posSortedDupFile, $outDupFile );
        }
    }

    sub sort_name {
        my $self           = shift;
        my $nameSortedFile = $outFile . ".nameSortedTmp.bam";
        $self->sort_bam( $outFile, $nameSortedFile, "name", $ncores );
        File::Copy::move( $nameSortedFile, $outFile );
    }

    sub output_stats {
        my $self = shift;
        ReadPair->output_stats();
    }

    sub output_cnts {
        my $self = shift;
        if ( !$doesOutBed ) { return 0; }
        $LOGGER->info("Writing focusTag counts into a BED file...\n");
        my @strands;
        if ($ignoreStrand) {
            @strands = qw( * );
        }
        else {
            @strands = qw( + - );
        }
        $LOGGER->debug("strands = @strands\n");
        for my $strand (@strands) {
            for my $i ( 0 .. $#SeqLevels ) {
                my $seqLevel  = $SeqLevels[$i];
                my $seqLength = $SeqLengths[$i];
                my $nBins     = POSIX::ceil( $seqLength / $mindist );
                $LOGGER->debug(
"Writing counts to BED: seqLevel = $seqLevel, strand = $strand, nBins = $nBins\n"
                );
                for my $binIdx ( 0 .. ( $nBins - 1 ) ) {
                    my $c = $cnts{$strand}{$seqLevel}{$binIdx};
                    $c = defined($c) ? $c : 0;
                    if ( $doesOutBedZeroCount || ( $c > 0 ) ) {
                        my $start = $mindist * $binIdx;
                        my $end   = $mindist * ( $binIdx + 1 ) - 1;
                        $end = $end > $seqLength ? $seqLength : $end;
                        print $outBedFh
                          "$seqLevel\t$start\t$end\t*\t$c\t$strand\n";
                    }
                }
            }
        }
    }
}

my $samReader = SamReader->new(
    inFile              => $inFile,
    outFile             => $outFile,
    outDupFile          => $outDupFile,
    outBedFile          => $outBedFile,
    outStatsFile        => $outStatsFile,
    doesOutBedZeroCount => $doesOutBedZeroCount,
    nameSorted          => $nameSorted,
    priority            => $priority,
    endType             => $endType,
    scheme              => $scheme,
    as_per_mate         => $as_per_mate,
    nm_per_mate         => $nm_per_mate,
    penalty             => $penalty,
    flags               => $flags,
    ncores              => $ncores,
    ncells              => $ncells,
    mindist             => $mindist,
    ignoreStrand        => $ignoreStrand,
    addSoftClipped      => $addSoftClipped,
    focusTag            => $focusTag,
    doesOutNoFocus      => $doesOutNoFocus,
);
$samReader->check_args();
$samReader->parse_qualfocus();
$samReader->init_fh();
$samReader->init_class();
$samReader->init_seqinfo();
$samReader->init_cnts();
$samReader->slurp_reads();
$samReader->sort_quals();
$samReader->detect_dup();
$samReader->output_cnts();
$samReader->fin_fh();
$samReader->postproc();
$samReader->output_stats();
$LOGGER->info("All done.\n");
