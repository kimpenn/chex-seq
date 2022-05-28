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
use Carp;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Log::Log4perl::Layout::PatternLayout;

my ( $inFile, $outFile );
my $endType             = "PE";
my $flags               = "";
my $qualities           = "A,B";
my $outputInPair        = "read";
my $doesOutputSingleton = 0;
my $version             = 0;
my $debug               = "info";
my $help                = 0;
our $VERSION = "0.31";
our $LOGGER  = get_logger(__PACKAGE__);

{

    package SamReader;

    my ( $endType, @quals, $outputInPair, $doesOutputSingleton );

    sub uniq {
        my %s;
        return grep { !$s{$_}++ } @_;
    }

    sub initClass {
        my $class = shift;
        ( $endType, $qualities, $outputInPair, $doesOutputSingleton ) = @_;
        @quals = split( ',', $qualities );
        @quals = map {
                ( $_ eq 'A' || $_ eq 'B' || $_ eq 'C' )
              ? ( "${_}1", "${_}2" )
              : ($_)
        } @quals;
        @quals = &uniq(@quals);
        $LOGGER->debug("\@quals = (@quals)\n");
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

    sub init {
        my $self = shift;
        $LOGGER->info("Opening file handles...\n");
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

    sub getPrimerQual {
        my ($readName) = @_;
        my @f          = split( ':', $readName, -1 );
        my $focusTag   = pop @f;
        my $qual       = pop @f;
        return ( $qual, $focusTag );
    }

    sub getFocusInPair {
        my ($focusTag) = @_;
        return "R1"  if $endType eq "SE";
        return "R$1" if ( $focusTag =~ /F([\d])/ );
        return undef;
    }

    sub getReadInPair {
        my ($flag) = @_;
        return "R1" if $endType eq "SE";
        return $flag & 0x40 ? "R1" : "R2";
    }

    sub isMateMapped {
        my ($flag) = @_;
        return undef if $endType eq "SE";
        return !( $flag & 0x08 );
    }

    sub iter {
        my $self = shift;
        $LOGGER->info("Iterating records...\n");
        my ( $inFh, $outFh ) = ( $self->{inFh}, $self->{outFh} );
        while (<$inFh>) {
            my $rec = $_;
            if (/^@/) {
                print $outFh $_;
                next;
            }
            my @f = split( "\t", $rec );
            my ( $readName, $flag )     = ( $f[0], $f[1] );
            my ( $qual,     $focusTag ) = &getPrimerQual($readName);
            my $focusInPair  = &getFocusInPair($focusTag);
            my $readInPair   = &getReadInPair($flag);
            my $isMateMapped = &isMateMapped($flag);
            for my $q (@quals) {
                next if $qual ne $q;
                if ( $endType eq "PE" ) {
                    if ( $q ne 'D' ) {
                        if ( !$doesOutputSingleton ) {
                            ( ( print $outFh $rec ) && next )
                              if $outputInPair eq "read"
                              && $focusInPair eq $readInPair;
                            ( ( print $outFh $rec ) && next )
                              if $outputInPair eq "mate"
                              && $focusInPair ne $readInPair;
                            ( ( print $outFh $rec ) && next )
                              if $outputInPair eq "both" && $isMateMapped;
                            ( ( print $outFh $rec ) && next )
                              if $outputInPair eq "either";
                        }
                        else {
                            ( ( print $outFh $rec ) && next )
                              if !$isMateMapped
                              && $outputInPair eq "read"
                              && $focusInPair eq $readInPair;
                            ( ( print $outFh $rec ) && next )
                              if !$isMateMapped
                              && $outputInPair eq "mate"
                              && $focusInPair ne $readInPair;
                            ( ( print $outFh $rec ) && next )
                              if !$isMateMapped && $outputInPair eq "either";
                        }
                    }
                    else
                    { # D reads don't have focus tag, so read/mate does not apply to D. That is, none to be output given qual is D and outputInPair is read or mate.
                        if ( !$doesOutputSingleton ) {
                            ( ( print $outFh $rec ) && next )
                              if $isMateMapped && $outputInPair eq "both";
                            ( ( print $outFh $rec ) && next )
                              if $outputInPair eq "either";
                        }
                        else {
                            ( ( print $outFh $rec ) && next )
                              if !$isMateMapped
                              && ( $outputInPair eq "either" );
                        }
                    }
                }
                else {    # SE
                    ( ( print $outFh $rec ) && next )
                      if $focusInPair eq $readInPair;
                }
            }
        }
    }

    sub fin {
        my $self = shift;
        $LOGGER->info("Closing file handles...\n");
        close $self->{inFh}
          or $LOGGER->warn("Cannot close filehandle for read!\n");
        close $self->{outFh}
          or $LOGGER->warn("Cannot close filehandle for write!\n");
    }

}

sub usage {
    print STDERR <<DOC;
Summary: 
    Filter a SAM/BAM file with annotated barcode/primer quality in the read names for the specified quality class (--qual) and read/mate criteria (--outputInPair). 

Usage:
    perl $0 --inFile Sample_XXX.star.primaryNoDup.PrimerAnnotated.bam --outFile Sample_XXX.star.primaryNoDup.A1read.bam [--endType PE] [--flags ''] [--qualities A,B] [--outputInPair read] [--version] [--help] [--debug info]

Options:
    --inFile, -i            input file, BAM/SAM format;
    --outFile, -o           output file, BAM/SAM format;
    --endType               choose from 'PE' or 'SE' (default: 'PE');
    --flags                 additional flags to be passed to samtools, e.g. '-f 0x100' (default: '');
    --qualities             the quality category (or categories, comma separated) to output. Options: A1, A2, B1, B2, C1, C2, D, A (= A1 or A2), B (= B1 or B2), C (= C1 or C2) (default: 'A,B');
    --outputInPair          for PE only. If the insert meets given quality criteria, which in a read pair to output. Options: read (the focus-carrying read only), mate (the mate of the focus-carrying read only, namely, the nonfocus-carrying read only), both (both the focus- and nonfocus-carrying read; both have to be mapped), either (either the focus- or nonfocus-carrying read, whichever is mapped) (default: read). Note, since D reads have neither focus nor nonfocus read, so read and mate do not apply to D reads. In other reads, nothing will be output if quality is D and outputInPair is read or mate;
    --doesOutputSingleton   for PE only. If an insert meets given quality criteria, output this read if its mate is not mapped (i.e. a singleton) (default: false). Note, since this option contradicts --outputInPair both, nothing will be output in this case. 
DOC
}

GetOptions(
    "inFile|i=s"          => \$inFile,
    "outFile|o=s"         => \$outFile,
    "endType=s"           => \$endType,
    "flags=s"             => \$flags,
    "qualities=s"         => \$qualities,
    "outputInPair=s"      => \$outputInPair,
    "doesOutputSingleton" => \$doesOutputSingleton,
    "help|h"              => \$help,
    "version|v"           => \$version,
    "debug=s"             => \$debug,
) or ( &usage() && exit(-1) );

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

$LOGGER->info(
"{ inFile = $inFile, outFile = $outFile, endType = $endType, flags = $flags, qualities = $qualities, outputInPair = $outputInPair, doesOutputSingleton = $doesOutputSingleton, debug = $debug, VERSION = $VERSION }\n"
);

SamReader->initClass( $endType, $qualities, $outputInPair,
    $doesOutputSingleton );
my $samReader = SamReader->new( $inFile, $outFile, $flags );
$samReader->init();
$samReader->iter();
$samReader->fin();
$LOGGER->info("All done.\n");
