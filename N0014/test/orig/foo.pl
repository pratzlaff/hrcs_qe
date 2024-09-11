#! /usr/bin/perl

use strict;
use warnings;

use PDL;
use PDL::Graphics::PGPLOT;
use Chandra::Tools::Common 'read_bintbl_cols';

my ($e,$qe)=read_bintbl_cols(
			     '/data/legs/rpete/flight/hrcs_qe/hrcsD2012-03-29qeN0013.fits',
			     qw/ energy qe /, { extname => "AXAF_QE1" }
			    );

my ($e2,$qe2)=read_bintbl_cols(
			     '/usr/local/ciao/CALDB/data/chandra/hrc/qe/hrcsD2012-03-29qeN0013.fits',
			     qw/ energy qe /, { extname => "AXAF_QE1" }
			    );

print "Min max caldb/N0014 QE ratios:";
print join(", ", minmax($qe2/$qe)),"\n";


$_ = $_->slice(',(0)')->sever for $e, $qe;

my $ediff=$e->slice("1:-1") - $e->slice("0:-2");
print "$_\n" for $e->dims;
dev "qe_ediff.ps/cps", { hardlw => 2, hardch => 1.5 };
line $e->slice('1:-1'), $ediff, 
  {
   border => 1,
   xtitle => 'Energy (keV)',
   ytitle => 'Bin Width (keV)'
   };



