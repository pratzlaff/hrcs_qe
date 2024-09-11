#! /usr/bin/perl -w
use strict;

use Chandra::Tools::Common qw(read_bintbl_cols);
use PDL;
use PDL::Graphics::PGPLOT;

my $f = @ARGV ? shift : 'hrcsD1999-07-22qeN0012.fits';

my ($r1, $e1, $qe1) = read_bintbl_cols($f, qw( regionid energy qe ), { extname => 'axaf_qe1' });
my ($r2, $e2, $qe2) = read_bintbl_cols($f, qw( regionid energy qe ), { extname => 'axaf_qe2' });
my ($r3, $e3, $qe3) = read_bintbl_cols($f, qw( regionid energy qe ), { extname => 'axaf_qe3' });

dev '/xs', 1, 2;
line $e1->slice(',(0)'), $qe1->slice(',(0)'), { color => 'orange' };
hold;
line $e1->slice(',(1)'), $qe1->slice(',(1)'), { color => 'blue' };
release;

if (
  sum(abs($e1->slice(',(0)') - $e3->slice(',(0)'))) == 0 and
  sum(abs($e1->slice(',(1)') - $e3->slice(',(1)'))) == 0
  ) {
  print "QE chips 1 and 3 ARE identical\n";
}
else {
  print "QE chips 1 and 3 ARE NOT identical\n";
}

line $e2->slice(',(0)'), $qe2->slice(',(0)');
hold;
line $e2->slice(',(1)'), $qe2->slice(',(1)'), { color => 'red' };
line $e2->slice(',(2)'), $qe2->slice(',(2)'), { color => 'blue' };
line $e2->slice(',(3)'), $qe2->slice(',(3)'), { color => 'orange' };

exit 0;
