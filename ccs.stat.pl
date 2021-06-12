use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Statistics::Basic qw(mean stddev);
GetOptions(
    "code"      => \ my $code,  # use CodecV1 instead of frame for calculation. Default option is to convert the code to frame number.
    "sam=s"     => \ my $sam,
    "window=i"  => \(my $win = 10),
    "out=s"     => \ my $out
) or die "Invalid options passed to $0\n";


if($sam=~/\.bam$/) {
    open IN, "samtools view -\@5 $sam|";
}
elsif($sam=~/\.sam$/) {
    open IN, $sam;
}
else {
    die "$sam is no a bam/sam file!\n";
}

my @c_ipd = ();
my @g_ipd = ();
my @a_ipd = ();
my @t_ipd = ();
my @c_pw = ();
my @g_pw = ();
my @a_pw = ();
my @t_pw = ();
my @rev_c_ipd = ();
my @rev_g_ipd = ();
my @rev_a_ipd = ();
my @rev_t_ipd = ();
my @rev_c_pw = ();
my @rev_g_pw = ();
my @rev_a_pw = ();
my @rev_t_pw = ();
while(<IN>)
{
    chomp;
    next if /^@/;
    my $line=$_;
    my @F=split /\t/, $_, 12;
    my ($query, $flag, $ref, $mapPos, $mapQ, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = @F;
    my %tags=split /:\S:|\s+/, $F[11];

    my $rg = $tags{"RG"};
    my $ec = $tags{"ec"};  #  Effective coverage
    my $fn = $tags{"fn"};  # Forward number of complete passes
    my $np = $tags{"np"};  # number of passes.
    my $rn = $tags{"rn"};  # Reverse number of complete passes
    my $rq = $tags{"rq"};  # Float in [0, 1] encoding expected accuracy
    my $zm = $tags{"zm"};
    my @fi = split /,/, $tags{"fi"};
    my @fp = split /,/, $tags{"fp"};
    my @ri = split /,/, $tags{"ri"};
    my @rp = split /,/, $tags{"rp"};
    my @sn = split /,/, $tags{"sn"};  # 4 floats for the average signal-to-noise ratio of A, C, G, and T (in that order) over the HQRegion
    shift @fi;
    shift @fp;
    shift @ri;
    shift @rp;
    unless ($code) {
        @fi = map {load_frame_codes($_)} @fi;
        @fp = map {load_frame_codes($_)} @fp;
        @ri = map {load_frame_codes($_)} @ri;
        @rp = map {load_frame_codes($_)} @rp;
    }
    @ri = reverse @ri;
    @rp = reverse @rp;
    shift @sn; 
    my @c = indexes($seq, "C");
    my @g = indexes($seq, "G");
    my @a = indexes($seq, "A");
    my @t = indexes($seq, "T");
    # print Dumper(@cpgSites);
    push @c_ipd, @fi[@c];
    push @c_ipd, @fi[@c];
    push @g_ipd, @fi[@g];
    push @a_ipd, @fi[@a];
    push @t_ipd, @fi[@t];
    push @c_pw, @fp[@c];
    push @g_pw, @fp[@g];
    push @a_pw, @fp[@a];
    push @t_pw, @fp[@t];
    push @rev_c_ipd, @ri[@c];
    push @rev_g_ipd, @ri[@g];
    push @rev_a_ipd, @ri[@a];
    push @rev_t_ipd, @ri[@t];
    push @rev_c_pw, @rp[@c];
    push @rev_g_pw, @rp[@g];
    push @rev_a_pw, @rp[@a];
    push @rev_t_pw, @rp[@t];

    # foreach my $i(@cpgSites) {
    #     my @fi_c = 
    #     my $subfi = join ",", @fi[$start..$end];
    #     my $subfp = join ",", @fp[$start..$end];
    #     my $subri = join ",", @ri[$start..$end];
    #     my $subrp = join ",", @rp[$start..$end];
    #     # print Dumper($substring);
    #     print join("\t", $query, $start, $end, "CG", $length, $np, $substring, $subfi, $subfp, $subri, $subrp), "\n";
    # }
}
close IN;
my $mean_c_ipd =  mean(@c_ipd);
my $mean_g_ipd =  mean(@g_ipd);
my $mean_a_ipd =  mean(@a_ipd);
my $mean_t_ipd =  mean(@t_ipd);
my $mean_c_pw =  mean(@c_pw);
my $mean_g_pw =  mean(@g_pw);
my $mean_a_pw =  mean(@a_pw);
my $mean_t_pw =  mean(@t_pw);
my $mean_c_rev_ipd =  mean(@rev_c_ipd);
my $mean_g_rev_ipd =  mean(@rev_g_ipd);
my $mean_a_rev_ipd =  mean(@rev_a_ipd);
my $mean_t_rev_ipd =  mean(@rev_t_ipd);
my $mean_c_rev_pw =  mean(@rev_c_pw);
my $mean_g_rev_pw =  mean(@rev_g_pw);
my $mean_a_rev_pw =  mean(@rev_a_pw);
my $mean_t_rev_pw =  mean(@rev_t_pw);

my $std_c_ipd =  stddev(@c_ipd);
my $std_g_ipd =  stddev(@g_ipd);
my $std_a_ipd =  stddev(@a_ipd);
my $std_t_ipd =  stddev(@t_ipd);
my $std_c_pw =  stddev(@c_pw);
my $std_g_pw =  stddev(@g_pw);
my $std_a_pw =  stddev(@a_pw);
my $std_t_pw =  stddev(@t_pw);
my $std_c_rev_ipd =  stddev(@rev_c_ipd);
my $std_g_rev_ipd =  stddev(@rev_g_ipd);
my $std_a_rev_ipd =  stddev(@rev_a_ipd);
my $std_t_rev_ipd =  stddev(@rev_t_ipd);
my $std_c_rev_pw =  stddev(@rev_c_pw);
my $std_g_rev_pw =  stddev(@rev_g_pw);
my $std_a_rev_pw =  stddev(@rev_a_pw);
my $std_t_rev_pw =  stddev(@rev_t_pw);

print $mean_c_ipd, ":", $std_c_ipd,"\t", $mean_g_ipd, ":", $std_g_ipd,"\t", $mean_a_ipd, ":", $std_a_ipd,"\t", $mean_t_ipd, ":", $std_t_ipd,"\n";
print $mean_c_pw,":",$std_c_pw,"\t", $mean_g_pw,":",$std_g_pw,"\t", $mean_a_pw,":",$std_a_pw,"\t", $mean_t_pw,":",$std_t_pw, "\n";
print $mean_c_rev_ipd,":",$std_c_rev_ipd,"\t", $mean_g_rev_ipd,":",$std_g_rev_ipd,"\t", $mean_a_rev_ipd,":",$std_a_rev_ipd,"\t", $mean_t_rev_ipd,":",$std_t_rev_ipd,"\n";
print $mean_c_rev_pw,":",$std_c_rev_pw,"\t", $mean_g_rev_pw,":",$std_g_rev_pw,"\t", $mean_a_rev_pw,":",$std_a_rev_pw,"\t", $mean_t_rev_pw,":",$std_t_rev_pw, "\n";

sub indexes {
    my $string = shift;
    my $sub = shift;
    my $offset = 0;
    my @out = ();
    my $index = index($string, $sub, $offset);
    while ($index != -1) {
        push @out, $index;
        $offset = $index + 1;
        $index = index($string, $sub, $offset);
    }
    return @out;
}

sub load_frame_codes {
    my $code = shift;
    my $frame;
    if ($code < 64) {
        $frame = $code;
    }
    elsif ($code < 128) {
        $frame = 64 + ($code - 64)*2;
    }
    elsif ($code < 192) {
        $frame = 192 + ($code - 128)*4;
    }
    elsif ($code < 255) {
        $frame = 448 + ($code - 192)*8;
    }
    else {
        die "ip or pw frame code wrong";
    }
    return $frame;
}

