use strict;
use warnings;
use 5.010;


my $filename = 'insta2.txt';
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Não foi possível abrir o arquivo '$filename' $!";
 
while (my $row = <$fh>) {
 system("./m5 $row");
}

exit;

