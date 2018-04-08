use strict;
use warnings;
use 5.010;

my @names =("j9030_9.sm", "j9041_4.sm", "j909_8.sm");

my $n_1 = 1000; 
my $n_2 = 0;
my $n_3 = 3000;
my $n_4 = 0;
my $n_5 = 0;
my $conc = "";
my $cont=0;
foreach my $n (@names){
       for($n_3 = 50000;$n_3<=250000;$n_3+=50000){
	for($cont = 0; $cont<=15 ;$cont++){
	  for($n_2 = 4; $n_2<=16; $n_2*=2){
	    for($n_4 = 10 ; $n_4 <= 100; $n_4 += 90){	
		$conc = $n."_s";
		$conc .= $cont;
		$conc .=".txt";
		if($n_2 == 4){
		   $n_5 = 3;
		   system("./md $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>novox2/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		   $n_5 = 2;
		   system("./md $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>novox2/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		} 
		if($n_2 == 8){
		   $n_5 = 4;
		   system("./md $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>novox2/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");		
		   $n_5 = 6;
		   system("./md $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>novox2/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		}
		 
		if($n_2 == 16){
		   $n_5 = 10;
		   system("./md $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>novox2/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");		
		   $n_5 = 12;
		   system("./md $conc $n_1 $n_2 $n_3 $n_4 1 $n_5 >>novox2/Resultado_$conc"."_$n_1"."_$n_2"."_$n_3"."_$n_4"."_$n_5.txt");	
		}
	    }
	  }
	}
	}
}

exit;

