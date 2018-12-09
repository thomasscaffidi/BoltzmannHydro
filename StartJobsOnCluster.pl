print "Hello \n";

my $cluster="lr2";
my $n=1;

my $DirName = $ARGV[0];
my $BatchFileName = "BaseForScript_cluster_${cluster}_n_${n}.sh";

system("mkdir ".$DirName);
chdir($DirName);
system("cp ../*.m .");
system("cp ../" . $BatchFileName ." .");

$filename = "BatchScript_${DirName}.sh";
system("cp " . $BatchFileName ." ".  $filename);
open(my $fh, '>>', $filename) or die "Could not open file '$filename' $!";
print $fh "\n module load matlab ";
print $fh "\n /global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/matlab -nodesktop -nosplash -r 'GenerateGoodData_SAS'";
system("sbatch ".$filename);
