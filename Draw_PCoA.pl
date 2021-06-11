#! /usr/bin/perl
=head1 Description: plot PCoA 

=head1 Author: huangyf@genomics.cn

=head1 Version

 Version 1.0: 2020-05-19

=head1 Usage

 perl Draw_PCoA.pl 

 [normal parameters]
 -map	[s]	Sample group information file
 -file	[s]	Species, gene or function relatvie abundance file, support .gz file. Beta diversity will be calculate by vegan package.
 -beta	[s]	Beta diversity distance matrices file
 -name	[s]	Name of group, or head name of column, only use the given groups to draw. For example: Group1:Group2,Description
 -pc	[s]	Coordinate display [default 1:2]
 -method[s]	Beta distrance method, [default bray]
 -prefix[s]	Output prefix [default All]
 -outdir[s]	Output directory [default ./]
 -resdir[s]	Copy result to directory
 -Rpath	[s]	R software pathway [default /hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1]
 
 [plot parameters]
 -draw_sa	Display the sample names or not, and do not recommend displaying if there are many samples
 -ellipse	Display the inertia ellipse
 -cex_sa[f]	Size of sample name [default 2]

=head1 Example

 perl Draw_PCoA.pl -map Mapping.txt -file species_profile.txt -name Control:Case,Description
 perl Draw_PCoA.pl -map Mapping.txt -beta unweighted_unifrac_Beta_diversity.txt -name Control:Case,Description

=head1 Note

 1. The accepted distances methods are: manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao or mahalanobis.

=cut

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Cwd;

my ($map,$file,$name,$ellipse,$draw_sa,$result_dir,$beta_file);
my ($Prefix,$outdir,$cex_sa,$Rpath,$pc,$method) = ("All","./",2,"/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1","1:2","bray");

GetOptions (
	'map:s'		=>	\$map,
	'file:s'	=>	\$file,
	'beta:s'	=>	\$beta_file,
	'name:s'	=>	\$name,
	'pc:s'		=>	\$pc,
	'method:s'	=>	\$method,
	'prefix:s'	=>	\$Prefix,
	'outdir:s'	=>	\$outdir,
	'resdir:s'	=>	\$result_dir,
	'draw_sa'	=>	\$draw_sa,
	'ellipse'	=>	\$ellipse,
	'cex_sa:f'	=>	\$cex_sa,
	'Rpath:s'	=>	\$Rpath
);

die(`pod2text $0`) if (!defined $map && !defined $name);
if (!defined $file && !defined $beta_file)
{
	die "Please provide beta distance file or relative abundance file!\n";
}
if (defined $file && defined $beta_file)
{
	die "-file and -beta cannot appear at the same time.\n";
}

#----------------- make directory -----------------
my $pwd = cwd();
$outdir = ($outdir eq "./" || $outdir eq ".") ? $pwd : ($outdir =~ /^\//) ? $outdir : $pwd."/".$outdir;
if (defined $result_dir)
{
	$result_dir = ($result_dir eq ".") ? $result_dir : ($result_dir =~ /^\//) ? $result_dir : $pwd."/".$result_dir;
}

$beta_file = "$pwd/$beta_file" if ($beta_file !~ /^\//);

#----------------- run subroutine -----------------
my @group_info;
if($name =~ /,/)
{
	$name =~ s/\s+//g;
	@group_info = split /\,/, $name;
} else {
	@group_info = $name;
}

open MAP, "$map";
chomp (my $map_head = <MAP>);
my @map_head = split /\t+/, $map_head;
close MAP;
	
my $group_flag = 0;
foreach my $group_info (@group_info)
{
	foreach my $map_head (@map_head)
	{
		if ($map_head eq $group_info) {
			$group_flag = 1;
			&plot ($group_flag, $group_info);
		}
	} 
	if ($group_info =~ /:/)
	{
		$group_flag = 2;
		&plot ($group_flag, $group_info);
	}
}

#----------------- Subroutine for drawing plot -----------------
sub plot {
my ($Group_flag, $Group_info) = @_;

open MAP, "$map";
my (%sample, @Map_head,@group);
if ($Group_flag == 1)
{
	chomp (my $Map_head = <MAP>);
	@Map_head = split /\t/, $Map_head;
} elsif ($Group_flag == 2)
{
	@group = ($Group_info =~ /:/) ? (split /\:/,$Group_info) : $Group_info;
}

while (<MAP>)
{
	chomp;
	next if (/^#/);
	my @array = split /\t/;
	my $flag = 0;
	for (my $i = 1; $i <= $#array; $i++)
	{
		if ($Group_flag == 2)
		{
			foreach my $group (@group)
			{
				if ($group eq $array[$i])
				{
					$sample{$array[0]} = $array[$i];
					$flag = 1;
				}
			}
		} elsif ($Group_flag == 1)
		{
			if ($Map_head[$i] eq $Group_info)
			{
				$sample{$array[0]} = $array[$i];
				push @group, $array[$i];
				$flag = 1
			}
		}
		last if ($flag == 1);
	}
}
close MAP;

if ($Group_flag == 1)
{
	my %hash;
	@group = grep {++$hash{$_} < 2} @group;
	@group = sort @group;
}
my $join_groups = (@group < 10) ? join "-",@group : join "-",$group[0], $group[-1],scalar(@group)."Groups";
system ("mkdir -m 755 -p $outdir/$join_groups") unless (-d "$outdir/$join_groups");
my $outdir = "$outdir/$join_groups";
my $prefix = "$Prefix\_$join_groups";
my $Result_dir;
if (defined $result_dir)
{
	system ("mkdir -m 755 -p $result_dir/$join_groups") unless (-d "$outdir/$join_groups");
	$Result_dir = "$result_dir/$join_groups";
}

open GR, ">$outdir/$prefix.group.xls";
print GR "SampleID\tgroupname\n";
foreach my $sample (sort keys %sample)
{
	print GR "$sample\t$sample{$sample}\n";
}
close GR;

#----------------- input profile file -----------------
my @pc = split /\:/, $pc;
my $dist_code = "";
if (defined $beta_file)
{
	$dist_code = "data <- read.table('$beta_file',header=TRUE,check.names=FALSE,row.names=1,sep='\t')
group <- read.table('$outdir/$prefix.group.xls',head=TRUE,check.name=FALSE,row.names=1,sep='\t')
group.names <- rownames(group)
beta.names <- colnames(data)
intersection <- intersect(group.names,beta.names)
beta.dist <- data[pmatch(intersection,beta.names),pmatch(intersection,beta.names)]
beta.dist <- as.dist(as.matrix(beta.dist))
group\$SampleID <- group.names";
} elsif (defined $file)
{
	($file =~ /\.gz$/) ? open IN, "gzip -dc $file |" : open IN, "$file";
	open OUT, ">$outdir/$prefix.profile.xls";
	chomp (my $head = <IN>);
	my @head = split /\t/, $head;
	my $line = 0;
	while (<IN>)
	{
		chomp;
		next if (/^#/);
		$line ++;
		my @array = split /\t/;
		my ($print,$total,$flag,$sample_title);
		for(my $i=1; $i<=$#array; $i++)
		{
			if (exists $sample{$head[$i]})
			{
				$sample_title .= "\t$head[$i]";
				$print .= "\t$array[$i]";
				$total += $array[$i];
				$flag = 1;
			}
		}
		print OUT "$sample_title\n" if ($line == 1);
		print OUT "$array[0]"."$print\n" if ($total >0 && $flag == 1);
	}
	close IN;
	close OUT;
	
	$dist_code = "data <- t(read.table('$outdir/$prefix.profile.xls',header=TRUE,check.name=FALSE))
group <- read.table('$outdir/$prefix.group.xls',head=TRUE,check.name=FALSE)
beta.dist <- vegdist(data, method = \"$method\")";
}

my $pcoa_code = "
PCOA <- pcoa(beta.dist, correction='none', rn=NULL)
PC.eig <- as.numeric(sprintf ('%0.2f',PCOA\$values[,'Relative_eig']*100))
x = PCOA\$vectors
sample_names = rownames(x)
pc <- as.data.frame(PCOA\$vectors)
pc\$SampleID = sample_names 
Merge.result <- merge(pc,group,by='SampleID',all=TRUE)
write.table(Merge.result,file='$outdir/$prefix.PC_coordinate.xls',quote=FALSE,sep='\t',row.names=F,append=F)
write.table(PC.eig,file='$outdir/$prefix.PC_eig.xls',quote=FALSE,sep='\t')
data.plot <- read.table(\"$outdir/$prefix.PC_coordinate.xls\",head=T)
PC <- read.table(\"$outdir/$prefix.PC_eig.xls\",head=T)
xlab=paste(\"PC$pc[0] (\",PC[$pc[0],1],\"%)\",sep='')
ylab=paste(\"PC$pc[1] (\",PC[$pc[1],1],\"%)\",sep='') 
";

my ($groupname) = ("",0,"");
foreach my $group (@group)
{
	$groupname .= "\"$group\",";
}

my $group_num = scalar (@group);

chop ($groupname);
my $draw_sa_code = (defined $draw_sa) ? "\n\tgeom_text_repel(data = data.plot,aes(x=Axis.$pc[0],y=Axis.$pc[1],label=SampleID),size=$cex_sa,segment.size=0.2,segment.color=\"grey\") +" : "";
my $ellipse_code = (defined $ellipse) ? "\n\tstat_ellipse(data = data.plot,aes(x=Axis.$pc[0],y=Axis.$pc[1],color=groupname),show.legend = FALSE) +" : "";

my $width = ($group_num%14==0) ? 1*($group_num/14) : 1*(int($group_num/14)+1);
my $legend_col = ($group_num%14==0) ? 1*($group_num/14) : 1*(int($group_num/14)+1);

my $Rscript = <<R;
library("ggplot2",lib.loc="$Rpath/library/")
library("RColorBrewer",lib.loc="$Rpath/library")
library("vegan",lib.loc="$Rpath/library/")
library("ggrepel",lib.loc="$Rpath/library/")
library("ggpubr",lib.loc="$Rpath/library/")
library('cowplot',lib.loc='$Rpath/library')
library("ape",lib.loc='$Rpath/library')

$dist_code
$pcoa_code
data.plot\$groupname <- factor(data.plot\$groupname,levels=c($groupname))
plot <- ggplot() + 
	geom_hline(yintercept=0,color="grey",linetype=4) +
	geom_vline(xintercept=0,color="grey",linetype=4) +
	geom_point(data = data.plot,aes(x=Axis.$pc[0],y=Axis.$pc[1],color=groupname),size=1,alpha=0.6) + $draw_sa_code $ellipse_code
	labs(x=xlab, y=ylab,color="",fill="") +
	scale_color_manual(values=colorRampPalette(brewer.pal(9, "Set1"))($group_num)) +
	theme(axis.text = element_text(color = "black",size = 8),
		axis.title = element_text(color = "black",size = 8),
		axis.ticks = element_line(color = "black"),
		axis.line = element_line(color = "black"),
#		legend.position = c(0,1),
#		legend.justification = c(0,1),
		legend.text = element_text(size=8),
		legend.key.size = unit(0.1,"in"),
		legend.key = element_blank(),
		legend.background = element_blank(),
		legend.spacing.x = unit(0.01,"in"),
		panel.border = element_rect(colour = "black", fill=NA),
		panel.grid = element_blank(),
		panel.background = element_blank()
#		panel.background = element_rect(colour = "black")
	)
pdf("$outdir/$prefix.PCoA_plot.pdf",width = 2.2 + $width,height= 2)
plot_grid(plot_grid(plot + theme(legend.position ="none")),plot_grid(get_legend(plot+ guides(col = guide_legend(ncol=$legend_col)) + theme(legend.position = c(0,0.6),legend.justification = 'left'))),rel_widths=c(2.2,$width))
png("$outdir/$prefix.PCoA_plot.png",width = 2.2 + $width,height= 2, units='in',res=300)
plot_grid(plot_grid(plot + theme(legend.position ="none")),plot_grid(get_legend(plot+ guides(col = guide_legend(ncol=$legend_col)) + theme(legend.position = c(0,0.6),legend.justification = 'left'))),rel_widths=c(2.2,$width))
R

open RS, ">$outdir/$prefix.PCoA_plot.R";
print RS "$Rscript";
close RS;
`source /hwfssz4/BC_PUB/Software/03.Soft_ALL/SourceMe2.sh; $Rpath/bin/R -f $outdir/$prefix.PCoA_plot.R`;

if (defined $result_dir)
{
	`cp $outdir/$prefix.PCoA_plot.pdf $Result_dir`;
	`cp $outdir/$prefix.PCoA_plot.png $Result_dir`;
	`cp $outdir/$prefix.PC_coordinate.xls $Result_dir`;
}
}
