###################################################
#                                                 #
#           Bioinformatic challenges              #
#             Rodrigo Álvarez Pardo               #
#                 Assignment 3                    #
#                                                 #
#                                                 #
###################################################

# Requirements:

require './Program.rb'

####################################################

# Global variables:

$target = "cttctt"                # Global variable contains the sequence to search in the exons
$length_target = $target.length() # Global variable with the target's length

#####################################################

puts "Working on the task...\n"
#We create the files that will contain:
#the positions referred to the genes and the chromosomes,
#in addition to the file that will contain the genes that do not have exonic sequences
$gff_genes = create_open_file("genes.gff3")
$gff_chr = create_open_file("chromosomes.gff3")
$no_targets = create_open_file("genes_without_target.txt")

######################################################

# We print the headers of each file
$gff_genes.puts "##gff-version 3"
$gff_chr.puts "##gff-version 3"
$no_targets.puts "Genes without #{$target.upcase} in exons\n\n"


genes_ids = load_gene_list("GeneList.txt") 

#######################################################



genes_ids.each do |gene|                          # Part of the program to search the targets in the exons
  
  seq_obj = gene_information(gene)                # Create the Bio:EMBL object from each gene
  
  unless seq_obj.nil?
    target_hash = get_exon_target(seq_obj)        # We get the targets inside exons of the gene
    
    
    if target_hash.empty? 
      $no_targets.puts gene                       # If the gene has no targets in exons, we add it to the file genes_without_target.txt
     
      
    else  
      add_features(gene, target_hash, seq_obj)    # We create new features and add them to each seq_obj
      chr = get_chromosome(gene, seq_obj)         # We return the chromosome number and postions
      convert_chromosome(gene, target_hash, chr)  # We convert the positions to the ones that correspond in the chromosome
    
    end
    
  end
  
end

###################################################

puts "You can chechk the output in the files: "
puts "\t- genes.gff3"
puts "\t- chromosomes.gff3"
puts "\t- genes_without_target.txt"
