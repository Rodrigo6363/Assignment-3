###################################################
#                                                 #
#           Bioinformatic challenges              #
#             Rodrigo Álvarez Pardo               #
#                 Assignment 3                    #
#                                                 #
#                                                 #
###################################################

# Requirements:

require 'bio'
require 'net/http'
require 'rest-client'


###################################################

def load_gene_list(filename)    # Method to read the genes from the file
  
  ld = File.open(filename, "r")
  genes = Array.new
  
  ld.each_line do |line|
    genes << line.delete("\n") 
  end
  
  ld.close
  return genes.uniq             # We apply uniq preventing posible duplications 

end

#######################################################



def gene_information(gene_id)   # Function to get the information about each gene
  
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&id=#{gene_id}&format=embl&sid=raw&Retrieve=Retrieve"
  
  response = RestClient::Request.execute(
  method: :get,
  url: address)
  
  return nil unless response
  record = response.body
  data = Bio::EMBL.new(record)
  bioseq = data.to_biosequence
  
  return bioseq

end
  
#######################################################

def find_target_exon(exon_id,target_seq,length_seq,exon_position,strand) # Method for obtaining exonic sequences in both positive and negative orientation
  
  target_exon = Hash.new # Hash that contains the exon id whit the mathches positions (init and end)
  
  case strand
    
  when '+' # Foward
    
    target_seq.each do |match_start|
      
      match_end = match_start + $length_target - 1
      
      if (match_start >= exon_position[0].to_i) && (match_start <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
        target_exon[[match_start, match_end]] = [exon_id, '+']
      end
    
    end
    
  when '-' # Reverse
    
    target_seq.each do |match_start|
      
      match_end = match_start + $length_target - 1
      
      if (match_start >= exon_position[0].to_i) && (match_start <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
        
        # To work will whit forward strand, we need to convert the positions as follows
        m_end = length_seq - match_end
        m_start = length_seq - match_start
        target_exon [[m_end, m_start]] = [exon_id, '-']
      
      end
      
    end
    
  end
  
  if not target_exon.empty?
    return target_exon
  end
  
end
        
######################################################

def get_exon_target (bio_seq_obj) # Method that given a Bio::EMBL object returns a hash whit the coordinates of target´s matches inside exons
  
  length_bio_seq = bio_seq_obj.length() # Length of nucleotide sequence
  target_positions_exon = Hash.new      # Hash that will contain the positions targeted inside exons
  
  # We get the target matches in both directions
  
  target_match_foward = bio_seq_obj.gsub(/#{$target}/).map{Regexp.last_match.begin(0)}
  target_match_reverse = bio_seq_obj.reverse.tr('atgc','tacg').gsub(/#{$target}/).map{Regexp.last_match.begin(0)}
  
  bio_seq_obj.features.each do |feature|
    
    position = feature.position
    
    next unless (feature.feature == 'exon' && (not position =~ /[A-Z]/)) # We look for the feature type exon
    
    exon_id = feature.qualifiers[0].value.gsub('exon_id =', '')
    
    if position =~ /complement/ # Exon in reverse strand
      
      position = position.tr('complement()', '').split('..')
      position_reverse = []
      
      # Getting a 2 elements array containing initial and end position
      
      position.each do |pst|
        position_reverse.insert(0, length_bio_seq - pst.to_i)
      end
      
      target_pst_exon = find_target_exon(exon_id, target_match_reverse, length_bio_seq, position_reverse, '-')
      # We call the method "find_target_in_exon" to determine which matches are inside of the exon.
      if not target_pst_exon.nil? # If we don't get a response we put the objectives in a hash
        target_positions_exon = target_positions_exon.merge(target_pst_exon)
      end
      
    else # Exon in forward strand
      
      position = position.split('..') 
      
      target_pst_exon = find_target_exon(exon_id, target_match_reverse, length_bio_seq, position, '+')
      # We call the method "find_target_in_exon" to determine which matches are inside of the exon.
      if not target_pst_exon.nil? # If we don't get a response we put the objectives in a hash
        target_positions_exon = target_positions_exon.merge(target_pst_exon)
      end
      
    end
    
  end
  
  return target_positions_exon
  
end

#######################################################

def add_features(gene_id, targets, bioseq) # Method that iterates over the target´s matches exons to add them as new features to the Bio: EMBL objects.
  
  exon_features = Array.new
  
  targets.each do |target, exon_id_strand|
    
    feat = Bio::Feature.new("#{$target.upcase}", "#{target[0]}...#{target[1]}")
    
    feat.append(Bio::Feature::Qualifier.new('nucleotide_motif', "#{$target.upcase}_in_#{exon_id_strand[0]}"))
    
    # New feature call Nucleotide motif that contains a region of nucleotide sequence corresponding to a known motif
    
    feat.append(Bio::Feature::Qualifier.new('strand', exon_id_strand[1]))
    
    $gff_genes.puts "#{gene_id}\t.\t#{feat.feature}\t#{target[0]}\t#{target[1]}\t.\t#{exon_id_strand[1]}\t.\tID=#{exon_id_strand[0]}"
    
    # Prints the feature in the GFF3 gene file 
    
    exon_features << feat
    
  end
  
  bioseq.features.concat(exon_features) # We add the new features created to the existing ones

  
end

####################################################

def get_chromosome(gene_id, bio_seq_obj) # Routine that given a Bio:Sequence object returns the chromosome and positions to which the sequence belongs
  
  bsopa = bio_seq_obj.primary_accession
  
  return false unless bsopa
  
  chromosome_array = bsopa.split(":")
  
  $gff_chr.puts "#{chromosome_array[2]}\t.\tgene\t#{chromosome_array[3]}\t#{chromosome_array[4]}\t.\t+\t.\tID=#{gene_id}"
  
  # Prints the information of the gene in the GFF
  
  return chromosome_array[2], chromosome_array[3], chromosome_array[4] # Return the chromosome number, gene start position and gene end position
end

####################################################

def create_open_file(file) # Method that checks whether a given file exists, in which case it deletes it and opens it.
  
  if File.exists?(file)
    File.delete(file)
  end
  
  return File.open(file, "a+")

end

######################################################

def convert_chromosome(gene,targets,chr) # This method translates the coordinates to the ones referring to the chromosome and it prints them on the GFF3 chromosome fil

  targets.each do |positions, exon_strand|
    position_start = chr[1].to_i + positions[0].to_i
    position_end = chr[1].to_i + positions[1].to_i
  
    $gff_chr.puts "#{chr[0]}\t.\tnucleotide_motif\t#{position_start }\t#{position_end}\t.\t#{exon_strand[1]}\t.\tID=#{exon_strand[0]};parent=#{gene}"
  
  end

#######################################################

end
