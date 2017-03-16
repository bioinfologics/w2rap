#!/usr/bin/env python

from Bio import SeqIO
import sys

def load_contig_pos_in_scaff(filename):
  """ {scaffold: (contig_id:, contig_pos_initial:, contig_pos_final:, orientation:)} """

  contig_pos_in_scaff = {}
   
  scaffold_id = None
  file = open(filename, "r")
  temp_list = []
  for line in file:
    line = line.strip()
    if line[0] == ">":
      if scaffold_id is not None:
        # we've hot a new one so add previous to dict
        contig_pos_in_scaff[scaffold_id] = temp_list
        temp_list = []
      scaffold_id = line[1:].strip()

    else:
      sline = line.strip().split("\t")
      temp_list.append({"contig_id": int(sline[0]), "contig_pos_initial": int(sline[1]),  "orientation": sline[2], "contig_pos_final": int(sline[3])})

    # add the final one
    contig_pos_in_scaff[scaffold_id] = temp_list

  return contig_pos_in_scaff

def get_scaffolds(scaffolds_file):
  """ """
  scaffolds = {scaff.id: str(scaff.seq) for scaff in SeqIO.parse(scaffolds_file, 'fasta')}
  return scaffolds

def get_contigs(contigs_file):
  """ """
  contigs = {int(contig.id): str(contig.seq) for contig in SeqIO.parse(contigs_file, 'fasta')}
  return contigs

def align_contigs(scaffold, contigs_data, contigs_seq):
  """ Align contigs to the scaffolds to remap Ns """

  #print "scaffold:", scaffold
  #print "contigs_data:", contigs_data
  #print "contigs_seq:", contigs_seq

  scaffold_list = list(scaffold)
  for cd in contigs_data:
    remapped_Ns = 0
    #print cd

    sequence = contig_seq[cd["contig_id"]]
    pos_initial = cd["contig_pos_initial"]
    pos_final = cd["contig_pos_final"]
    orientation = cd["orientation"]

    if orientation == '+':
      #print "orientacion +"
      contig_position = len(sequence)-1
      scaffold_position = pos_initial + pos_final - 1
      while scaffold_position > pos_initial:
        if sequence[contig_position] == "N":
          scaffold_list[scaffold_position] = "N"
          remapped_Ns += 1
        contig_position -= 1
        scaffold_position -= 1

    elif orientation == '-':
      #print "orientacion -"
      contig_position = 0
      scaffold_position = pos_initial + pos_final - 1
      while scaffold_position > pos_initial: 
        if sequence[contig_position] == "N":
          scaffold_list[scaffold_position] = "N"
          remapped_Ns += 1
        scaffold_position -= 1
        contig_position += 1

  return "".join(scaffold_list)

if __name__ == "__main__":

  if len(sys.argv) == 1:
    print "Usage: {0} contig_pos_in_scaff scaffolds_file contigs_file scaffolds_output".format(sys.argv[0])
    sys.exit()

  scaff_index = load_contig_pos_in_scaff(sys.argv[1])
  print "References loaded: %s" %(len(scaff_index.keys()), )

  scaffolds = get_scaffolds(sys.argv[2])
  print "Scaffolds loaded: %s" %(len(scaffolds.keys()), )

  contigs = get_contigs(sys.argv[3])
  print "Contigs loaded: %s" %(len(contigs.keys()), )
  
  scaff_index_set = set(scaff_index.keys())
  scaffolds_complete_set = set(scaffolds.keys())

  output_scaffolds = open(sys.argv[4], "w")
  nscaffolds = 0
  modified = 0
  not_modified = 0
  for scaff in scaffolds_complete_set:
    if scaff in scaff_index_set and scaff is not None:
      modified += 1
      contigs_in_scaff = scaff_index[scaff]
      contig_seq = {contig["contig_id"]: contigs[contig["contig_id"]] for contig in contigs_in_scaff}
      nmapped_scaffold = align_contigs(scaffolds[scaff], contigs_in_scaff, contig_seq)
      output_scaffolds.write(">%s\n%s\n" %(scaff, nmapped_scaffold))

    elif scaff is not None:
      # these are sequences in the final scaffolds that don't appear in the contigPosInscaff file, ie. singleton contigs
      # so write the contig from the .contigs file rather than the one from .scafSeq which has the C/G problem
      # contig >123 in .contig will be sequence >C123 in .scafSeq
      contig_id = int(scaff[1:])
      not_modified += 1
      
      output_scaffolds.write(">%s\n%s\n" %(scaff, contigs[contig_id]))

    else:
      pass

    if nscaffolds % 1000 == 0 and nscaffolds > 0:
      print "Scaffolds processed: {0}, ({1} modified, {2} not modified)".format(nscaffolds, modified, not_modified)
    nscaffolds += 1

  print "Completed - Scaffolds processed: {0}, ({1} modified, {2} not modified)".format(nscaffolds, modified, not_modified)

  output_scaffolds.close()

