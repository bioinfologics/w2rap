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
    if line[0] != ">": # and scaffold_id is not None:
      sline = line.strip().split("\t")
      temp_list.append({"contig_id": int(sline[0]), "contig_pos_initial": int(sline[1]),  "orientation": sline[2], "contig_pos_final": int(sline[3])})

    else:
      contig_pos_in_scaff[scaffold_id] = temp_list
      scaffold_id = line[1:].strip()
      temp_list = []

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

  scaffold_list = list(scaffold)
  for cd in contigs_data:
    remapped_Ns = 0
    #print cd

    sequence = contig_seq[cd["contig_id"]]
    posicion_inicial = cd["contig_pos_initial"]
    posicion_final = cd["contig_pos_final"]
    orientation = cd["orientation"]

    if orientation == '+':
      #print "orientacion +"
      contig_position = len(sequence)-1
      scaffold_position = posicion_inicial+posicion_final-1
      while scaffold_position > posicion_inicial:
        if sequence[contig_position] == "N":
          scaffold_list[scaffold_position] = "N"
          remapped_Ns += 1
        contig_position -= 1
        scaffold_position -= 1

    elif orientation == '-':
      #print "orientacion -"
      scaffold_position = posicion_inicial+posicion_final-1
      contig_position = 0
      while scaffold_position > posicion_inicial: 
        if sequence[contig_position] == "N":
          scaffold_list[scaffold_position] = "N"
          remapped_Ns += 1
        scaffold_position -= 1
        contig_position += 1

    else:
      print "non ident orientation"
      raise 

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
  modificados = 0
  n_modificados = 0
  for scaff in scaffolds_complete_set:
    if scaff in scaff_index_set and scaff is not None:
      #print contigs_in_scaff
      modificados += 1
      contigs_in_scaff = scaff_index[scaff]
      contig_seq = {contig["contig_id"]: contigs[contig["contig_id"]] for contig in contigs_in_scaff}
      Nmaped_scaffold = align_contigs(scaffolds[scaff], contigs_in_scaff, contig_seq)
      output_scaffolds.write(">%s\n%s\n" %(scaff, Nmaped_scaffold))

    elif scaff is not None:
      n_modificados += 1
      output_scaffolds.write(">%s\n%s\n" %(scaff, scaffolds[scaff]))
    else:
      pass

    if nscaffolds % 1000 == 0:
      print "Scaffolds processed: ", nscaffolds, ",", modificados,",",n_modificados
    nscaffolds += 1
  #raw_input()
  output_scaffolds.close()

