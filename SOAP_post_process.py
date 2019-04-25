#!/usr/bin/env python

from Bio import SeqIO
import sys
import re
import tempfile
import os

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

    sequence = contigs_seq[cd["contig_id"]]
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

def remap_ns(contig_pos_in_scaff, scaf_seq, contigs, fd):

  print "** Re-mapping N-stretches **"

  scaff_index = load_contig_pos_in_scaff(contig_pos_in_scaff)
  print "References loaded: %s" %(len(scaff_index.keys()), )

  scaffolds = get_scaffolds(scaf_seq)
  print "Scaffolds loaded: %s" %(len(scaffolds.keys()), )

  contigs = get_contigs(contigs)
  print "Contigs loaded: %s" %(len(contigs.keys()), )

  scaff_index_set = set(scaff_index.keys())
  scaffolds_complete_set = set(scaffolds.keys())

  output_scaffolds = os.fdopen(fd, "w")
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

  print "DONE - Scaffolds processed: {0}, ({1} modified, {2} not modified)\n".format(nscaffolds, modified, not_modified)

  output_scaffolds.close()

def collapse(before_gap_seq, after_gap_seq, gap_length):
  """ Work out whether to collapse a sequence over a gap """

  # rtn = 0 : don't collapse (default)
  # rtn = 1 : collapse

  rtn = 0

  #print "gap_length=", gap_length

  if int(gap_length) < 200:
    #print "before", before_gap_seq
    #print "after", after_gap_seq

    repeat_start = after_gap_seq[0:30]
    #print "seq to look for before gap",repeat_start

    found_before_gap = before_gap_seq.find(repeat_start)
    if found_before_gap > -1:
      #print "found at",found_before_gap
      repeat_seq_before = before_gap_seq[found_before_gap:]
      #print "before",repeat_seq_before

      repeat_seq_after = after_gap_seq[:len(repeat_seq_before)]
      #print "after",repeat_seq_after

      if repeat_seq_before == repeat_seq_after and len(repeat_seq_before) < 200:
        #print "repeat_length=",len(repeat_seq_before)
        rtn = 1

  return rtn


def collapse_gaps(tmp_file, output):
  """ Process a FASTA fie and collapse gaps generated by the SOAP scaffolder where we have the same seqeunce before and after the gap.
  To be collapsed, the gap must be <200bp and the perfect overlap is <=199bp """

  print "** Collapsing repeats around gaps **"

  seq_count = 0
  collapse_count = 0
  not_collapse_count = 0

  # open output file
  fout = open(output, 'w')

  seqiter = SeqIO.parse(open(tmp_file), 'fasta')
  for seq in seqiter:
    #print "checking", seq.id, "length", len(seq.seq)

    seq_count = seq_count + 1
    new_seq = ""
    prev_gap_end = 0

    # find gaps and get start and end co-ords
    p = re.compile("N+")
    for m in p.finditer(str(seq.seq)):
      #print "start=", m.start(), "end=", m.end()
      gap_start = m.start()
      gap_end = m.end()

      #print "first N at", gap_start + 1
      #print "last N at", gap_end

      gap_length = int(gap_end) - int(gap_start)

      # get 200 bases before and after the gap
      before_gap_seq = seq.seq[gap_start - 200:gap_start - 1]
      after_gap_seq = seq.seq[gap_end:gap_end + 200]
      if collapse(before_gap_seq, after_gap_seq, gap_length) == 1:	# collapse
        # record seq from end of prev gap to start of current gap (which includes the collapsed repeat)
        new_seq = new_seq + seq.seq[prev_gap_end:gap_start]
        collapse_count = collapse_count + 1
      else:	# don\t collapse
        # record seq from end of prev gap to end of current gap
        new_seq = new_seq + seq.seq[prev_gap_end:gap_end]
        not_collapse_count = not_collapse_count + 1

      # record the prev gap end
      prev_gap_end = gap_end

    # add the sequence after the final gap
    new_seq = new_seq + seq.seq[prev_gap_end:]

    # write the new seq to a file
    fout.write(">{0}\n{1}\n".format(seq.id, new_seq))

  fout.close

  print "DONE - {0} sequences processed, {1} collapsed, {2} not collapsed".format(seq_count, collapse_count, not_collapse_count)

if __name__ == "__main__":

  if len(sys.argv) != 5:
    print "Usage: {0} contig_pos_in_scaff scaffolds_file contigs_file scaffolds_output".format(sys.argv[0])
    sys.exit()

  # get temp file for intermediate file
  (fd, filename) = tempfile.mkstemp()
  remap_ns(sys.argv[1], sys.argv[2], sys.argv[3], fd)

  collapse_gaps(filename, sys.argv[4])

  # delete temp file
  os.remove(filename)

