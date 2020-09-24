#### Arguments to be passed ####
# sys.argv[1] = file of blast output
import sys
import argparse
import subprocess
import os
import datetime
import pickle
import Bio
from Bio.SubsMat import MatrixInfo

def run_blast(fasta, blastdb, db, outputdir, name):
    ''' Runs blastp
        
        fasta: path fasta containing sequences to blast
        db: path to database to blast against
        outputdir: path to directory in which to write blast output
        name: sample name to distinguish file
    '''

    outfile = outputdir + "/" + name + "." + "blast.out"
    if os.path.isfile(outfile) == False:
        subprocess.call([blastp, "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", db, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", outfile])
    else:
        print("Blast file for " + name + " already exists - skipping")

def score_match(pair, matrix):
    ''' Gives a score from a matrix for a given pair of sequences

        pair: pair to score (tuple)
        matrix: scoring matrix (matrix)

        Return value: score
    '''
    # If the pair is not in the matrix, reverse it and get score
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
        # If the pair is in the matrix, get score
    else:
        return matrix[pair]

def score_pairwise(seq1, seq2, matrix):
    ''' For two sequences of equal length, scores them at non-anchor positions given a matrix
        seq1: first sequence
        seq2: second sequence
        matrix: scoring matrix
    '''
    score = 0
    last_ep = len(seq1) -1
    # Add score to running total for all non-anchor positions
    for i in range(len(seq1)):
        if i != 1 and i != last_ep:
            pair = (seq1[i], seq1[i])
            score += score_match(pair, matrix)
        return score

def process_blast(blast_results, matrix):
    ''' Processes results of blastp to obtain dictionary of best results
        Keys are neoepitopes sequences
        Values are lists of associated blast data: E value, sequence of match, raw protein similarity score
        Transcript and gene of match peptide are also stored
        blast_results: path to file containing results of blastp
        matrix: scoring matrix to use for comparing sequences

        return value: dictionary
    '''
    dict = pickle.load(open("humanDict.pickle", "rb"))

    blast_dict = {}
    with open(blast_results, "r") as fh:
        for line in fh:
            # Obtain relevant data - epitope sequence, alignment length, E value, matching sequence, peptide name
            line = line.strip("\t").split("\t")
            epitope = line[0]#.split("=")[1]
            #print(epitope)
            length = int(line[2])
            eval = float(line[6])
            match_seq = line[5]
            match_pep = line[1].replace("ref", "").replace("|", "")
            match_transcript = dict[match_pep][1]
            match_gene = dict[match_pep][0]
            # Check for presence of invalid characters in match seq
            invalids = ["B", "J", "O", "U", "X", "Z", "*"]
            invalid_matches = []
            for char in invalids: 
                if char in match_seq:
                    invalid_matches.append(char)

            # Check if epitope is not already in dictionary, the alignment is the right length, and there are no invalid characters
            if epitope not in blast_dict and length == len(epitope) and invalid_matches == []:
                # Find BLOSUM score between neoepitope and peptide match, then store data
                match_ps = score_pairwise(epitope, match_seq, matrix)
                blast_dict[epitope] = [epitope, eval, match_transcript, match_gene, match_seq, match_ps]
            # If epitope is in dictionary, but E value for this entry is better, replace data
            elif epitope in blast_dict and length == len(epitope) and eval < blast_dict[epitope][1] and invalid_matches == []:
                # Fine BLOSUM score between neoepitope and peptide match, then store data
                match_ps = score_pairwise(epitope, match_seq, matrix)
                blast_dict[epitope] = [epitope, eval, match_transcript, match_gene, match_seq, match_ps]
            # If epitope is in dictionary and E value for this entry is equivalent, compare further
            elif epitope in blast_dict and length == len(epitope) and eval == blast_dict[epitope][1] and invalid_matches == []:
                # If the match sequence is te same as previous entry, skip
                if match_seq == blast_dict[epitope][4] and match_transcript not in blast_dict[epitope][2] and match_gene not in blast_dict[epitope][3]:
                    blast_dict[epitope][2] = blast_dict[epitope][2] + "," + match_transcript
                    blast_dict[epitope][3] = blast_dict[epitope][3] + "," + match_gene
                # If the match sequence is a better match to neoepitope, replace data
                else: 
                    match_ps = score_pairwise(epitope, match_seq, matrix)
                    if match_seq == epitope or match_ps > blast_dict[epitope][5]:
                        blast_dict[epitope] = [epitope, eval, match_transcript, match_gene, match_seq, match_ps]
    return blast_dict

def make_output(dict_file, name, outputdir):
    outfile = outputdir + "/" + name + ".epitopes.annotated.tsv"
    with open(outfile, "w") as out:
        # Write header to outfile
        out.write("New_seq\te-value\tMatch_ps\n")
        for x in dict_file:
            seq = dict_file[x][0]
            new_seq = dict_file[x][4]
            e = dict_file[x][1]
            match_ps = dict_file[x][5]
            out.write(seq + "\t" + new_seq + "\t" + str(e) + "\t" + str(match_ps) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input', type=str, required=True,
    #        help='path input epitope data tsv file'
    #    )
    parser.add_argument('-o', '--outdir', type=str, required=True, 
            help='path to output directory'
        )
    parser.add_argument('-s', '--sample', type=str, required=True, 
            help='sample name'
        )
    args = parser.parse_args()

    print('Timestamp: {:%Y-%m-%d %H:%m:%S}'.format(datetime.datetime.now()) + " Running sequence_similariy.py for " + args.sample)

    # Set BLOSUM62 as blosum
    blosum = MatrixInfo.blosum62

    # Set paths to blast databases
    humanDB = "/home/eknodel/blastdb/Homo_sapiens.GRCh38.pep.all.fa"

    # Produce fasta file containing neoepitopes
    fasta_path = "/data/CEM/shared/public_data/cohen_melanoma/netMHC/patient_3466/test.in"
    #fasta_path = args.outdir + "/" + args.sample + "peptides.fasta"

    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Running blast now...")

    blastp = "blastp"
    # Run blast comparing neoepitopes to human peptides
    run_blast(fasta_path, blastp, humanDB, args.outdir, args.sample)

    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Processing blast results...")

    # Process blast data and save to dictionaries
    hum_dict = process_blast(args.outdir+"/"+args.sample+".blast.out", blosum)
    #print(hum_dict["MMSWKTRGK"])
    make_output(hum_dict, args.sample, args.outdir)
