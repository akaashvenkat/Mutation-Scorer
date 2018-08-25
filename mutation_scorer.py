#!/usr/bin/env python

# File:             mutation_scorer.py
# Authors:          Akaash Venkat, Yevgeniy Sazhnyev, Professor Jie Zheng
# Description:      Calculates change in free energy (DDG) of mutation of proteins.

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.toolbox import *
from time import localtime, gmtime, strftime
import os, struct, sys, time

AMINO_START_POSITION = 243
FILE_OUTPUTS = []
TIMESTAMP = strftime("%d%b%Y_%H.%M.%S_", localtime())

INITIAL_RELAX_PDB_SCORE = 0
RELAX_POSE = 0
RELAX_RELAX = 0
RELAX_SCOREFXN = 0
RELAX_OUTPUT_TEXT_FILE = TIMESTAMP + "output_relax.txt"

INITIAL_DDG_PDB_SCORE = 0
DDG_POSE = 0
DDG_RELAX = 0
DDG_SCOREFXN = 0
DDG_OUTPUT_TEXT_FILE = TIMESTAMP + "output_ddg.txt"


def main():

    relax_output_file = open(RELAX_OUTPUT_TEXT_FILE, "w+")
    ddg_output_file = open(DDG_OUTPUT_TEXT_FILE, "w+")

    mutation_relax_list = []
    mutation_relax_score = []
    mutation_relax_str_score = []
    diff_relax_str_score = []
    
    mutation_ddg_list = []
    mutation_ddg_score = []
    mutation_ddg_str_score = []
    diff_ddg_str_score = []

    os.system('clear')
    print("------------------------")
    print("How many mutations would you like to calculate score functions for?")
    print("------------------------")
    cnt = input("Number of Mutations:   ")
    count = int(float(cnt))

    print("------------------------")
    print("For each mutation, enter Amino Acid Position, followed by Amino Acid, as one word.")
    print("------------------------")
    
    for counter in range(0, count):
        temp_index = str(counter + 1)
        temp_line = input("Mutation #" + temp_index + ":   ")
        mutation_relax_list.append(temp_line)
        
    os.system('clear')

    initRelaxPyrosetta()
    run_pos = 1
    for counter in range(0, count):
        if (counter > 0):
            if (mutation_relax_list[counter] != mutation_relax_list[counter - 1]):
                run_pos = 1
                mutation_relax_score.append(calculateRelaxScore(mutation_relax_list[counter], run_pos))
            else:
                run_pos = run_pos + 1
                mutation_relax_score.append(calculateRelaxScore(mutation_relax_list[counter], run_pos))
        else:
            mutation_relax_score.append(calculateRelaxScore(mutation_relax_list[counter], run_pos))

    initDDGPyrosetta()
    for counter in range(0, len(FILE_OUTPUTS)):
        mutation_ddg_score.append(calculateDDGScore(FILE_OUTPUTS[counter]))

    os.system('clear')

    relax_output_file.write("------------------------" + "\n")
    relax_output_file.write("Initial Relax PDB Score: " + str(INITIAL_RELAX_PDB_SCORE) + "\n")
    relax_output_file.write(" " + "\n")

    for y in range(0, len(mutation_relax_list)):
        mutation_relax_list[y] = "Mutation " + str(y + 1) + ": " + mutation_relax_list[y]
        mutation_relax_str_score.append("Relax Score: " + str(mutation_relax_score[y]))
        diff_relax_str_score.append("(Difference: " + str(mutation_relax_score[y] - INITIAL_RELAX_PDB_SCORE) + ")")

    print_list = []
    print_list.append(mutation_relax_list)
    print_list.append(mutation_relax_str_score)
    print_list.append(diff_relax_str_score)
    row_format = '{:<35}' * len(print_list)
    for v in zip(*print_list):
        relax_output_file.write(row_format.format(*v) + "\n")
    relax_output_file.write("------------------------" + "\n")

    ddg_output_file.write("------------------------" + "\n")
    ddg_output_file.write("Initial DDG PDB Score: " + str(INITIAL_DDG_PDB_SCORE) + "\n")
    ddg_output_file.write(" " + "\n")

    for y in range(0, len(mutation_ddg_score)):
        mutation_ddg_list.append(mutation_relax_list[y])
        mutation_ddg_str_score.append("DDG Score: " + str(mutation_ddg_score[y]))
        diff_ddg_str_score.append("(Difference: " + str(mutation_ddg_score[y] - INITIAL_DDG_PDB_SCORE) + ")")

    print_list = []
    print_list.append(mutation_ddg_list)
    print_list.append(mutation_ddg_str_score)
    print_list.append(diff_ddg_str_score)
    row_format = '{:<35}' * len(print_list)
    for v in zip(*print_list):
        ddg_output_file.write(row_format.format(*v) + "\n")
    ddg_output_file.write("------------------------" + "\n")

    relax_output_file.close()
    ddg_output_file.close()

    calculateRelaxAveDiff()
    calculateDDGAveDiff()
    print("Relax and DDG PDB scores and their average differences can be found in output .txt files.")




def initRelaxPyrosetta():
    
    global INITIAL_RELAX_PDB_SCORE
    global RELAX_POSE
    global RELAX_RELAX
    global RELAX_SCOREFXN

    init(extra_options = "-beta_nov16_cart -in:file:s 4wxs.clean.pdb -packing:use_input_sc -relax:constrain_relax_to_start_coords -in:ignore_unrecognized_res -relax:coord_constrain_sidechains -relax:ramp_constraints false -relax:cartesian -relax:min_type lbfgs_armijo_nonmonotone")
    
    RELAX_POSE = Pose()
    RELAX_POSE = pose_from_pdb("4wxs.clean.pdb")
    RELAX_SCOREFXN = ScoreFunction()
    RELAX_SCOREFXN = get_fa_scorefxn()
    RELAX_SCOREFXN(RELAX_POSE)
    
    RELAX_RELAX = pyrosetta.rosetta.protocols.relax.FastRelax()
    RELAX_RELAX.constrain_relax_to_start_coords(True)
    RELAX_RELAX.coord_constrain_sidechains(True)
    RELAX_RELAX.ramp_down_constraints(False)
    RELAX_RELAX.cartesian(True)
    RELAX_RELAX.min_type("lbfgs_armijo_nonmonotone")
    RELAX_RELAX.set_scorefxn(RELAX_SCOREFXN)
    RELAX_RELAX.apply(RELAX_POSE)
    
    RELAX_POSE.dump_pdb(TIMESTAMP + "4wxs.relaxed.pdb")
    INITIAL_RELAX_PDB_SCORE = RELAX_SCOREFXN(RELAX_POSE)



def calculateRelaxScore(mutation_string, mutation_no):

    global RELAX_POSE
    global RELAX_RELAX
    global RELAX_SCOREFXN

    num_amino_acids = 0
    amino_pos = []
    amino_acid = []

    parsed_list = mutation_string.split(' ')
    for x in range(0, len(parsed_list)):
        y = parsed_list[x]
        amino_pos.append(int(y[0:len(y)-1]))
        amino_acid.append(y[-1:])
        num_amino_acids = num_amino_acids + 1

    RELAX_POSE = pose_from_pdb("4wxs.clean.pdb")
    for x in range(0, num_amino_acids): 
        mutate_residue(RELAX_POSE, (amino_pos[x] - AMINO_START_POSITION), amino_acid[x])

    RELAX_RELAX.apply(RELAX_POSE)
    file_output = TIMESTAMP + mutation_string[:4] + "_" + str(mutation_no) + ".pdb"
    FILE_OUTPUTS.append(file_output)
    RELAX_POSE.dump_pdb(file_output)

    return RELAX_SCOREFXN(RELAX_POSE)



def calculateRelaxAveDiff():
    
    file_content = []
    with open(RELAX_OUTPUT_TEXT_FILE) as file:
        file_content = file.readlines()
    file_content = [valid_content[:-2] for valid_content in file_content]
    for counter in range (0, 3):
        file_content.pop(0)
    file_content.pop(len(file_content)-1)

    simplified_list = []
    for counter in range (0, len(file_content)):
        col1 = ""
        if counter < 9:
            col1 = file_content[counter][12:16]
        elif counter < 99 and counter >= 10:
            col1 = file_content[counter][13:17]
        else:
            col1 = file_content[counter][14:18]
        col2 = file_content[counter][83:]
        col2 = col2.replace(' ', '')
        col2 = col2.replace(')', '')
        
        row = []
        row.append(col1)
        row.append(col2)
        simplified_list.append(row)

    mut_count = []
    mut_count.append(1)
    for counter in range(1, len(simplified_list)):
        if (simplified_list[counter][0] == simplified_list[counter - 1][0]):
            mut_count[len(mut_count) - 1] = mut_count[len(mut_count) - 1] + 1
        else:
            mut_count.append(1)

    final_output_mut = []
    final_output_diff = []
    parse_pos = 0
    for x in range(0, len(mut_count)):
        score_sum = 0
        for y in range(parse_pos, parse_pos + mut_count[x]):
            score_sum = score_sum + float(simplified_list[y][1])
        final_output_mut.append(simplified_list[parse_pos][0])
        score_ave = score_sum / mut_count[x]
        final_output_diff.append(score_ave)
        parse_pos = parse_pos + mut_count[x]

    for x in range(0, len(final_output_mut)):
        final_output_mut[x] = "Mutation: " + final_output_mut[x]
        final_output_diff[x] = "Average Relax PDB Score Difference: " + str(final_output_diff[x])
    
    final_out = []
    final_out.append(final_output_mut)
    final_out.append(final_output_diff)
    
    relax_output_file = open(RELAX_OUTPUT_TEXT_FILE, "a+")
    row_format = '{:<35}' * len(final_out)
    relax_output_file.write("")
    for v in zip(*final_out):
        relax_output_file.write(row_format.format(*v) + "\n")
    relax_output_file.write("------------------------" + "\n")
    relax_output_file.close()



def initDDGPyrosetta():
    
    global INITIAL_DDG_PDB_SCORE
    global DDG_POSE
    global DDG_RELAX
    global DDG_SCOREFXN
    
    output_pdb = TIMESTAMP + "4wxs.relaxed.pdb"
    
    init(extra_options = "-beta_nov16_cart -in:file:s output_pdb -packing:use_input_sc -relax:constrain_relax_to_start_coords -in:ignore_unrecognized_res -relax:coord_constrain_sidechains -relax:ramp_constraints false -relax:cartesian -relax:min_type lbfgs_armijo_nonmonotone -ddg:mut_file -ddg:iterations 5 -max_cycles 200 -relax:min_type lbfgs_armijo_nonmonotone -fa_max_dis 9.0")
    
    DDG_POSE = Pose()
    DDG_POSE = pose_from_pdb(output_pdb)
    DDG_SCOREFXN = ScoreFunction()
    DDG_SCOREFXN = get_fa_scorefxn()
    DDG_SCOREFXN(DDG_POSE)
    
    DDG_RELAX = pyrosetta.rosetta.protocols.relax.FastRelax()
    DDG_RELAX.constrain_relax_to_start_coords(True)
    DDG_RELAX.coord_constrain_sidechains(True)
    DDG_RELAX.ramp_down_constraints(False)
    DDG_RELAX.cartesian(True)
    DDG_RELAX.min_type("lbfgs_armijo_nonmonotone")
    DDG_RELAX.set_scorefxn(DDG_SCOREFXN)
    INITIAL_DDG_PDB_SCORE = DDG_SCOREFXN(DDG_POSE)



def calculateDDGScore(input_file):
    
    global DDG_POSE
    global DDG_RELAX
    global DDG_SCOREFXN
    
    init(extra_options = "-beta_nov16_cart -in:file:s input_file -packing:use_input_sc -relax:constrain_relax_to_start_coords -in:ignore_unrecognized_res -relax:coord_constrain_sidechains -relax:ramp_constraints false -relax:cartesian -relax:min_type lbfgs_armijo_nonmonotone -ddg:mut_file -ddg:iterations 5 -max_cycles 200 -relax:min_type lbfgs_armijo_nonmonotone -fa_max_dis 9.0")
    
    DDG_POSE = Pose()
    DDG_POSE = pose_from_pdb(input_file)
    DDG_SCOREFXN = ScoreFunction()
    DDG_SCOREFXN = get_fa_scorefxn()
    DDG_SCOREFXN(DDG_POSE)
    
    DDG_RELAX = pyrosetta.rosetta.protocols.relax.FastRelax()
    DDG_RELAX.constrain_relax_to_start_coords(True)
    DDG_RELAX.coord_constrain_sidechains(True)
    DDG_RELAX.ramp_down_constraints(False)
    DDG_RELAX.cartesian(True)
    DDG_RELAX.min_type("lbfgs_armijo_nonmonotone")
    DDG_RELAX.set_scorefxn(DDG_SCOREFXN)
    return DDG_SCOREFXN(DDG_POSE)



def calculateDDGAveDiff():

    file_content = []
    with open(DDG_OUTPUT_TEXT_FILE) as file:
        file_content = file.readlines()
    file_content = [valid_content[:-2] for valid_content in file_content]
    for counter in range (0, 3):
        file_content.pop(0)
    file_content.pop(len(file_content)-1)
    
    simplified_list = []
    for counter in range (0, len(file_content)):
        col1 = ""
        if counter < 9:
            col1 = file_content[counter][12:16]
        elif counter < 99 and counter >= 10:
            col1 = file_content[counter][13:17]
        else:
            col1 = file_content[counter][14:18]
        col2 = file_content[counter][83:]
        col2 = col2.replace(' ', '')
        col2 = col2.replace(')', '')
        
        row = []
        row.append(col1)
        row.append(col2)
        simplified_list.append(row)
    
    mut_count = []
    mut_count.append(1)
    for counter in range(1, len(simplified_list)):
        if (simplified_list[counter][0] == simplified_list[counter - 1][0]):
            mut_count[len(mut_count) - 1] = mut_count[len(mut_count) - 1] + 1
        else:
            mut_count.append(1)

    final_output_mut = []
    final_output_diff = []
    parse_pos = 0
    for x in range(0, len(mut_count)):
        score_sum = 0
        for y in range(parse_pos, parse_pos + mut_count[x]):
            score_sum = score_sum + float(simplified_list[y][1])
        final_output_mut.append(simplified_list[parse_pos][0])
        score_ave = score_sum / mut_count[x]
        final_output_diff.append(score_ave)
        parse_pos = parse_pos + mut_count[x]

    for x in range(0, len(final_output_mut)):
        final_output_mut[x] = "Mutation: " + final_output_mut[x]
        final_output_diff[x] = "Average DDG PDB Score Difference: " + str(final_output_diff[x])
    
    final_out = []
    final_out.append(final_output_mut)
    final_out.append(final_output_diff)
    
    ddg_output_file = open(DDG_OUTPUT_TEXT_FILE, "a+")
    row_format = '{:<35}' * len(final_out)
    ddg_output_file.write("")
    for v in zip(*final_out):
        ddg_output_file.write(row_format.format(*v) + "\n")
    ddg_output_file.write("------------------------" + "\n")
    ddg_output_file.close()



if __name__ == '__main__':
    main()
