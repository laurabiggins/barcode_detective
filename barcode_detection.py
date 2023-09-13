#!/usr/bin/env python
import re
import argparse

# TODO: check for Ns

def main():
    options = read_options()

    # read in codes
    all_codes = read_codes(options.found_codes)
    
    #setup output dictionaries
    [observed_dict,mia_dict] = create_output_dicts(all_codes)

    #check for correct matches
    observed_dict = check_correct_codes(observed_dict)
    
    # check for polyG
    observed_dict = checkPolyG(observed_dict)

    # check for PhiX
    observed_dict = checkPhiX(observed_dict)

    # print(observed_dict)
    # print("=====================")
    # print(mia_dict)

    for observed in observed_dict:
        if observed_dict[observed]["explained"] == False:
            info_message = None
            obs_barcode = observed_dict[observed]["barcode_1"]

            for i in mia_dict:
                
                exp_barcode = mia_dict[i]["barcode_1"]

                # go through the different checks
                if isReverseComplemented(obs_barcode, exp_barcode) is not None:
                    print(f"found one that's reverse complemented, obs {obs_barcode}, exp {exp_barcode}")
                    info_message = isReverseComplemented(obs_barcode, exp_barcode)
                    break

                if isReverse(obs_barcode, exp_barcode) is not None:
                    print(f"found one that's reversed, obs {obs_barcode}, exp {exp_barcode}")
                    info_message = isReverse(obs_barcode, exp_barcode)
                    break

                if isComplemented(obs_barcode, exp_barcode) is not None:
                    print(f"found one that's complemented, obs {obs_barcode}, exp {exp_barcode}")
                    info_message = isComplemented(obs_barcode, exp_barcode)
                    break

                if isOneMismatch(obs_barcode, exp_barcode) is not None:
                    print(f"found one that's a potential typo, obs {obs_barcode}, exp {exp_barcode}")
                    info_message = isOneMismatch(obs_barcode, exp_barcode)
                    break

            if info_message is not None:
                #print(f"found an explanation: {info_message}, MIA seq i = {i}")
                observed_dict[observed]["explained"] = True
                observed_dict[observed]["sierra"] = True
                observed_dict[observed]["info"] = info_message
                observed_dict[observed]["name"] = mia_dict[i]["name"]
                # remove entry from MIA dictionary
                del mia_dict[i]

            else:
                observed_dict[observed]["info"] = "Couldn't find any matches in the checks that have been run"
                #print(f"Couldn't find any matches in the checks that have run {observed_dict[observed]}")

    print(observed_dict)
    print("=====================")
    print(mia_dict)

    # Write Output to files
    write_sierra(observed_dict)
    write_info(observed_dict,mia_dict)



# Write output for new sierra barcode file
    #order= Sample name \t barcode1 \t barcode2
def write_sierra(observed_dict):
    with open("new_barcodes.txt",mode="w",encoding="UTF8") as out_sierra:
        for key,value in observed_dict.items():
            if value["sierra"] == True:
                out_line = [value["name"],value["barcode_1"]]
                if value["dual"] == True:
                    out_line.append(value["barcode_2"])
                out_line = "\t".join(out_line)

                print(out_line, file = out_sierra)

#Write output for info file
def write_info(observed_dict,mia_dict):    
    with open("barcode_detective_info.txt",mode="w",encoding="UTF8") as out_info:

        #For the barcodes we observe:
        print(f"The following barcodes were observed\n",file = out_info)
        out_header = "\t".join(["barcode","match","info"])
        print(out_header,file = out_info)

        for key,value in observed_dict.items():
            if value["name"] != None:
                out_line = "\t".join([value["full_barcode"],value["name"],value["info"]])
            elif value["explained"] == True:
                out_line = "\t".join([value["full_barcode"],"miscode",value["info"]])
            else:
                out_line = "\t".join([value["full_barcode"],"unexplained",value["info"]])
            
            print(out_line,file = out_info)

        #For barcodes which are still m.i.a

        print(f"\nThe following expected barcodes are still missing\n",file = out_info)
        print("\t".join(["barcode","name"]),file =out_info)
        # check if this breaks if given an empty dictionary
        for key,value in mia_dict.items():
            out_line = "\t".join([value["full_barcode"],value["name"]])
            print(out_line,file = out_info)

def checkPhiX(observed_dict):
    pattern = 'GTATGCCG'
    p = re.compile(pattern)

    for key, value in observed_dict.items():
        observed_seq = value["full_barcode"]
        m = p.search(observed_seq)

        if m is None:
            # Also check rev comp of sequence
            pattern = 'CGGCATAC'
            p = re.compile(pattern)
            m = p.search(observed_seq)

        if m is not None:
            value["explained"] = True
            value["sierra"] = False
            value["info"] = f"found phiX: {observed_seq}"

    return(observed_dict)


def checkPolyG(observed_dict):

    for key, value in observed_dict.items():
        observed_seq = value["full_barcode"]

        if isDual(observed_seq):
            seq_length = min(getDualLengths(observed_seq)) # use minimum length of the dual barcodes
        else:
            seq_length = len(observed_seq)
        pattern = 'G{' + str(seq_length-1) + ',}' # searching for Gs the full length of barcode -1 
        #print(f"pattern: {pattern}") 
        p = re.compile(pattern)
        m = p.search(observed_seq)

        if m is not None:
            value["explained"] = True
            value["sierra"] = False
            value["info"] = f"found polyG: {observed_seq}"

    return(observed_dict)

def check_correct_codes(observed_dict):

    for key, value in observed_dict.items():
        if value["name"] != None:
            value["explained"] = True
            value["sierra"] = True
            value["info"] = f"For sample {value['name']} the provided barcode ({value['full_barcode']}) is correct"

    return(observed_dict)    

def create_output_dicts(all_codes):
    observed_dict= {}
    mia_dict = {}
    for count, line in enumerate(all_codes):
        line_no = "line_"+str(count)

        if line[1] > 0:
            observed_dict[line_no] = create_obs_dict(line)
        else:
            mia_dict[line_no] = create_mia_dict(line)

    return([observed_dict,mia_dict])


def create_mia_dict(line):
    mia_dict = split_if_Dual(line)

    mia_dict.update({"name" : line[2],
                     }
                     )
    return(mia_dict)


def create_obs_dict(line):
    observed_dict = split_if_Dual(line)

    if len(line) == 2:
        line.append(None)

    observed_dict.update({
                    "name" : line[2],
                    "freq":line[1],
                    "explained": False,
                    "sierra": False,
                    "info": ""
                    }
                    )
    
    return(observed_dict)


def split_if_Dual(line):
    observed_dict = {"full_barcode" : line[0]}

    if isDual(line[0]):
        barcodes = line[0].split("_")
        observed_dict.update({"barcode_1" : barcodes[0],
                         "barcode_2" : barcodes[1],
                         "dual":True})
    else:
        observed_dict.update({"barcode_1" : line[0],
                              "barcode_2" : "",
                         "dual":False})
    return(observed_dict)

def isDual(seq):
    if "_" in seq:
        return True
    
def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '_': '_'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def reverseSeq(seq):
    return seq[::-1]

def complementSeq(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '_': '_'}
    complemented_seq = "".join(complement.get(base, base) for base in seq)
    return complemented_seq

def isReverse(observed_seq, expected_seq):
    reversed_seq = reverseSeq(observed_seq)
    if expected_seq == reversed_seq:
        return (f"seq is reversed")

def isComplemented(observed_seq, expected_seq):
    complemented_seq = complementSeq(observed_seq)
    if expected_seq == complemented_seq:
        return (f"seq is complemented")

def isReverseComplemented(observed_seq, expected_seq):
    rev_comp = reverseComplement(observed_seq)
    if rev_comp == expected_seq:
        return (f"seq is reverse complemented") 

# checking for a single typo
def isOneMismatch(observed_seq, expected_seq):
    obs = list(observed_seq)
    exp = list(expected_seq)
    if len(obs) == len(exp):
        diff = False
        for i in range(0, len(obs), 1):
            if (obs[i] != exp[i]):

                # If first mismatch
                if (diff == False):
                    diff = True

                # Second mismatch
                else:
                    diff = False
                    break
        if diff == True:
            return(f"One mismatch between the 2 sequences - potential typo")


# get a list of the 2 barcode lengths
def getDualLengths(seq):
    if not isDual(seq):
        print(f"This isn't a dual barcode - why are you running this function?")
    else:
        seq_lengths = list(map(len, seq.split("_")))
        return seq_lengths


def read_codes(file):
    # define empty list to read in codes to
    all_codes = []
    with open(file,mode="r",encoding="UTF8") as infile:
        for count, line in enumerate(infile):
            #ignore the header line
            if count == 0:
                continue
            line = line.strip()
            line = line.split("\t")

            #convert freq to numeric
            line[1] = float(line[1])
            #print(line)
            all_codes.append(line)

    return(all_codes)
      
def read_options(): #command line prompt
    parser = argparse.ArgumentParser(description="""
    Help for function""")

    parser.add_argument("-found_codes", help="A tsv file containing 3 columns: codes, frequency, name", type=str)

    options = parser.parse_args()
    return options    

if (__name__ == "__main__"):
    main()
