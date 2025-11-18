# Nimagen provide the bed file in a format that is not directly translatable to use in tools like iVar or the ARTIC pipeline.
# This bit of code will convert that file into one usable by those.
# What we have been given is the coordinates of the insert, so need to add or subtract the len of the primer as needed.
# New bed file has been written to current ARTIC specs

def pooler(primer_name):
    number = int(primer_name.split("_")[1])
    if number % 2 == 0:
        return "2"
    else:
        return "1"

with open("resources/Nimagen_primer_details/SARS_CoV_2_V5_insert.bed", "r") as input_primers:
    with open ("resources/NimagenV5_covid19.bed","w") as outfile:
        for line in input_primers:
            if line.split("\t")[0] == "Name":
                continue
            primer_name = line.split("\t")[0].rstrip("pAB")
            primer_left_start = int(line.split("\t")[1]) - 1 - len(line.split("\t")[3])
            primer_left_end = int(line.split("\t")[1]) - 1
            primer_right_start = int(line.split("\t")[2])
            primer_right_end = int(line.split("\t")[2]) + len(line.split("\t")[4].rstrip("\n"))
            left_primer_seq = line.split("\t")[3]
            right_primer_seq = line.split("\t")[4].rstrip("\n")

            outfile.write(f"MN908947.3\t{primer_left_start}\t{primer_left_end}\t{primer_name}LEFT_0\t{pooler(primer_name)}\t+\t{left_primer_seq}\n")
            outfile.write(f"MN908947.3\t{primer_right_start}\t{primer_right_end}\t{primer_name}RIGHT_0\t{pooler(primer_name)}\t-\t{right_primer_seq}\n")

# First iteration of this tool highlighted that some of the positions given by Nimagen don't exactly match the sequence of the primer when checked against MN908947
# This section is going to check the sequece of the regions in the bedfile against primer seq

with open("resources/NimagenV5_covid19_V2.bed","r") as bed_file:
    position_dict = {}
    primer_dict = {}
    for line in bed_file:
        pos_list = []
        pos_list.append(int(line.split()[1]))
        pos_list.append(int(line.split()[2]))
        position_dict[line.split()[3]] = pos_list
        primer_dict[line.split()[3]] = line.split()[6]
    
    with open("resources/MN908947.flat.fa", "r") as ref:
        seq_dict = {}
        for line in ref:
            if line == ">MN908947.3":
                continue
            else:
                sequence = line
            for k,v in position_dict.items():
                if k.split("_")[2] == "RIGHT":
                    ref_seq = sequence[v[0]:v[1]]
                    rev_seq = ref_seq[::-1]
                    rev_comp = ""
                    for i in rev_seq:
                        if i == "A":
                            j = "T"
                        if i == "G":
                            j = "C"
                        if i == "C":
                            j = "G"
                        if i == "T":
                            j = "A"
                        rev_comp += j
                    seq_dict[k] = rev_comp
                else:
                    seq_dict[k] = sequence[v[0]:v[1]]
        
        print(f"Primer\tprimer_seq\tref_seq\tmatch")
        for k,v in primer_dict.items():
            print(f"{k}\t{v}\t{seq_dict[k]}\t{v == seq_dict[k]}")           
