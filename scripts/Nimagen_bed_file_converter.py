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
            primer_right_start = int(line.split("\t")[2]) + 1
            primer_right_end = int(line.split("\t")[2]) + 1 + len(line.split("\t")[4].rstrip("\n"))
            left_primer_seq = line.split("\t")[3]
            right_primer_seq = line.split("\t")[4].rstrip("\n")

            outfile.write(f"MN908947.3\t{primer_left_start}\t{primer_left_end}\t{primer_name}LEFT_0\t{pooler(primer_name)}\t+\t{left_primer_seq}\n")
            outfile.write(f"MN908947.3\t{primer_right_start}\t{primer_right_end}\t{primer_name}RIGHT_0\t{pooler(primer_name)}\t-\t{right_primer_seq}\n")