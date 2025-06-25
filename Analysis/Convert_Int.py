import sys

def convert_bed_scores(input_bed, output_bed):
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            if line.strip() == "" or line.startswith("#"):
                continue  # skip blank lines or comments
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                try:
                    fields[4] = str(int(float(fields[4])))
                except ValueError:
                    pass  # leave the original value if it can't be converted
            outfile.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_score_to_int.py input.bed output.bed")
        sys.exit(1)

    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    convert_bed_scores(input_bed, output_bed)
