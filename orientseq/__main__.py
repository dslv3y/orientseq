import argparse
import pysam
from Bio.Seq import Seq

def parse_arguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', help='input file to orient', required=True)
  parser.add_argument('-o', '--output', help='output file to write', required=True)
  parser.add_argument('-f', '--format', choices=['fq', 'fa', 'sam', 'bam'], help='file format, by default reads the input file extension')
  parser.add_argument('-t', '--threshold', type=int, default=5, help='maximum possible difference between polyA and polyT to count as ambigous')
  parser.add_argument('-a', '--ambigous', help='output file for ambigous reads, not printed anywhere by default')
  return parser.parse_args()

def detect_format(file_name):
  fmt = file_name.split(".")[-1]
  if fmt in ['fq', 'fastq', 'fa', 'fasta']:
    return 0, 'fx'
  elif fmt in ['sam', 'bam']:
    return 0, fmt
  else:
    return 1, fmt

def find_poly(sequence, max_mismatches):
  max_poly = {'T': 0, 'A': 0}
  current_poly_length = current_mismatches = 0
  current_poly_base = 'T'
  for i,base in enumerate(sequence):
    if i == len(sequence) // 2:
      current_poly_length = current_mismatches = 0
      current_poly_base = 'A'
    if base == current_poly_base:
      current_poly_length += 1
      current_mismatches = 0
      max_poly[current_poly_base] = max(max_poly[current_poly_base], current_poly_length)
    else:
      if current_mismatches < max_mismatches:
        current_mismatches += 1
        current_poly_length += 1
      else:
        max_poly[current_poly_base] = max(max_poly[current_poly_base], current_poly_length)
        current_poly_length = current_mismatches = 0
  return max_poly

def read_and_write_alignment(f, out, threshold):
  stats = {'fwd_cnt': 0, 'fwd_len': 0, 'revcomp_cnt': 0, 'revcomp_len': 0, 'ambigous_cnt': 0, 'ambigous_len': 0, 'polyA_len': 0, 'polyT_len': 0}
  for read in f:
    max_poly = find_poly(read.query_sequence, 0)
    diff = max_poly['A'] - max_poly['T']
    if abs(diff) < threshold:
      stats['ambigous_cnt'] += 1
      stats['ambigous_len'] += len(read.query_sequence)
    elif diff < 0:
      read.query_sequence = str(Seq(read.query_sequence).reverse_complement())
      stats['revcomp_cnt'] += 1
      stats['revcomp_len'] += len(read.query_sequence)
      stats['polyT_len'] += max_poly['T']
    else:
      stats['fwd_cnt'] += 1
      stats['fwd_len'] += len(read.query_sequence)
      stats['polyA_len'] += max_poly['A']
    out.write(read)
  return stats

def read_and_write_fastx(f, out, threshold):
  stats = {'fwd_cnt': 0, 'fwd_len': 0, 'revcomp_cnt': 0, 'revcomp_len': 0, 'ambigous_cnt': 0, 'ambigous_len': 0, 'polyA_len': 0, 'polyT_len': 0}
  for read in f:
    max_poly = find_poly(read.sequence, 0)
    diff = max_poly['A'] - max_poly['T']
    if abs(diff) < threshold:
      stats['ambigous_cnt'] += 1
      stats['ambigous_len'] += len(read.sequence)
    elif diff < 0:
      read.sequence = str(Seq(read.sequence).reverse_complement())
      stats['revcomp_cnt'] += 1
      stats['revcomp_len'] += len(read.sequence)
      stats['polyT_len'] += max_poly['T']
    else:
      stats['fwd_cnt'] += 1
      stats['fwd_len'] += len(read.sequence)
      stats['polyA_len'] += max_poly['A']
    out.write(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n")
  return stats

def print_stats(stats, input_file_name):
  print(f"""Stats for orienting {input_file_name}:
  ========================================================

  Normal reads count:           {stats['fwd_cnt']}
  Normal reads average length:  {stats['fwd_len']/(stats['fwd_cnt'] + 0.001)}
  Average polyA length:         {stats['polyA_len']/(stats['fwd_cnt'] + 0.001)}
  ========================================================

  Revcomp reads count:           {stats['revcomp_cnt']}
  Revcomp reads average length:  {stats['revcomp_len']/(stats['revcomp_cnt'] + 0.001)}
  Average polyT length:          {stats['polyT_len']/(stats['revcomp_cnt'] + 0.001)}
  ========================================================

  Ambigous reads count:          {stats['ambigous_cnt']}
  Ambigous reads average length: {stats['ambigous_len']/(stats['ambigous_cnt'] + 0.001)}
  ========================================================
  """)

def main():
  args = parse_arguments()
  fmt_status, fmt = detect_format(args.input)
  fmt = args.format if args.format else fmt
  if fmt_status == 1:
    raise ValueError(f"Can't parse .{fmt} files!")
  Path(args.output).touch()
  if fmt == 'fx':
    with pysam.FastxFile(args.input, "r") as f, open(args.output, "w") as out:
      stats = read_and_write_fastx(f, out, args.threshold)
      print_stats(stats, args.input)
  elif fmt == 'sam':
    with pysam.AlignmentFile(args.input, "r", check_sq=False) as f, pysam.AlignmentFile(args.output, "w", header=f.header) as out:
      stats = read_and_write_alignment(f, out, args.threshold)
      print_stats(stats, args.input)
  elif fmt == 'bam':
    with pysam.AlignmentFile(args.input, "rb", check_sq=False) as f, pysam.AlignmentFile(args.output, "wb", header=f.header) as out:
      stats = read_and_write_alignment(f, out, args.threshold)
      print_stats(stats, args.input)

if __name__ == '__main__':
  main()