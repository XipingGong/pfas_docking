#!/usr/bin/env python3

import argparse
import os
import math

def split_file(input_file, output_name, num_parts):
    # Read all lines
    with open(input_file, "r") as f:
        lines = f.readlines()

    total_lines = len(lines)
    lines_per_file = math.ceil(total_lines / num_parts)

    print(f"Total lines: {total_lines}")
    print(f"Splitting into {num_parts} parts, about {lines_per_file} lines each.")

    for i in range(num_parts):
        start = i * lines_per_file
        end = min(start + lines_per_file, total_lines)
        output_file = f"{output_name}{i+1}.log"

        with open(output_file, "w") as f_out:
            f_out.writelines(lines[start:end])

        print(f"Created {output_file} with {end - start} lines.")

    print(f"\n✅ Successfully split '{input_file}' into {num_parts} smaller files.")

def main():
    parser = argparse.ArgumentParser(description="Split a large file into multiple smaller files.")
    parser.add_argument("input_file", help="Path to the input file to split.")
    parser.add_argument("--output_name", default="output", help="Base name for output files (default: output)")
    parser.add_argument("--num_parts", type=int, default=50, help="Number of parts to split into (default: 50)")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"❌ Error: Input file '{args.input_file}' does not exist.")
        return

    split_file(args.input_file, args.output_name, args.num_parts)

if __name__ == "__main__":
    main()

