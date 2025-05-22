#!/usr/bin/env python3

import argparse
import os

def split_file(input_file, output_name, target_parts):
    with open(input_file, "r") as f:
        lines = f.readlines()

    total_lines = len(lines)
    if total_lines == 0:
        print("❌ Error: Input file is empty.")
        return

    # Calculate how many files get one extra line (to balance the split)
    lines_per_file = total_lines // target_parts
    extra_files = total_lines % target_parts

    print(f"📄 Total lines: {total_lines}")
    print(f"🎯 Target parts: {target_parts}")
    print(f"📦 Actual parts: {target_parts}")
    print(f"📏 {extra_files} files will have {lines_per_file + 1} lines, and {target_parts - extra_files} will have {lines_per_file} lines\n")

    start = 0
    for i in range(target_parts):
        current_lines = lines_per_file + 1 if i < extra_files else lines_per_file
        end = start + current_lines
        output_file = f"{output_name}{i+1}.log"

        with open(output_file, "w") as f_out:
            f_out.writelines(lines[start:end])
        print(f"✅ Created {output_file} with {current_lines} lines.")

        start = end

    print(f"\n🎉 Successfully split '{input_file}' into {target_parts} balanced parts.")

def main():
    parser = argparse.ArgumentParser(description="Split a file into equal-sized parts with exact line counts.")
    parser.add_argument("input_file", help="Path to the input file to split.")
    parser.add_argument("--output_name", default="output", help="Base name for output files (default: output)")
    parser.add_argument("--num_parts", type=int, default=50, help="Exact number of output files (default: 50)")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"❌ Error: File '{args.input_file}' does not exist.")
        return

    if args.num_parts <= 0:
        print(f"❌ Error: num_parts must be greater than 0.")
        return

    if args.num_parts > len(open(args.input_file).readlines()):
        print(f"❌ Error: num_parts exceeds number of lines in the file.")
        return

    split_file(args.input_file, args.output_name, args.num_parts)

if __name__ == "__main__":
    main()

