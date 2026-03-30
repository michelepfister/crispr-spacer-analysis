import argparse
import sys

from pathlib import Path
from process import process

def main():
    """Main CLI function for parsing FASTA spacer files."""
    parser = argparse.ArgumentParser(
        description="Parse FASTA files containing CRISPR spacer sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py input.fasta
  python main.py /path/to/spacers.fasta --output results
  python main.py input.fasta
        """
    )
    
    parser.add_argument(
        "file_path",
        help="Path to the FASTA file to parse"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="Output file path (default: print to stdout)",
        default=None
    )    
    args = parser.parse_args()
    
    # Check if input file exists
    input_path = Path(args.file_path)
    if not input_path.exists():
        print(f"Error: File '{args.file_path}' not found.", file=sys.stderr)
        sys.exit(1)
    
    if not input_path.is_file():
        print(f"Error: '{args.file_path}' is not a file.", file=sys.stderr)
        sys.exit(1)
    
    try:
        process(args)

    
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Parse error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

