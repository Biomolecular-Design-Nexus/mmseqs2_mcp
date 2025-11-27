"""
Model Context Protocol (MCP) for mmseqs2

MMseqs2 is an ultra-fast and sensitive sequence search and clustering suite for protein and nucleotide sequences.
It provides tools for large-scale sequence analysis including homology search, clustering, and taxonomic assignment.
MMseqs2 can run 10000 times faster than BLAST while maintaining high sensitivity for biological sequence analysis.

This MCP Server provides tools for generating Multiple Sequence Alignments (MSA) using MMseqs2.
"""

import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Optional, Literal
from fastmcp import FastMCP

# Server definition
mcp = FastMCP(name="mmseqs2")

# Database path - load from environment variable or use default
MMSEQS2_DB_PATH = os.path.expanduser(
    os.environ.get(
        "MMSEQS2_DB_PATH",
        "~/.db/protein/uniref100/uniref100.fasta.db_padded"
    )
)

# Find mmseqs2 binary
# Try to find mmseqs in the environment or system PATH
MMSEQS_BIN = None
for possible_path in [
    os.path.join(os.path.dirname(os.path.dirname(__file__)), "env", "bin", "mmseqs"),
    shutil.which("mmseqs")
]:
    if possible_path and os.path.exists(possible_path):
        MMSEQS_BIN = possible_path
        break

if MMSEQS_BIN is None:
    MMSEQS_BIN = "mmseqs"  # Fallback to PATH lookup


def _generate_msa_impl(
    sequence: Optional[str] = None,
    fasta_file: Optional[str] = None,
    sequence_name: str = "query",
    output_dir: Optional[str] = None,
    database_path: Optional[str] = None,
    gpu: bool = True,
    threads: int = 64,
    sensitivity: float = 7.5,
    num_iterations: int = 10,
    e_value: float = 0.001,
    max_seqs: int = 100000,
    return_format: Literal["a3m", "path"] = "a3m"
) -> str:
    """
    Generate Multiple Sequence Alignment (MSA) for a protein sequence using MMseqs2.

    This tool runs the complete MMseqs2 pipeline to search against a protein database
    and generate a multiple sequence alignment in a3m format.

    Args:
        sequence: Protein sequence as a string (one-letter amino acid codes).
                 Either sequence or fasta_file must be provided.
        fasta_file: Path to a FASTA file containing the query sequence(s).
                   Either sequence or fasta_file must be provided.
        sequence_name: Name/identifier for the sequence (used if sequence is provided).
        output_dir: Directory to store output files. If None, uses a temporary directory.
        database_path: Path to the MMseqs2 database. If None, uses MMSEQS2_DB_PATH
                      environment variable or defaults to UniRef100 padded database.
        gpu: Use GPU acceleration for search (default: True).
        threads: Number of CPU threads to use (default: 64).
        sensitivity: Sensitivity parameter for search (default: 7.5, higher = more sensitive).
        num_iterations: Number of search iterations (default: 10).
        e_value: E-value threshold (default: 0.001).
        max_seqs: Maximum number of sequences to return (default: 100000).
        return_format: Return format - "a3m" returns the MSA content as string,
                      "path" returns the path to the output file.

    Returns:
        Either the MSA content as a string (a3m format) or the path to the output file,
        depending on return_format parameter.

    Example:
        >>> msa = generate_msa(
        ...     sequence="MKTFIFLALLGAAVAFPVDDDDKIVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCASAGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAGYLEGGKDSCQGDSGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN",
        ...     sequence_name="DHFR"
        ... )
    """
    # Use global database path if not specified
    if database_path is None:
        database_path = MMSEQS2_DB_PATH

    # Validate inputs
    if sequence is None and fasta_file is None:
        raise ValueError("Either 'sequence' or 'fasta_file' must be provided")

    if sequence is not None and fasta_file is not None:
        raise ValueError("Only one of 'sequence' or 'fasta_file' should be provided")

    # Check if database exists
    if not os.path.exists(database_path):
        raise FileNotFoundError(f"Database not found at: {database_path}")

    # Determine if we need to use a temporary directory
    use_temp_dir = output_dir is None
    if use_temp_dir:
        temp_dir = tempfile.mkdtemp(prefix="mmseqs2_")
        work_dir = temp_dir
    else:
        work_dir = output_dir
        os.makedirs(work_dir, exist_ok=True)

    try:
        # Create query FASTA file if sequence is provided
        if sequence:
            query_fasta = os.path.join(work_dir, f"{sequence_name}.fasta")
            with open(query_fasta, 'w') as f:
                f.write(f">{sequence_name}\n{sequence}\n")
        else:
            query_fasta = fasta_file
            # Extract sequence name from fasta file for output naming
            with open(fasta_file, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('>'):
                    sequence_name = first_line[1:].split()[0]

        # Define file paths
        query_db = os.path.join(work_dir, f"{sequence_name}_db")
        result_db = os.path.join(work_dir, f"{sequence_name}_result_db")
        msa_db = os.path.join(work_dir, f"{sequence_name}_msa_db")
        tmp_dir = os.path.join(work_dir, "tmp")
        msa_dir = os.path.join(work_dir, f"{sequence_name}_msa")
        output_a3m = os.path.join(work_dir, f"{sequence_name}.a3m")

        # Step 1: Create query database
        print(f"Step 1/4: Creating query database from {query_fasta}...")
        print(f"  Using mmseqs binary: {MMSEQS_BIN}")
        result = subprocess.run([
            MMSEQS_BIN, "createdb", query_fasta, query_db
        ], check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)

        # Step 2: Search against database
        print(f"Step 2/4: Searching against database (this may take a while)...")
        search_cmd = [
            MMSEQS_BIN, "search",
            query_db, database_path, result_db, tmp_dir,
            "--threads", str(threads),
            "-s", str(sensitivity),
            "--num-iterations", str(num_iterations),
            "-e", str(e_value),
            "--max-seqs", str(max_seqs)
        ]

        if gpu:
            search_cmd.insert(3, "--gpu")
            search_cmd.insert(4, "1")

        print(f"  Command: {' '.join(search_cmd)}")
        result = subprocess.run(search_cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)

        # Step 3: Convert result to MSA
        print(f"Step 3/4: Converting search results to MSA format...")
        result = subprocess.run([
            MMSEQS_BIN, "result2msa",
            query_db, database_path, result_db, msa_db,
            "--msa-format-mode", "6"
        ], check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)

        # Step 4: Unpack MSA database
        print(f"Step 4/4: Unpacking MSA database...")
        result = subprocess.run([
            MMSEQS_BIN, "unpackdb",
            msa_db, msa_dir,
            "--unpack-suffix", ".a3m"
        ], check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)

        # Step 5: Concatenate a3m files
        print(f"Step 5/5: Concatenating a3m files...")
        a3m_files = sorted(Path(msa_dir).glob("*.a3m"))
        print(f"  Found {len(a3m_files)} a3m file(s)")
        with open(output_a3m, 'w') as outfile:
            for a3m_file in a3m_files:
                with open(a3m_file, 'r') as infile:
                    outfile.write(infile.read())

        print(f"\nMSA generation complete! Output saved to: {output_a3m}")

        # Read and return the result
        if return_format == "a3m":
            with open(output_a3m, 'r') as f:
                result = f.read()
            print(f"Returning MSA content ({len(result)} characters)")
        else:  # return_format == "path"
            result = output_a3m
            print(f"Returning path: {result}")

        return result

    except subprocess.CalledProcessError as e:
        error_msg = f"MMseqs2 command failed: {e.cmd}\n"
        if e.stderr:
            error_msg += f"Error output: {e.stderr.decode()}\n"
        raise RuntimeError(error_msg)

    finally:
        # Clean up temporary directory if used
        if use_temp_dir:
            shutil.rmtree(temp_dir, ignore_errors=True)


@mcp.tool()
def generate_msa(
    sequence: Optional[str] = None,
    fasta_file: Optional[str] = None,
    sequence_name: str = "query",
    output_dir: Optional[str] = None,
    database_path: Optional[str] = None,
    gpu: bool = True,
    threads: int = 64,
    sensitivity: float = 7.5,
    num_iterations: int = 10,
    e_value: float = 0.001,
    max_seqs: int = 100000,
    return_format: Literal["a3m", "path"] = "a3m"
) -> str:
    """MCP tool wrapper for generate_msa"""
    return _generate_msa_impl(
        sequence=sequence,
        fasta_file=fasta_file,
        sequence_name=sequence_name,
        output_dir=output_dir,
        database_path=database_path,
        gpu=gpu,
        threads=threads,
        sensitivity=sensitivity,
        num_iterations=num_iterations,
        e_value=e_value,
        max_seqs=max_seqs,
        return_format=return_format
    )


@mcp.tool()
def generate_msa_from_file(
    fasta_file: str,
    output_dir: str,
    database_path: Optional[str] = None,
    gpu: bool = True,
    threads: int = 64,
    sensitivity: float = 7.5,
    num_iterations: int = 10,
    e_value: float = 0.001,
    max_seqs: int = 100000
) -> str:
    """
    Generate Multiple Sequence Alignment (MSA) from a FASTA file and save to output directory.

    This is a convenience wrapper around generate_msa that always saves results to disk.
    Useful when you want to preserve all intermediate files and the final MSA.

    Args:
        fasta_file: Path to a FASTA file containing the query sequence(s).
        output_dir: Directory to store output files (will be created if it doesn't exist).
        database_path: Path to the MMseqs2 database. If None, uses MMSEQS2_DB_PATH
                      environment variable or defaults to UniRef100 padded database.
        gpu: Use GPU acceleration for search (default: True).
        threads: Number of CPU threads to use (default: 64).
        sensitivity: Sensitivity parameter for search (default: 7.5).
        num_iterations: Number of search iterations (default: 10).
        e_value: E-value threshold (default: 0.001).
        max_seqs: Maximum number of sequences to return (default: 100000).

    Returns:
        Path to the generated .a3m file containing the MSA.

    Example:
        >>> output_path = generate_msa_from_file(
        ...     fasta_file="examples/DHFR.fasta",
        ...     output_dir="examples/output"
        ... )
    """
    return _generate_msa_impl(
        fasta_file=fasta_file,
        output_dir=output_dir,
        database_path=database_path,
        gpu=gpu,
        threads=threads,
        sensitivity=sensitivity,
        num_iterations=num_iterations,
        e_value=e_value,
        max_seqs=max_seqs,
        return_format="path"
    )


if __name__ == "__main__":
    mcp.run()
