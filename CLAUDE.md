# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MMseqs2 MCP is a Model Context Protocol server that wraps MMseqs2 (ultra-fast protein sequence search) to provide MSA (Multiple Sequence Alignment) generation as MCP tools. Built with the FastMCP framework, it integrates with Claude Code and Gemini CLI.

## Key Commands

```bash
# Setup environment (creates conda env, installs mmseqs2 GPU binary + fastmcp)
bash quick_setup.sh

# Run the MCP server directly
./env/bin/python src/server.py

# Debug with MCP inspector
fastmcp run src/server.py:mcp --transport http --port 8001
npx @modelcontextprotocol/inspector

# Docker build and run (mount database volume)
docker build -t mmseqs2-mcp .
docker run --gpus all -v /path/to/uniref100:/data/uniref100 mmseqs2-mcp
```

## Architecture

**Single-file server** (`src/server.py`): The entire MCP server lives in one file using FastMCP. It exposes two tools:
- `generate_msa()` — accepts a raw sequence string or FASTA file path, returns MSA content or file path
- `generate_msa_from_file()` — convenience wrapper that always saves to disk

Both tools call `_generate_msa_impl()` which orchestrates a 5-step MMseqs2 subprocess pipeline: `createdb → search → result2msa → unpackdb → concatenate`. Intermediate files are auto-cleaned after completion.

**External dependency**: Requires a pre-built MMseqs2 padded database (UniRef100). Path is configured via `MMSEQS2_DB_PATH` env var (default: `~/.db/protein/uniref100/uniref100.fasta.db_padded`).

**Binary resolution**: The server looks for `mmseqs` binary first in `../env/bin/mmseqs` (relative to src/), then falls back to system PATH.

## Key Environment Variables

- `MMSEQS2_DB_PATH` — path to the MMseqs2 padded reference database (required at runtime)
- `CUDA_VISIBLE_DEVICES` — set automatically per-query when `cuda_device` parameter is specified

## Dependencies

Only Python dependency is `fastmcp` (see `requirements.txt`). MMseqs2 is a standalone binary downloaded separately (GPU version from `https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz`).
