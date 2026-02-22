# MMseqs2 MCP

**Ultra-fast protein sequence search and MSA generation using MMseqs2 via Docker**

An MCP (Model Context Protocol) server for sequence analysis with core tools:
- Search protein sequences against large databases (UniRef100, etc.)
- Generate multiple sequence alignments (MSA) in A3M format
- GPU-accelerated search for high-throughput workflows
- Sequence clustering and database management

## Quick Start with Docker

### Approach 1: Pull Pre-built Image from GitHub

The fastest way to get started. A pre-built Docker image is automatically published to GitHub Container Registry on every release.

```bash
# Pull the latest image
docker pull ghcr.io/macromnex/mmseqs2_mcp:latest

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add mmseqs2 -- docker run -i --rm --user `id -u`:`id -g` -v `pwd`:`pwd` ghcr.io/macromnex/mmseqs2_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker
- Claude Code installed

That's it! The MMseqs2 MCP server is now available in Claude Code.

---

### Approach 2: Build Docker Image Locally

Build the image yourself and install it into Claude Code. Useful for customization or offline environments.

```bash
# Clone the repository
git clone https://github.com/MacromNex/mmseqs2_mcp.git
cd mmseqs2_mcp

# Build the Docker image
docker build -t mmseqs2_mcp:latest .

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add mmseqs2 -- docker run -i --rm --user `id -u`:`id -g` -v `pwd`:`pwd` mmseqs2_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker
- Claude Code installed
- Git (to clone the repository)

**About the Docker Flags:**
- `-i` — Interactive mode for Claude Code
- `--rm` — Automatically remove container after exit
- `` --user `id -u`:`id -g` `` — Runs the container as your current user, so output files are owned by you (not root)
- `-v` — Mounts your project directory so the container can access your data

---

## Verify Installation

After adding the MCP server, you can verify it's working:

```bash
# List registered MCP servers
claude mcp list

# You should see 'mmseqs2' in the output
```

---

## Next Steps

- **Detailed documentation**: See [detail.md](detail.md) for comprehensive guides on:
  - Available MCP tools and parameters
  - Local Python environment setup (alternative to Docker)
  - Database setup (UniRef100, etc.)
  - Example workflows and use cases
  - Troubleshooting

---

## Usage Examples

Once registered, you can use the MMseqs2 tools directly in Claude Code. Here are some common workflows:

### Example 1: Generate MSA for a Protein Sequence

```
Can you create an MSA for the DHFR protein sequence using the mmseqs2 MCP? The sequence is: MISLIAALAVDRVIGMENAMPWNLPADLAWFKRNTLNKPVIMGRHTWESIGRPLPGRKNIILSSQPGTDDRVTWVKSVDEAIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFEILERR. Save the output A3M file to /path/to/output/DHFR.a3m.
```

### Example 2: Sequence Similarity Search

```
I have a FASTA file at /path/to/query.fasta. Can you search it against the UniRef100 database using MMseqs2 and save the alignment results to /path/to/results/? Use sensitivity 7.5 for high-quality hits.
```

### Example 3: Batch MSA Generation

```
I have multiple protein sequences in /path/to/sequences.fasta that I need MSAs for. Can you run MMseqs2 search for each sequence and generate A3M files in /path/to/msas/ directory?
```

---

## Troubleshooting

**Docker not found?**
```bash
docker --version  # Install Docker if missing
```

**Claude Code not found?**
```bash
# Install Claude Code
npm install -g @anthropic-ai/claude-code
```

**Database not found?**
- Set the `MMSEQS2_DB_PATH` environment variable to point to your database
- Add `-e MMSEQS2_DB_PATH=/path/to/db` to the docker run command

---

## License

GPL — Based on [MMseqs2](https://github.com/soedinglab/MMseqs2) by Söding Lab.
