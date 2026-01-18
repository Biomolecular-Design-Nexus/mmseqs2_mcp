#!/bin/bash
#===============================================================================
# MMseqs2 MCP Quick Setup Script
#===============================================================================
# This script sets up the complete environment for MMseqs2 MCP server.
# Sequence search, clustering, and MSA generation using MMseqs2.
#
# After cloning the repository, run this script to set everything up:
#   cd mmseqs2_mcp
#   bash quick_setup.sh
#
# Once setup is complete, register in Claude Code with the config shown at the end.
#
# Options:
#   --skip-env        Skip conda environment creation
#   --skip-mmseqs     Skip MMseqs2 binary download
#   --help            Show this help message
#===============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/env"
PYTHON_VERSION="3.10"

# Print banner
echo -e "${BLUE}"
echo "=============================================="
echo "      MMseqs2 MCP Quick Setup Script         "
echo "=============================================="
echo -e "${NC}"

# Helper functions
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check for conda/mamba
check_conda() {
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
        info "Using mamba (faster package resolution)"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
        info "Using conda"
    else
        error "Neither conda nor mamba found. Please install Miniconda or Mambaforge first."
        exit 1
    fi
}

# Parse arguments
SKIP_ENV=false
SKIP_MMSEQS=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-env) SKIP_ENV=true; shift ;;
        --skip-mmseqs) SKIP_MMSEQS=true; shift ;;
        -h|--help)
            echo "Usage: ./quick_setup.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-env        Skip conda environment creation"
            echo "  --skip-mmseqs     Skip MMseqs2 binary download"
            echo "  -h, --help        Show this help message"
            exit 0
            ;;
        *) warn "Unknown option: $1"; shift ;;
    esac
done

# Check prerequisites
info "Checking prerequisites..."
check_conda
success "Prerequisites check passed"

# Step 1: Create conda environment
echo ""
echo -e "${BLUE}Step 1: Setting up conda environment${NC}"

if [ "$SKIP_ENV" = true ]; then
    info "Skipping environment creation (--skip-env)"
elif [ -d "$ENV_DIR" ] && [ -f "$ENV_DIR/bin/python" ]; then
    info "Environment already exists at: $ENV_DIR"
else
    info "Creating conda environment with Python ${PYTHON_VERSION}..."
    $CONDA_CMD create -p "$ENV_DIR" python=${PYTHON_VERSION} -y

    info "Installing fastmcp..."
    "${ENV_DIR}/bin/pip" install fastmcp
fi

# Step 2: Download MMseqs2
echo ""
echo -e "${BLUE}Step 2: Installing MMseqs2${NC}"

if [ "$SKIP_MMSEQS" = true ]; then
    info "Skipping MMseqs2 download (--skip-mmseqs)"
elif [ -f "${ENV_DIR}/bin/mmseqs" ]; then
    info "MMseqs2 already installed"
else
    if command -v wget &> /dev/null; then
        info "Downloading MMseqs2 GPU version..."
        wget -q https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz
        tar xzf mmseqs-linux-gpu.tar.gz
        cp mmseqs/bin/mmseqs "${ENV_DIR}/bin/"
        rm -rf mmseqs mmseqs-linux-gpu.tar.gz
        success "MMseqs2 installed"
    else
        warn "wget not found, skipping mmseqs2 GPU installation. Install manually or install wget."
    fi
fi

# Step 3: Verify installation
echo ""
echo -e "${BLUE}Step 3: Verifying installation${NC}"

"${ENV_DIR}/bin/python" -c "import fastmcp; print('Core packages OK')" && success "Core packages verified" || error "Package verification failed"

if [ -f "${ENV_DIR}/bin/mmseqs" ]; then
    "${ENV_DIR}/bin/mmseqs" version
    success "MMseqs2 installed"
else
    warn "MMseqs2 binary not found"
fi

# Print summary
echo ""
echo -e "${GREEN}=============================================="
echo "           Setup Complete!"
echo "==============================================${NC}"
echo ""
echo "Environment: $ENV_DIR"
echo ""
echo -e "${YELLOW}Note:${NC} Set MMSEQS2_DB_PATH environment variable to point to your database."
echo ""
echo -e "${YELLOW}Claude Code Configuration:${NC}"
echo ""
cat << EOF
{
  "mcpServers": {
    "mmseqs2": {
      "command": "${ENV_DIR}/bin/python",
      "args": ["${SCRIPT_DIR}/src/mmseqs2_mcp.py"],
      "env": {
        "MMSEQS2_DB_PATH": "/path/to/your/database"
      }
    }
  }
}
EOF
echo ""
echo "To add to Claude Code:"
echo "  claude mcp add mmseqs2 -- ${ENV_DIR}/bin/python ${SCRIPT_DIR}/src/mmseqs2_mcp.py"
echo ""
