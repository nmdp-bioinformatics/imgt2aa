#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_URL="https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip"
VERSION="3.65.0"
EXON3=0
DOWNLOAD=1
URL="$DEFAULT_URL"

usage() {
  cat <<EOF
Usage: $(basename "$0") [-v VERSION] [-e] [--skip-download] [--url URL]

Download hla.xml, build the .db files, and run quarterly output generation.

Options:
  -v VERSION        HLA database version to write under data/ (default: $VERSION)
  -e                Enable exon3 mode for both Perl steps
  --skip-download   Reuse existing hla.xml instead of downloading it
  --url URL         Override the IMGT XML zip URL
  -h, --help        Show this help
EOF
}

while (($# > 0)); do
  case "$1" in
    -v)
      VERSION="${2:?missing value for -v}"
      shift 2
      ;;
    -e)
      EXON3=1
      shift
      ;;
    --skip-download)
      DOWNLOAD=0
      shift
      ;;
    --url)
      URL="${2:?missing value for --url}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

cd "$ROOT_DIR"

if (( DOWNLOAD )); then
  echo "Downloading hla.xml.zip from $URL"
  curl -L "$URL" -o hla.xml.zip
  echo "Extracting hla.xml"
  unzip -o hla.xml.zip hla.xml
fi

IMGTT_ARGS=(-v "$VERSION")
QUARTERLY_ARGS=(-v "$VERSION")

if (( EXON3 )); then
  IMGTT_ARGS=(-e "${IMGTT_ARGS[@]}")
  QUARTERLY_ARGS=(-e "${QUARTERLY_ARGS[@]}")
fi

echo "Running imgt2aa.pl ${IMGTT_ARGS[*]}"
perl "$ROOT_DIR/imgt2aa.pl" "${IMGTT_ARGS[@]}"

echo "Running quarterly.pl ${QUARTERLY_ARGS[*]}"
perl "$ROOT_DIR/quarterly.pl" "${QUARTERLY_ARGS[@]}"

echo "Done"
