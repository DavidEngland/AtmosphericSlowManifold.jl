#!/usr/bin/env bash
# scripts/audit_seb.sh — scan a SpectralBL-Analytics/data tree for surface energy
# budget (SEB) fields: ground/soil heat flux (G), soil/snow/ice thermistors,
# skin temperature, net radiation, soil moisture, snow depth.
#
# Usage:
#   chmod +x audit_seb.sh
#   ./audit_seb.sh ~/Documents/GitHub/SpectralBL-Analytics/data
#
# Optionally install `ncdump` (part of netcdf-bin / libnetcdf-dev) to also
# scan .nc variable tables — most raw campaign data here is NetCDF.

set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
DEFAULT_ROOT="$SCRIPT_DIR/../../SpectralBL-Analytics/data"
ROOT="${1:-$DEFAULT_ROOT}"

if [[ ! -d "$ROOT" ]]; then
  echo "error: data root is not a directory: $ROOT" >&2
  echo "usage: $0 [DATA_ROOT]" >&2
  exit 2
fi

# --- filenames that hint at SEB data ---------------------------------------
FNAME_PAT='soil|snow|ice_phase|thermistor|flux.?plate|skin|ground.?heat|SEB|surface_flux|aws'

# --- column/header tokens (AmeriFlux/FLUXNET-style + generic) --------------
#   G, G_1_1_1        -> ground heat flux
#   TS, TS_1_1_1      -> soil temperature
#   SWC, SWC_1_1_1    -> soil water content
#   SHF                -> soil heat flux (alt naming)
#   Tskin, Tsfc         -> skin/surface temperature
#   Rn, NETRAD          -> net radiation
#   LWd/LWu, SWd/SWu  -> radiation components
COL_PAT='^(g(_[0-9]+)*|ts(_?[0-9]+)*|swc(_?[0-9]+)*|shf|t_?skin|t_?sfc|tice|tsn(ow|w)|rn|rnet|netrad|lwd|lwu|swd|swu|qsub|soil_?temp(erature)?|soil_?moist(ure)?|snow_?depth|ground_?heat(_?flux)?)$'

print_header_matches() {
  awk -v pattern="$COL_PAT" '
    NF {
      gsub(/\r/, "")
      count = split($0, fields, /[,;[:space:]]+/)
      for (i = 1; i <= count; i++) {
        field = fields[i]
        gsub(/^"|"$/, "", field)
        normalized = tolower(field)
        if (normalized ~ pattern ||
            normalized ~ /(soil|snow|thermistor|ground.*heat|skin.*temp|net.*rad|sensible_heat_flux)/) {
          print "    " field
        }
      }
      exit
    }
  ' "$1"
}

echo "############################################"
echo "# 1) Filenames matching SEB-related keywords"
echo "############################################"
find "$ROOT" -type f \
  \( -iname "*.csv" -o -iname "*.txt" -o -iname "*.dat" -o -iname "*.nc" \
     -o -iname "*.res" -o -iname "*.awstxt" -o -iname "*.pdf" -o -iname "*.md" \) \
  | grep -Ei "$FNAME_PAT" || echo "  (none found)"

echo
echo "############################################"
echo "# 2) Column headers matching SEB variables"
echo "#    (tokenizes the first non-empty line)"
echo "############################################"
find "$ROOT" -type f \( -iname "*.csv" -o -iname "*.txt" -o -iname "*.dat" -o -iname "*.awstxt" -o -iname "*.res" \) -print0 \
| while IFS= read -r -d '' f; do
    matches=$(print_header_matches "$f" 2>/dev/null || true)
    if [[ -n "$matches" ]]; then
      echo "$f"
      echo "$matches"
    fi
  done

echo
echo "############################################"
echo "# 3) NetCDF variable names matching SEB terms"
echo "############################################"
if command -v ncdump >/dev/null 2>&1; then
  find "$ROOT" -type f -iname "*.nc" -print0 \
  | while IFS= read -r -d '' f; do
      vars=$(ncdump -h "$f" 2>/dev/null \
        | grep -Ei 'soil|snow|ground.?heat|skin|Tsfc|Tskin|SHF|Rnet|NETRAD|SWC|TS_' || true)
      if [ -n "$vars" ]; then
        echo "-- $f --"
        echo "$vars"
        echo
      fi
    done
else
  echo "  ncdump not found — install netcdf-bin (Linux) or 'brew install netcdf' (macOS)"
  echo "  to scan .nc variable tables. Falling back to filename-only matches above."
fi

echo
echo "Done. Cross-check any hits against the campaign's *_BIF.xlsx or"
echo "README/readme_*.txt files for exact variable definitions and depths."